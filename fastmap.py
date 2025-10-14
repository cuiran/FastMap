#! /usr/bin/python3

import pandas as pd
import numpy as np
import argparse
import scipy.stats
from functools import reduce
import warnings
from scipy.stats import chi2
import rpy2.robjects as ro
warnings.filterwarnings("ignore")

def outer_merge_on_index(df1, df2):
    return pd.merge(df1, df2, how='outer', left_index=True, right_index=True)

def share_cohort_check(row_name, col_name):
    # check if the combined names have shared cohort names
    # If there are shared, return True, if not return False
    # The purpose of this is to not combine different alpha columns from the same cohort
    # row_name and colname are of the form EUR_sim4_chr3_Omni25_1000G_alpha1 or EUR_sim4_chr3_Omni25_1000G_alpha1,EUR_sim3_chr3_Omni25_HRC_alpha1
    row_cohorts = ['_'.join(x.split('_')[:-1]) for x in row_name.split(',')]
    col_cohorts = ['_'.join(y.split('_')[:-1]) for y in col_name.split(',')]
    shared = set(row_cohorts)&set(col_cohorts)
    return len(shared) > 0

def compute_coverage_and_jaccard(df, reverse_colname_map = dict(), prev_pip_coverage_df = pd.DataFrame(),
                                prev_jaccard_df = pd.DataFrame(), top_n = 100):
    names = df.columns
    jaccard_df = pd.DataFrame(index = names, columns = names)
    pip_coverage_df = pd.DataFrame(index = names, columns = names)
    for i,row in enumerate(names):
        for j,col in enumerate(names):
            if j <= i:
                continue
            if share_cohort_check(reverse_colname_map[row], reverse_colname_map[col]):
                # if the alphas are from the same cohort, skip it
                continue
            if row in prev_pip_coverage_df.index and col in prev_pip_coverage_df.columns:
                jaccard_df.loc[row, col] = prev_jaccard_df.loc[row, col]
                pip_coverage_df.loc[row, col] = prev_pip_coverage_df.loc[row, col]
                pip_coverage_df.loc[col, row] = prev_pip_coverage_df.loc[col, row]
                continue
            sel_no_na = df[[row, col]].dropna().copy()
            sel_no_na[f"rank_{row}"] = sel_no_na[row].rank(ascending=False, method='first').astype(int)
            sel_no_na[f"rank_{col}"] = sel_no_na[col].rank(ascending=False, method='first').astype(int)
            # Get top-N variant indices for each ancestry
            top_row = set(sel_no_na[sel_no_na[f"rank_{row}"] <= top_n].index)
            top_col = set(sel_no_na[sel_no_na[f"rank_{col}"] <= top_n].index)
            union_idx = list(top_row | top_col)
            intersect_idx = list(top_row & top_col)
            # Compute Jaccard index, allowing non-overlapping sets
            if not union_idx:
                jaccard_df.loc[row, col] = 0
                continue
            jaccard = len(intersect_idx) / len(union_idx) if union_idx else 0
            jaccard_df.loc[row, col] = jaccard
            # Compute PIP coverages
            a, b = sel_no_na.loc[intersect_idx, [row, col]].to_numpy().T
            pip_coverage_df.loc[row, col] = a.sum()
            pip_coverage_df.loc[col, row] = b.sum()
    return pip_coverage_df, jaccard_df

def compute_scaling_factors(df):
    df_remove_na_rows = df.dropna(how = 'all').copy()
    # scale PIPs based on how many missing variants there are
    scaling_factors = (len(df_remove_na_rows) - df_remove_na_rows.isnull().sum())/len(df_remove_na_rows)
    return scaling_factors

def combine_pips(df):
    df_remove_na_rows = df.dropna(how = 'all').copy()
    # scale PIPs based on how many missing variants there are
    scaling_factors = (len(df_remove_na_rows) - df_remove_na_rows.isnull().sum())/len(df_remove_na_rows)
    df_scaled = df_remove_na_rows*scaling_factors
    print(f'Scaling {df.columns.to_list()} by ')
    print(f'{scaling_factors}')
    df_scaled.fillna(1/len(df_scaled), inplace = True)
    result_df = pd.DataFrame(index = df_scaled.index)
    product = df_scaled.prod(axis = 1)
    result_df['prob_fastmap'] = product / product.sum()
    result_df = result_df.reindex(df.index)
    return result_df, scaling_factors

def get_components(sorted_priority_score, jaccard_df, jaccard_thresh = 0.5, L = 10):
    used_keys = []
    sorted_key_list = list(sorted_priority_score.keys())
    used_keys.append(sorted_key_list[0])
    counter = 1
    while counter < L:
        used_df = jaccard_df[used_keys]
        jaccard_similar_idx = list(used_df[(used_df>jaccard_thresh).any(axis = 1)].index)
        use_key = next((key for key in sorted_key_list if key not in used_keys+jaccard_similar_idx))
        used_keys.append(use_key)
        counter += 1
    return used_keys

def get_unused_cohorts(used_keys, cohorts, reverse_colname_map):
    used_cohorts = set(['_'.join(x.split('_')[:-1]) for sublist in [reverse_colname_map[k].split(',') for k in used_keys] for x in sublist])
    print(f"All cohorts: {cohorts}")
    print(f"Used cohorts: {used_cohorts}")
    return set(cohorts).difference(used_cohorts)

def get_fname(region, cohort, results_files):
    l = [x for x in results_files if region in x and cohort in x]
    if len(l) == 1:
        return l[0]
    else:
        raise ValueError(f"No file name contain {region} and {cohort} in the list of results files provided.")

def ivw_meta_analysis(betas, ses):
    # Remove pairs where either beta or se is NaN
    mask = ~np.isnan(betas) & ~np.isnan(ses)
    betas = betas[mask]
    ses = ses[mask]
    
    if len(betas) == 0:
        return np.nan, np.nan, np.nan

    weights = 1 / ses**2
    meta_beta = np.sum(weights * betas) / np.sum(weights)
    meta_se = np.sqrt(1 / np.sum(weights))
    chi2_stat = (meta_beta / meta_se)**2
    p_value = chi2.sf(chi2_stat, df=1)
    return meta_beta, meta_se, p_value

def apply_ivw_to_dataframe(df, variant_colname):
    beta_cols = [x for x in df.columns if "beta" in x]
    se_cols = [x for x in df.columns if "se" in x]
    
    results = df.apply(lambda row: ivw_meta_analysis(
        row[beta_cols].values.astype(float),
        row[se_cols].values.astype(float)
    ), axis=1)
    
    # Convert results into separate columns
    meta_df = pd.DataFrame(results.tolist(), columns=['meta_beta', 'meta_se', 'meta_pval'])
    return pd.concat([df[[variant_colname]], meta_df], axis=1)

def get_all_results(regions = [], cohorts = [],
                    wdl = False, dirname_format = "", pheno = "",
                    results_files = [], z_files = [],
                    variant_colname = "variant",
                    beta_colname = "beta",
                    se_colname = "se",
                    prob_colname = "prob",
                    susie_L = 10):
    # store all alpha results in a dictionary called results
    # with key = region, value = dataframe for that region each L columns store the alpha values for a cohort
    # store all marginal pip values in a dictionary called pips with key = region, value = dataframe for that region each column store the pips for a cohort
    results = dict()
    pips = dict()
    info_dict = dict()
    alpha_cols = [f'alpha{i}' for i in range(1, susie_L+1)]
    for region in regions:
        results_list = []
        pip_list = []
        info_list = []
        for cohort in cohorts:
            fname = get_fname(region, cohort, results_files)
            z_fname = get_fname(region, cohort, z_files)
            print(f"Reading in results from cohort {cohort}, region {region}, file name {fname}...")
            r_data = ro.r['readRDS'](fname)
            z_df = pd.read_csv(z_fname, usecols = [variant_colname, beta_colname, se_colname], delimiter = '\s+')
            alpha_mat = np.asarray(r_data.rx2('alpha'))
            pip = np.asarray(r_data.rx2("pip"))
            prefix_without_pheno = cohort.replace(f'chr3_pheno{pheno}', 'chr3')
            # add cohort names (without the phenotype numbers) to column names
            cohort_specific_alpha_cols = [f"{prefix_without_pheno}_{x}" for x in alpha_cols]
            cohort_specific_prob_col = f"{prefix_without_pheno}_{prob_colname}"
            cohort_specific_beta_se_cols = [f"{prefix_without_pheno}_{x}" for x in [beta_colname, se_colname]]
            alpha_df = pd.DataFrame(index = z_df[variant_colname], columns = cohort_specific_alpha_cols)
            pip_df = pd.DataFrame(index = z_df[variant_colname], columns = [cohort_specific_prob_col])
            info_df = pd.DataFrame(index = z_df[variant_colname], columns = cohort_specific_beta_se_cols)
            alpha_df[cohort_specific_alpha_cols] = alpha_mat.T
            pip_df[cohort_specific_prob_col] = pip
            info_df[cohort_specific_beta_se_cols] = z_df[[beta_colname, se_colname]].values
            # store SuSiE alpha results and marginal PIP results in lists
            results_list.append(alpha_df)
            pip_list.append(pip_df)
            info_list.append(info_df)
        # merge all cohorts results in one region
        region_results = reduce(outer_merge_on_index, results_list)
        results[region] = region_results
        region_pips = reduce(outer_merge_on_index, pip_list)
        pips[region] = region_pips
        region_info = reduce(outer_merge_on_index, info_list)
        info_dict[region] = region_info.reset_index()
    return results, pips, info_dict, variant_colname

def fastmap_multi_causal(results_dict, pips_dict, info_dict, variant_colname = 'variant', max_thresh = 0.1, initial_rev_max_thresh = 0.001, 
                         max_div = 3, pip_thresh = 0.5, jaccard_thresh = 0, L=10, top_n = 100,
                         meta_sig_thresh = 5e-8, nonsig_pip_thresh = 0.1):
    fastmap_results_list = []
    for region in results_dict.keys():
        print(f'--------------------------- Region {region} ---------------------------')
        region_df = results_dict[region]
        colname_map = dict({colname: ind for ind, colname in enumerate(region_df.columns)}) # change column names to numbers, easier to parse
        reverse_colname_map = {colname_map[k]: k for k in colname_map.keys()}
        updated_region_df = region_df.copy().rename(columns = colname_map)
        all_cols_df = updated_region_df.copy() # this is to track all generated columns
        cohorts = list(set(['_'.join(x.split('_')[:-1]) for x in region_df.columns]))
        # compute initial pip coverage and Jaccard matrices
        initial_pip_coverage_df, initial_jaccard_df = compute_coverage_and_jaccard(updated_region_df, reverse_colname_map = reverse_colname_map, top_n=top_n)
        # each ancestry result gets a priority score based on their max PIP
        # the final result uses the priority score to determine which columns to prioritize
        priority_score = {x: all_cols_df[x].max() for x in all_cols_df.columns}
        
        # some initial settings
        div = 1 # divisor for the thresholds
        rev_max_thresh = initial_rev_max_thresh
        used_keys = []
        combined = False
        while div <= max_div:
            pip_coverage_df, jaccard_df = (initial_pip_coverage_df.copy(), initial_jaccard_df.copy())
            max_pip_coverage = pip_coverage_df.unstack().max()
            if max_pip_coverage < max_thresh/max_div: 
                print(f'max PIP coverage {max_pip_coverage} is too low')
                break
            if div > 1: # if div is already increased to 2, we need to find a rev_max_thresh that works for this data
                max_col_name, max_row_name = pip_coverage_df.infer_objects(copy = False).fillna(-2).unstack().idxmax()
                rev_max_thresh = pip_coverage_df.loc[max_col_name, max_row_name]/div
                print(f'############# Lowering thresholds to 1/{div} of previous #############')
            print(f'############# current max PIP threshold: {max_thresh/div}, reverse max threshold: {rev_max_thresh} #############')
            while max_pip_coverage >= max_thresh/div:
                max_col_name, max_row_name = pip_coverage_df.infer_objects(copy = False).fillna(-2).unstack().idxmax()
                print(max_col_name, max_row_name)
                print(f"Max PIP coverage {max_pip_coverage} occurs between {reverse_colname_map[max_row_name]} and {reverse_colname_map[max_col_name]}")
                print(f"Jaccard col,row={jaccard_df.loc[max_col_name, max_row_name]}, row,col = {jaccard_df.loc[max_row_name, max_col_name]}, reverse pip coverage = {pip_coverage_df.loc[max_col_name, max_row_name]}")
                new_name = f"{max_row_name},{max_col_name}"
                # A few conditions need to be satistifed before combining
                jaccard_condition = ((jaccard_df.loc[max_col_name, max_row_name] > jaccard_thresh) | (jaccard_df.loc[max_row_name, max_col_name] > jaccard_thresh))
                if not jaccard_condition: print("Jaccard condition not satistifed")
                rev_pip_coverage_condition = (pip_coverage_df.loc[max_col_name, max_row_name] >= rev_max_thresh)
                if not rev_pip_coverage_condition: print(f"Reverse PIP coverage below threshold {rev_max_thresh}")
                no_repeat_condition = (new_name not in all_cols_df.columns)
                if not no_repeat_condition: print(f"{new_name} already exists")

                if jaccard_condition & rev_pip_coverage_condition & no_repeat_condition:
                    selected_cols = [max_row_name, max_col_name]
                    print(f"combining {max_row_name} and {max_col_name}")
                    combined_pips_df,_ = combine_pips(updated_region_df.loc[:, selected_cols].copy())
                    combined = True
                    assert len(combined_pips_df) == len(updated_region_df)
                    assert all(combined_pips_df.index == updated_region_df.index)
                    # drop combined ancestries out of dataframe
                    updated_region_df.drop([max_row_name, max_col_name], axis = 1, inplace = True)
                    # assign the combined pips to a new column
                    reverse_colname_map.update({new_name: f"{reverse_colname_map[max_row_name]},{reverse_colname_map[max_col_name]}"})
                    updated_region_df[new_name] = combined_pips_df['prob_fastmap']
                    all_cols_df[new_name] = combined_pips_df['prob_fastmap']
                    # record priority score
                    # bonux 1 point for combination of ancestries
                    priority_score[new_name] = len(new_name.split(',')) + max_pip_coverage + np.nanmax([jaccard_df.loc[max_col_name, max_row_name], jaccard_df.loc[max_row_name, max_col_name]])
                    # compute new dfs
                    pip_coverage_df, jaccard_df = compute_coverage_and_jaccard(updated_region_df, reverse_colname_map = reverse_colname_map, 
                                                                               prev_pip_coverage_df = pip_coverage_df.copy(),
                                                                                   prev_jaccard_df = jaccard_df.copy(), top_n = top_n)
                    max_pip_coverage = pip_coverage_df.unstack().max()
                    print(updated_region_df[[new_name]].sort_values(new_name, ascending = False).head(3))
                    # after obtaining the final updated_region_df, 
                    # use the column with the highest priority score and that are the least similar as components for the marginal PIP computation
                    sorted_priority_score = {k: v for k, v in sorted(priority_score.items(), key=lambda item: item[1], reverse = True) if k in updated_region_df.columns}
                    used_keys = get_components(sorted_priority_score, jaccard_df = jaccard_df, jaccard_thresh=0.1, L=L)
                else:        
                    # assign -1 value to the location and look again
                    pip_coverage_df.loc[max_row_name, max_col_name] = -1
                    max_pip_coverage = pip_coverage_df.unstack().max()
            # check if combination has occured and if all signals have combined
            # get high-PIP columns
            high_pip_columns = updated_region_df.loc[:,(updated_region_df > pip_thresh).any()].columns
            mapped_original_cols = [colname_map[x] for x in region_df.columns]
            signal_not_combined = any([x in mapped_original_cols for x in high_pip_columns]) # check if any high pip columns are still left un-combined
            if (not combined) or signal_not_combined:
                div = div + 1
                updated_region_df = region_df.copy().rename(columns = colname_map)
            else:
                break
        if len(used_keys) == 0:
            sorted_priority_score = {k: v for k, v in sorted(priority_score.items(), key=lambda item: item[1], reverse = True) if k in updated_region_df.columns}
            used_keys = get_components(sorted_priority_score, jaccard_df = jaccard_df, jaccard_thresh=0.1, L=L)
        
        # perform post-hoc QC
        # find variants with moderate-to-high PIPs but do not pass GWAS significance in meta-analysis
        # assign such variants PIP=1/M, M being the number of SNPs in the region.
        '''
        print('Performing post-hoc QC using meta-analysis p-values')
        meta_df = apply_ivw_to_dataframe(info_dict[region], variant_colname)
        meta_sig_df = meta_df[meta_df['meta_pval']<meta_sig_thresh]
        meta_sig_variants = meta_sig_df[variant_colname].to_list()
        row_mask = ~all_cols_df.index.isin(meta_sig_variants)
        row_mask_df = pd.DataFrame(np.tile(row_mask[:, np.newaxis], (1, all_cols_df.shape[1])), 
                                index=all_cols_df.index, columns=all_cols_df.columns)
        mask = row_mask_df & (all_cols_df > nonsig_pip_thresh)
        unif_prior = 1/len(all_cols_df)
        all_cols_df[mask] = unif_prior
        print(f"############# Removed {len(np.where(mask)[0])} variants PIP >{nonsig_pip_thresh} and not meta significant with threshold {meta_sig_thresh} #############")
        '''
        
        # compute combined PIPs for variants that have at least one non-NaN value in the selected columns 
        df_no_na = all_cols_df[used_keys].dropna(how = 'all')
        non_nan_part = pd.DataFrame(index = df_no_na.index)
        print(f"Computing PIP from {len(used_keys)} components: {used_keys}, corresponding to {[reverse_colname_map[x] for x in used_keys]}")
        non_nan_part['prob_fastmap'] = 1 - np.prod(1 - df_no_na[used_keys], axis = 1) # prod ignores NaN, i.e. treats them as 1 in the multiplication
        non_nan_part['note'] = 'from_combining'
        # fill in NA with MAX of the un-used cohorts
        if len(non_nan_part) == len(all_cols_df):
            fastmap_df = non_nan_part
        else:
            df_all_na = all_cols_df[used_keys][all_cols_df[used_keys].isna().all(axis=1)]
            all_na_part = pd.DataFrame(index = df_all_na.index)
            unused_cohorts = get_unused_cohorts(used_keys, cohorts, reverse_colname_map)
            print(f"Unused cohorts: {unused_cohorts}")
            print(pips_dict[region].columns)
            print(unused_cohorts)
            all_na_part['prob_fastmap'] = pips_dict[region][[f"{x}_prob" for x in unused_cohorts]].max(axis = 1).reindex(all_na_part.index).values
            all_na_part['note'] = 'from_filling_max'
            fastmap_df = pd.concat([non_nan_part, all_na_part], axis = 0)
            assert len(fastmap_df) == len(all_cols_df)
            assert len(fastmap_df[fastmap_df['prob_fastmap'].isna()]) == 0
        fastmap_df['region'] = region
        fastmap_results_list.append(fastmap_df)
    fastmap_all_regions = pd.concat(fastmap_results_list, axis = 0)
    return fastmap_all_regions

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort-fnames', type = str)
    parser.add_argument('--region-names', type = str)
    parser.add_argument('--config-name', type = str)
    parser.add_argument('--config', type = str)
    parser.add_argument('--pheno', type = str)
    parser.add_argument('--results-files', type = str, help = "comma separated file names")
    parser.add_argument('--z-files', type = str, help = "comma separated file names")
    parser.add_argument('--wdl', action = 'store_true')
    parser.add_argument('--dirname-format', type = str)
    parser.add_argument('--variant-colname', type = str)
    parser.add_argument('--prob-colname', type = str)
    parser.add_argument('--susie-l', type = int, help = "Number of components output by SuSiE. Usually 10.")
    parser.add_argument('--max-threshold', type = float, help = "Lowest threshold for maximum PIP coverage. Lowering this to less than 0.2 may lower the quality of the results.")
    parser.add_argument('--rev-max-threshold', type = float, help = "Set the minimum reverse pip coverage (if the max PIP coverage occurs at row 3 col 4, the reverse PIP coverage is row 4 col 3). Lowering this to less than 0.01 may lower the quality of the results.")
    parser.add_argument('--jaccard-threshold', type = float, help = "Set the threshold for the Jaccard similarity index, allowing combining if two components Jaccard similarity is above this threshold.")    
    parser.add_argument('--top-n', type = int, help = "The top n variants to use in order to check for similarity and PIP coverages.")
    parser.add_argument('--max-div', type = int, help = "Maximum divisor allowed for the thresholds to divide by. Increasing this allows for more combinations of PIPs but also may lower the quality of the results.")
    parser.add_argument('--pip-thresh', type = float, help = "Lower max_threshold to encourage all components with PIP higher than this threshold to combine.")
    parser.add_argument('--L', type = int, help = "Maximum allowed number of causal SNPs for FastMap")
    args = parser.parse_args()

    for key, value in vars(args).items():
        print(f"{key}: {value}")

    region_names_df = pd.read_csv(args.region_names, header = None)
    region_list = region_names_df[0].values.flatten()
    cohort_names_df = pd.read_csv(args.cohort_fnames, header = None)
    cohort_list = cohort_names_df[0].values.flatten()

    results, pips, info_dict, variant_colname = get_all_results(regions = region_list,
                                    cohorts = cohort_list,
                                    wdl = args.wdl,
                                    dirname_format = args.dirname_format,
                                    pheno = args.pheno,
                                    results_files = args.results_files.split(","),
                                    z_files = args.z_files.split(","),
                                    variant_colname = args.variant_colname,
                                    prob_colname = args.prob_colname,
                                    susie_L = args.susie_l)
    
    fastmap_result_df = fastmap_multi_causal(results, pips, info_dict, variant_colname= variant_colname,
                                             max_thresh = args.max_threshold, 
                                             initial_rev_max_thresh = args.rev_max_threshold, 
                                             jaccard_thresh = args.jaccard_threshold,
                                             max_div = args.max_div, L = args.L, top_n = args.top_n,
                                             pip_thresh = args.pip_thresh)
    
    out_prefix = f"{args.config_name}.config{args.config}.pheno{args.pheno}"
    print(f"Writing FastMap results to files {out_prefix}.fastmap.snp.feather")
    fastmap_result_df.reset_index().to_feather(f"{out_prefix}.fastmap.snp.feather")

