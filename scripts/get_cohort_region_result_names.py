import pandas as pd 
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config-file', type = str)
    parser.add_argument('--config', type = str)
    parser.add_argument('--pheno', type = str)
    parser.add_argument('--bed-file', type = str)
    parser.add_argument('--results-path-format', type = str)
    parser.add_argument('--z-path-format', type = str)

    args = parser.parse_args()

    config_df = pd.read_csv(args.config_file, delimiter = "\s+", header = None)
    config_list = config_df.iloc[int(args.config) - 1, :].to_list()
    config_list = [x.replace("chr3", f"chr3_pheno{args.pheno}") for x in config_list]
    results_prefixes = [x.replace('chr3', f"chr3_pheno{args.pheno}") for x in config_list]
    bed_df = pd.read_csv(args.bed_file, header = None, delimiter = "\t")
    bed_df['region_fname'] = [f"chr{x[0]}.{x[1]}-{x[2]}" for x in bed_df.values]
    # write out cohort names and region names
    cohort_names = pd.DataFrame(data = config_list, columns = ['cohort'])
    cohort_names.to_csv("cohort_names.txt", header = False, index = False)
    bed_df.loc[:, "region_fname"].to_csv("region_names.txt", index = False, header = False)
    # get all the results file names
    # results-path-format will be in the form of gs://{gcloud_path}/simulations/sims/assoc_1causal/phenoXXphenoXX/XXcohortXX/XXcohortXX.XXregionXX.abf.snp
    results_name_list_of_lists = [[args.results_path_format.replace("XXphenoXX", args.pheno).replace("XXcohortXX", cohort).replace("XXregionXX", region) for cohort in config_list] for region in bed_df['region_fname']]
    results_name_list = [x for sublist in results_name_list_of_lists for x in sublist]
    results_names = pd.DataFrame(data = results_name_list, columns = ['results_fname'])
    results_names.to_csv("results_fnames.txt", header = False, index = False)
    if args.z_path_format is not None:
        z_name_list_of_lists = [[args.z_path_format.replace("XXphenoXX", args.pheno).replace("XXcohortXX", cohort).replace("XXregionXX", region) for cohort in config_list] for region in bed_df['region_fname']]
        z_name_list = [x for sublist in z_name_list_of_lists for x in sublist]
        z_names = pd.DataFrame(data = z_name_list, columns = ['z_fname'])
        z_names.to_csv("z_fnames.txt", header = False, index = False)

