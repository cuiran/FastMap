import pandas as pd 
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config-file', type = str)
    parser.add_argument('--config', type = str)
    parser.add_argument('--pheno', type = str)
    parser.add_argument('--sumstats-format', type = str)
    args = parser.parse_args()

    config_df = pd.read_csv(args.config_file, delimiter = "\s+", header = None)
    config_list = config_df.iloc[int(args.config) - 1, :].to_list()
    # sumstats will be in the format gs://{gcloud_path}/simulations/sims/assoc_1causal/phenoXXpheno_numXX/XXfnameXX.PHENO.glm.linear.gz
    df = pd.DataFrame(columns = ['file_names'])
    df['file_names'] = [args.sumstats_format.replace("XXpheno_numXX", args.pheno).replace("XXfnameXX", config).replace("chr3", "chr3_pheno" + args.pheno) for config in config_list]
    df.to_csv("individual_cohorts.txt", index = False, header = False)