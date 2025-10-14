version 1.0
task get_individual_cohorts {
    input {
        String docker
        String zones

        Pair[String, String] config_pheno
        String config = config_pheno.left
        String pheno = config_pheno.right
        File config_file
        File script
        String sumstats_format
        String config_name
        String meta_bedfile_format
    }

    command <<<
        python3 ~{script} \
            --config ~{config} \
            --pheno ~{pheno} \
            --config-file ~{config_file} \
            --sumstats-format ~{sumstats_format}
    >>>

    output {
        Array[String] individual_cohorts = read_lines("individual_cohorts.txt")
        String bed_file = sub(sub(sub(meta_bedfile_format, "XXconfigXX", config), "XXphenoXX", pheno), "XXconfig_nameXX", config_name)
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "1 GB"
        disks: "local-disk 5 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task get_cohort_region_result_names {
    input {
        String zones
        String docker
        Pair[String, String] config_pheno
        String config = config_pheno.left
        String pheno = config_pheno.right
        File config_file
        File bed_file
        File script
        String results_path_format
        String z_path_format
    }

    command <<<
        python3 ~{script} \
            --config-file ~{config_file} \
            --config ~{config} \
            --pheno ~{pheno} \
            --bed-file ~{bed_file} \
            --results-path-format ~{results_path_format} \
            --z-path-format ~{z_path_format}
    >>>

    output {
        File cohort_names = "cohort_names.txt"
        File region_names = "region_names.txt"
        File results_fnames = "results_fnames.txt"
        Array[String] results_files = read_lines("results_fnames.txt")
        Array[String] z_files = read_lines("z_fnames.txt")
    }
    runtime {
        docker: "${docker}"
        zones: "${zones}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 5 HDD"
        preemptible: 2
        noAddress: false
    }
}

task fastmap_multi_causal {
    input {
        String zones
        String docker
        File cohort_names
        File region_names
        Array[File] results_files
        Array[File] z_files
        String config_name
        String config
        String pheno
        String variant_colname
        Int susie_l
        Int L
        Float max_threshold
        Float rev_max_threshold
        Float jaccard_threshold
        Float pip_threshold
        Int max_div
        Int top_n
        String folder_name
        String sim_folder_name
        String prob_colname
        File script
        Int cpu
        Int mem
    }

    command <<<
        python3 ~{script} \
            --cohort-fnames ~{cohort_names} \
            --region-names ~{region_names} \
            --config-name ~{config_name} \
            --config ~{config} \
            --pheno ~{pheno} \
            --variant-colname ~{variant_colname} \
            --results-files ~{sep=',' results_files} \
            --z-files ~{sep=',' z_files} \
            --max-threshold ~{max_threshold} \
            --rev-max-threshold ~{rev_max_threshold} \
            --jaccard-threshold ~{jaccard_threshold} \
            --max-div ~{max_div} \
            --L ~{L} \
            --susie-l ~{susie_l} \
            --top-n ~{top_n} \
            --pip-thresh ~{pip_threshold} \
            --prob-colname ~{prob_colname}

        gsutil -m cp ~{config_name}.config~{config}.pheno~{pheno}*.snp.feather gs://{gcloud_path}/sims/~{sim_folder_name}/fastmap/~{folder_name}/~{config_name}/
        cp stdout ~{config_name}.config~{config}.pheno~{pheno}.log
        gsutil -m cp ~{config_name}.config~{config}.pheno~{pheno}.log gs://{gcloud_path}/sims/~{sim_folder_name}/fastmap/~{folder_name}/~{config_name}/
    >>>

    output {
        File fastmap_result = config_name + '.config' + config + '.pheno' + pheno + '.fastmap.snp.feather'
    }
    runtime {
        docker: "${docker}"
        zones: "${zones}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 100 HDD"
        preemptible: 2
        noAddress: false
    }
}

workflow fastmap_meta_regions {
    input {
        String zones
        String docker

        File configlistfile
        File phenolistfile
        File config_file
        String config_name = basename(config_file, ".tsv")

        Array[String] configs = read_lines(configlistfile)
        Array[String] phenos = read_lines(phenolistfile)

        Array[Pair[String, String]] all_pairs = cross(configs, phenos)
    }

    meta {
        allowNestedInputs: true
    }

    scatter (config_pheno in all_pairs) {
        call get_individual_cohorts {
            input: zones = zones,
                docker = docker,
                config_pheno = config_pheno,
                config_file = config_file,
                config_name = config_name
        }

        call get_cohort_region_result_names {
            input: zones = zones,
                docker = docker,
                config_pheno = config_pheno,
                config_file = config_file,
                bed_file = get_individual_cohorts.bed_file
        }

        call fastmap_multi_causal {
            input: zones = zones,
                docker = docker,
                cohort_names = get_cohort_region_result_names.cohort_names,
                region_names = get_cohort_region_result_names.region_names,
                results_files = get_cohort_region_result_names.results_files,
                z_files = get_cohort_region_result_names.z_files,
                config_name = config_name,
                config = config_pheno.left,
                pheno = config_pheno.right
        }
    }
}