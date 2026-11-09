version 1.0

workflow fit_ancestry_model {
    input {
        File scores
        File pcs
    }

    call find_ancestry_coefficients {
        input:
            scores = scores,
            pcs = pcs
    }

    output {
        File mean_coef = find_ancestry_coefficients.mean_coef
        File var_coef = find_ancestry_coefficients.var_coef
    }
}


task find_ancestry_coefficients {
    input {
        File scores
        File pcs
        Int mem_gb = 16
    }

    Int disk_size = ceil(2.5*(size(scores, "GB") + size(pcs, "GB"))) + 10

    command <<<
        Rscript -e "
        library(readr); \
        source('/usr/local/pgsc_calc_wdl-aou/ancestry_adjustment.R'); \
        scores <- read_tsv('~{scores}'); \
        pcs <- read_tsv('~{pcs}'); \
        scores <- prep_scores(scores); \
        model <- fit_prs(scores, pcs); \
        write_tsv(model[['mean_coef']], 'mean_coef.txt'); \
        write_tsv(model[['var_coef']], 'var_coef.txt'); \
        "
    >>>

    output {
        File mean_coef = "mean_coef.txt"
        File var_coef = "var_coef.txt"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/uwgac/calc_scores_aou:0.0.1"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}
