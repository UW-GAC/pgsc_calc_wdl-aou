version 1.0

workflow adjust_scores {
    input {
        File scores
        File pcs
        File mean_coef
        File var_coef
    }

    call adjust_prs {
        input:
            scores = scores,
            pcs = pcs,
            mean_coef = mean_coef,
            var_coef = var_coef
    }

    output {
        File adjusted_scores = adjust_prs.adjusted_scores
    }
}


task adjust_prs {
    input {
        File scores
        File pcs
        File mean_coef
        File var_coef
        Int mem_gb = 16
    }

    Int disk_size = ceil(2.5*(size(scores, "GB") + size(pcs, "GB") + size(mean_coef, "GB") + size(var_coef, "GB"))) + 10

    command <<<
        Rscript -e "
        library(readr); \
        source('/usr/local/pgsc_calc_wdl-aou/ancestry_adjustment.R'); \
        scores <- read_tsv('~{scores}'); \
        pcs <- read_tsv('~{pcs}'); \
        mean_coef <- read_tsv('~{mean_coef}'); \
        var_coef <- read_tsv('~{var_coef}'); \
        scores <- prep_scores(scores); \
        adjusted_scores <- adjust_prs(scores, pcs, mean_coef, var_coef); \
        write_tsv(adjusted_scores, 'adjusted_scores.txt'); \
        "
    >>>

    output {
        File adjusted_scores = "adjusted_scores.txt"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/uwgac/calc_scores_aou:0.0.3"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}
