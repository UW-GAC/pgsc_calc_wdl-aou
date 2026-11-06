version 1.0

workflow subset_score_file {
    input {
        File scorefile
        File variants
    }

    call subset_scorefile {
        input:
            scorefile = scorefile,
            variants = variants
    }

    output {
        File scorefile_subset = subset_scorefile.scorefile_subset
    }

}


task subset_scorefile {
    input {
        File scorefile
        File variants
        Int mem_gb = 64
    }

    Int disk_size = ceil(5*(size(scorefile, "GB") + size(variants, "GB"))) + 10
    String filename = basename(scorefile, ".gz")

    command <<<
        Rscript -e " \
        library(tidyverse); \
        score_vars <- read_tsv('~{scorefile}'); \
        overlap_vars <- readLines('~{variants}'); \
        names(score_vars)[1] <- 'ID'; \
        score_vars <- filter(score_vars, is.element(ID, overlap_vars)); \
        write_tsv(score_vars, '~{filename}_subset.gz'); \
        "
    >>>

    output {
        File scorefile_subset = "~{filename}_subset.gz"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/rocker/tidyverse:4.6.1"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}
