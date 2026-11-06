version 1.0

import "fit_ancestry_model.wdl" as fit_ancestry
import "adjust_scores.wdl" as adjust_scores
import "subset_score_file.wdl" as subset

workflow calc_scores {
    input {
        Array[File] scorefiles
        File pgen
        File pvar
        File psam
        Boolean harmonize_scorefile = false
        Boolean add_chr_prefix = false
        Boolean ancestry_adjust = true
        File? pcs
        File? subset_variants
        String cohort = "cohort"
    }

    scatter (scorefile in scorefiles) {
        if (harmonize_scorefile) {
            call harmonize_score_file {
                input:
                    scorefile = scorefile
            }
        }

        if (add_chr_prefix) {
            call chr_prefix {
                input:
                    file = select_first([harmonize_score_file.scorefile_harmonized, scorefile])
            }
        }

        if (defined(subset_variants)) {
            call subset.subset_scorefile {
                input:
                    scorefile = select_first([chr_prefix.outfile, harmonize_score_file.scorefile_harmonized, scorefile]),
                    variants = select_first([subset_variants, ""])
            }
        }

        File scorefile_final = select_first([subset_scorefile.scorefile_subset, chr_prefix.outfile, harmonize_score_file.scorefile_harmonized, scorefile])

        call n_cols {
            input:
                file = scorefile_final
        }

        call plink_score {
            input:
                scorefile = scorefile_final,
                scorefile_ncols = n_cols.ncols,
                pgen = pgen,
                pvar = pvar,
                psam = psam,
                prefix = cohort
        }

        call compute_overlap {
            input:
                scorefile = select_first([chr_prefix.outfile, harmonize_score_file.scorefile_harmonized, scorefile]),
                variants = plink_score.variants
        }

        if (ancestry_adjust) {
            call fit_ancestry.find_ancestry_coefficients {
                input:
                    scores = plink_score.scores,
                    pcs = select_first([pcs, ""])
            }

            call adjust_scores.adjust_prs {
                input:
                    scores = plink_score.scores,
                    pcs = select_first([pcs, ""]),
                    mean_coef = find_ancestry_coefficients.mean_coef,
                    var_coef = find_ancestry_coefficients.var_coef
            }
        }
    }

    Array[File] scores = plink_score.scores
    Array[File?] maybe_adjusted_scores = adjust_prs.adjusted_scores
    Array[File] overlaps = compute_overlap.overlap

    # Remove null entries before aggregating
    Array[File] adjusted_scores = select_all(maybe_adjusted_scores)

    if (defined(adjusted_scores)) {
        call aggregate_results {
            input:
                raw_scores = scores,
	            adjusted_scores = adjusted_scores,
	            overlap_files = overlaps,
                prefix = cohort
         }
    }

    output {
        File? aggregate_scores_avg = select_first([aggregate_results.aggregate_raw_avg, scores])
        File? aggregate_scores_sum = select_first([aggregate_results.aggregate_raw_sum, scores])
        File? aggregate_adjusted_scores = aggregate_results.aggregate_adjusted
        File? score_overlap = select_first([aggregate_results.aggregate_overlap, overlaps])
    }
}


task harmonize_score_file {
    input {
        File scorefile
    }

    Int disk_size = ceil(5*(size(scorefile, "GB"))) + 10
    String filename = basename(scorefile, ".gz")

    command <<<
        set -e -o pipefail
        zcat ~{scorefile} | \
        awk '{
            if (FNR==1) { print $0; next}
            split($1, a, ":")
            comp['A']='T'; comp['T']='A'; comp['C']='G'; comp['G']='C'
            if (a[4] != $2 && a[4] != comp[$2]) {
                for (i=3; i<=NF; i++) {
                    $i=$i*-1
                }
            }
            $2=a[4]
            print $0
        }' OFS="\t" > ~{filename}_harmonized
        gzip ~{filename}_harmonized
    >>>

    output {
        File scorefile_harmonized = "~{filename}_harmonized.gz"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/rtibiocloud/plink:v2.0.0-a.6.15_18fccba"
        disks: "local-disk ~{disk_size} SSD"
        memory: "16G"
    }
}


task chr_prefix {
    input {
        File file
    }

    Int disk_size = ceil(5*(size(file, "GB"))) + 10
    String filename = basename(file, ".gz")

    command <<<
        set -e -o pipefail
        zcat ~{file} | awk '{$1="chr"$1; print $0}' OFS="\t" > ~{filename}_chrprefix
        gzip ~{filename}_chrprefix
    >>>

    output {
        File outfile = "~{filename}_chrprefix.gz"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/rtibiocloud/plink:v2.0.0-a.6.15_18fccba"
        disks: "local-disk ~{disk_size} SSD"
        memory: "16G"
    }
}


task n_cols {
    input {
        File file
    }

    Int disk_size = ceil(1.5*(size(file, "GB"))) + 5

    command <<<
        Rscript -e "dat <- readr::read_tsv('~{file}', n_max=10); writeLines(as.character(ncol(dat)), 'ncols.txt')"
    >>>

    output {
        Int ncols = read_int('ncols.txt')
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/rocker/tidyverse:4.6.1"
        disks: "local-disk ~{disk_size} SSD"
        memory: "16G"
    }
}


task plink_score {
    input {
        File scorefile
        Int scorefile_ncols
        File pgen
        File pvar
        File psam
        String prefix = "out"
        Int mem_gb = 16
        Int cpu = 2
    }

    Int disk_size = ceil(1.5*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB") + size(scorefile, "GB"))) + 10
    String score_cols = if (scorefile_ncols > 3) then "3-~{scorefile_ncols}" else "3"

    command <<<
        plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} --score ~{scorefile} \
            no-mean-imputation header-read list-variants cols=+scoresums --score-col-nums ~{score_cols} \
            --out ~{prefix}
    >>>

    output {
        File scores = "~{prefix}.sscore"
        File variants = "~{prefix}.sscore.vars"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/rtibiocloud/plink:v2.0.0-a.6.15_18fccba"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
        cpu: "~{cpu}"
    }
}


task compute_overlap {
    input {
        File scorefile
        File variants
        Int mem_gb = 64
    }

    Int disk_size = ceil(10*(size(scorefile, "GB") + size(variants, "GB"))) + 10

    command <<<
        set -e
        mv ~{scorefile} scorefile.txt.gz
        gunzip scorefile.txt.gz

        R << RSCRIPT
            library(tidyverse)
            score_file <- 'scorefile.txt'
            overlap_vars <- readLines('~{variants}')
            score_vars_head <- read_tsv(score_file, n_max=10)
            id <- names(score_vars_head)[1]
            pgs <- names(score_vars_head)[str_detect(names(score_vars_head), '^PGS')]
            batches <- split(pgs, ceiling(seq_along(pgs) / 100))
            overlap <- list()
            for (b in batches) {
                cols <- c(setNames("character", id), setNames(rep("numeric", length(b)), b))
                score_vars <- data.table::fread(score_file, select=cols)
                names(score_vars)[1] <- "ID"
                for (p in b) {
                    vars <- select(score_vars, ID, weight=!!p)
                    vars <- filter(vars, weight != 0)
                    vars_overlap <- filter(vars, ID %in% overlap_vars)
                    ov <- nrow(vars_overlap) / nrow(vars)
                    wt <- sum((vars_overlap[['weight']])^2) / sum((vars[['weight']])^2)
                    overlap[[p]] <- tibble(
                    score = p,
                    n_variants = nrow(vars),
                    n_variants_overlap = nrow(vars_overlap),
                    overlap = ov,
                    beta_fraction = wt
                    )
                }
            }
            overlap <- bind_rows(overlap)
            write_tsv(overlap, 'score_overlap.tsv')
        RSCRIPT
    >>>

    output {
        File overlap = "score_overlap.tsv"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/rocker/tidyverse:4.6.1"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}


task aggregate_results {
    input {
        Array[File] raw_scores
        Array[File] adjusted_scores
        Array[File] overlap_files
        String prefix = "all"
        Int mem_gb = 8
    }

    Int disk_size = ceil(3*(size(raw_scores, "GB") + size(adjusted_scores, "GB") + size(overlap_files, "GB"))) + 10

    command <<<
        Rscript /usr/local/pgsc_calc_wdl-aou/aggregate_scores.R \
            --raw_files "~{sep=',' raw_scores}" \
            --adjusted_files "~{sep=',' adjusted_scores}" \
            --overlap_files "~{sep=',' overlap_files}" \
            --out_prefix ~{prefix}
    >>>

    output {
        File aggregate_raw_avg = "~{prefix}_raw_avg_scores.tsv"
        File aggregate_raw_sum = "~{prefix}_raw_sum_scores.tsv"
        File aggregate_adjusted = "~{prefix}_adjusted_scores.tsv"
        File aggregate_overlap = "~{prefix}_overlap.tsv"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/uwgac/calc_scores_aou:0.0.1"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb} GB"
	}
}
