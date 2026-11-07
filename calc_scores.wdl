version 1.0

import "fit_ancestry_model.wdl" as fit_ancestry
import "adjust_scores.wdl" as adjust_scores
import "subset_score_file.wdl" as subset

workflow calc_scores {
    input {
        Array[File] scorefiles
        Array[Array[File]] plink_filesets
        Boolean harmonize_scorefile = false
        Boolean add_chr_prefix = false
        Boolean ancestry_adjust = true
        File? pcs
        File? subset_variants
        String cohort = "AOU"
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

        scatter (plink_fileset in plink_filesets) {
            call plink_score_by_chromosome {
                input:
                    scorefile = scorefile_final,
                    scorefile_ncols = n_cols.ncols,
                    plink_bed = plink_fileset[0],
                    plink_bim = plink_fileset[1],
                    plink_fam = plink_fileset[2]
            }
        }

        call combine_chromosome_scores {
            input:
                chr_score_files = plink_score_by_chromosome.chr_scores_file
        }

        call combine_chromosome_variants {
            input:
                chr_variants_files = plink_score_by_chromosome.chr_variants_file
        }


        call compute_overlap {
            input:
                scorefile = select_first([chr_prefix.outfile, harmonize_score_file.scorefile_harmonized, scorefile]),
                variants = combine_chromosome_variants.variants_file
        }

        if (ancestry_adjust) {
            call fit_ancestry.find_ancestry_coefficients {
                input:
                    scores = combine_chromosome_scores.scores_file,
                    pcs = select_first([pcs, ""])
            }

            call adjust_scores.adjust_prs {
                input:
                    scores = combine_chromosome_scores.scores_file,
                    pcs = select_first([pcs, ""]),
                    mean_coef = find_ancestry_coefficients.mean_coef,
                    var_coef = find_ancestry_coefficients.var_coef
            }
        }
    }

    Array[File] scores = combine_chromosome_scores.scores_file
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


task plink_score_by_chromosome {
    input {
        File scorefile
        Int scorefile_ncols
        File plink_bed
        File plink_bim
        File plink_fam
        Int mem_gb = 16
        Int cpu = 2
    }

    Int disk_size = ceil(1.5*(size(plink_bed, "GB") + size(plink_bim, "GB") + size(plink_fam, "GB") + size(scorefile, "GB"))) + 10
    String score_cols = if (scorefile_ncols > 3) then "3-~{scorefile_ncols}" else "3"

    command <<<
        plink2 \
            --bed ~{plink_bed} \
            --bim ~{plink_bim} \
            --fam ~{plink_fam} \
            --score ~{scorefile} \
            no-mean-imputation header-read list-variants cols=-scoreavgs,+scoresums,+denom \
            --score-col-nums ~{score_cols} \
            --out tmp

        # If the output sscore file does not exist, create an empty file.
        # If it exists, this command won't change it.
        touch tmp.sscore
        touch tmp.sscore.vars

        echo "Number of variants found:" `wc -l < tmp.sscore.vars`
    >>>

    output {
        File chr_scores_file = "tmp.sscore"
        File chr_variants_file = "tmp.sscore.vars"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/rtibiocloud/plink:v2.0.0-a.6.15_18fccba"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
        cpu: "~{cpu}"
    }
}

task combine_chromosome_scores {
    input {
        Array[File] chr_score_files
    }

    command <<<

        R --vanilla << RSCRIPT

        library(tidyverse)

        score_files = str_split("~{sep=' ' chr_score_files}", " ")[[1]]
        print(score_files)

        chr_levels = c(1:22, "X", "Y")

        x = score_files %>%
            lapply(read_tsv, show_col_types=F) %>%
            bind_rows()

        # Hopefully better memory usage.
        combined = x %>%
            group_by(`#IID`) %>%
            summarise(
                ALLELE_CT=sum(ALLELE_CT),
                DENOM=sum(DENOM),
                NAMED_ALLELE_DOSAGE_SUM=sum(NAMED_ALLELE_DOSAGE_SUM),
                across(starts_with("PGS"), sum, .names = "{.col}")
            ) %>%
            mutate(
                across(starts_with("PGS"), function(x) x / DENOM, .names="{.col}_AVG")
            )

        names(combined) = str_replace(names(combined), "SUM_AVG", "AVG")
        write_tsv(combined, "scores.txt")

        RSCRIPT
    >>>

    output {
        File scores_file = "scores.txt"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/rocker/tidyverse:4.6.1"
        disks: "local-disk 10 SSD"
        memory: "16G"
    }
}

task combine_chromosome_variants {
    input {
        Array[File] chr_variants_files
    }

    command <<<

        R --vanilla << RSCRIPT

        library(tidyverse)

        vars_files = str_split("~{sep=' ' chr_variants_files}", " ")[[1]]
        print(vars_files)

        x = vars_files %>%
            lapply(read_tsv, show_col_types=F) %>%
            bind_rows()

        write_tsv(x, "variants.txt")
        RSCRIPT
    >>>

    output {
        File variants_file = "variants.txt"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/rocker/tidyverse:4.6.1"
        disks: "local-disk 10 SSD"
        memory: "16G"
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
