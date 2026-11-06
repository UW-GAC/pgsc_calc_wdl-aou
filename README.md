# pgsc_calc_wdl

## pgsc_calc

WDL wrapper for the [pgsc_calc](https://pgsc-calc.readthedocs.io/en/latest/) workflow

The first step of the workflow [formats the input genomes](https://pgsc-calc.readthedocs.io/en/latest/how-to/prepare.html) in PLINK2 format.

The next step [runs the pgsc_calc nextflow workflow](https://pgsc-calc.readthedocs.io/en/latest/getting-started.html) inside a docker container.

Either pgs_id or scorefile should be specified, not both. If a scorefile is provided, its genome build should match the genotypes and be specified in `target_build`.

input | description
--- | ---
vcf | Array of VCF files. If provided, will be converted to pgen/pvar/psam. If not provided, use pgen/pvar/psam inputs instead.
pgen | Array of pgen files
pvar | Array of pvar files
psam | Array of psam files
chromosome | Array of chromosome strings (1-22, X, Y) corresponding to `vcf` or `pgen/pvar/psam`. If there is one file with multiple chromosomes, this input should be an empty string (`[""]`)
target_build | `"GRCh38"` (default) or `"GRCh37"`
pgs_id | PGS catalog IDs to calculate (e.g. `"PGS001229, PGS000802"`)
scorefile | Score file in the [required format](https://pgsc-calc.readthedocs.io/en/latest/how-to/calculate_custom.html).
ancestry_ref_panel | Google bucket path of a reference panel file (e.g. `"gs://fc-a8511200-791a-4375-bccf-fbe41ac3f9f6/pgsc_HGDP+1kGP_v1.tar.zst"`) to [perform ancestry adjustment](https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html). If not provided, no ancestry adjustment is performed.
sampleset_name | Name of the sampleset; used to construct output file names (default `"cohort"`)
arguments | [Additional arguments](https://pgsc-calc.readthedocs.io/en/latest/reference/params.html#param-ref) to pass to psgc_calc

Output files from pgsc_calc are described [here](https://pgsc-calc.readthedocs.io/en/latest/explanation/output.html#interpret).


## pgsc_calc_prepare_genomes

Standalone workflow to convert VCF to pgen/pvar/psam. The pvar generated from the workflow will have variant ids in the form chr:pos:ref:alt without the "chr" prefix.

input | description
--- | ---
vcf | Array of VCF files
merge_chroms | Boolean for whether to merge files to a single set of output files with all chromosomes
snps_only | Boolean for whether to keep only SNPs in output

output description
--- | ---
pgen | Array of pgen files
pvar | Array of pvar files
psam | Array of psam files


## calc_scores

Calculate scores without using Nextflow. Use pgsc_calc_prepare_genomes first to generate files.

The first two columns of the score file are expected to be "ID" (variant id) and "effect_allele". All subsequent columns are score weights, with the header of each column containing the score name.

If ancestry_adjust is true, a file with PCs must be supplied. The PC file is expected to have column "IID" for the sample ID and columns starting with "PC" for the PCs. For each score in the scorefile, the code fits a model regressing the score on the PCs. The coefficients from this model are used to adjust the score. Both unadjusted and adjusted scores are returned by the workflow.

If your pvar file contains variant ids starting with a "chr" prefix but your score file does not, set add_chr_prefix to "true" to add the "chr" prefix to the score file to match the pvar. pgsc_calc_prepare_genomes generates a pvar file without the "chr" prefix, so if you prepared your genomes with that workflow, "add_chr_prefix" should be set to "false" (default).

input | description
--- | ---
scorefile | score file
pgen | pgen file
pvar | pvar file
psam | psam file
harmonize_scorefile | Boolean for whether to harmonize scorefile to consistent effect allele (default false)
add_chr_prefix | Boolean for whether to add "chr" prefix to scorefile variant ids to match pvar (default false)
ancestry_adjust | Boolean for whether to adjust scores for ancestry using PCs
pcs | optional file with PCs to adjust for ancestry
subset_variants | optional list of variants to subset from the scorefile before running

output description
--- | ---
scores | sscore file output by plink
adjusted_scores | if ancestry_adjust is true, a file with adjusted scores is returned
variants | variants included in sscore
overlap | TSV file with fraction of overlapping variants and ratio of squared weights for each score


## combine_weighted_scores

Using harmonized scoring files (same effect allele for each variant) and a file of weights for a set of scores, weight each score and sum over the weighted scores.

input | description
--- | ---
score_files | Array of score files
weights | File with columns "score" and "weight" specifying the weight for each score. The scores should be a subset of those in score_files.

output | description
--- | ---
combined_scores | File with combined weighted scores for each variant. Will have three columns: "ID" (variant id), "effect_allele", and "score".

# Building the AOU-specific docker image

This docker image includes all the files that were originally obtained with wget or `source` within a command block, which is not allowed on the AOU platform.

1. Push all changes to the repository. Note that the Docker image will build off the "main" branch on GitHub.

1. Build the image. Make sure to include no caching, or else local scripts will not be updated.

    ```bash
    docker build -f calc_scores_aou.dockerfile --no-cache -t uwgac/calc_scores_aou:X.Y.Z .
    ```

1. Push the image to Docker Hub.

    ```bash
    docker push uwgac/calc_scores_aou:X.Y.Z
    ```
