## Path to your read counts file
## Format: transcript_name<tab>sample1<tab>sample2<tab>sample3<tab>sample4<tab>etc.
## Example: GeneA<tab>100<tab>200<tab>300<tab>400<tab>etc.
## Header in the first line is required with the sample names
counts = "counts.txt"

## Path to your samples file containing the condition and the sample name
## Format: sample_name<tab>condition<tab>sample_name
## Example: sample1<tab>control<tab>sample1
## Must label the condition column as "condition"; for RIP analysis, must inclue the "assay" column for Input or IP
samples = "samples.txt"

## Type of test to perform in DESeq2; Wald test, Likelihood ratio test (LRT), or IP/Input enrichment test (RIP)
analysis_type = "LRT"

rule all:
    input:
        "sizeFactors.txt",
        # "clusters.txt"
## Run DESeq2
rule run_deseq2:
    input:
        counts=counts,
        samples=samples
    output:
        "sizeFactors.txt"
    run:
        if analysis_type == "Wald":
            shell(
                "Rscript scripts/deseq.R {input.counts} {input.samples}"
            )
        elif analysis_type == "LRT":
            shell(
                "Rscript scripts/deseqLRT.R {input.counts} {input.samples}"
            )
        elif analysis_type == "RIP":
            shell(
                "Rscript scripts/deseqRIP.R {input.counts} {input.samples}"
            )
        else:
            raise ValueError("analysis_type must be Wald, LRT, or RIP")

## Find clusters -- turn this on by uncommenting the rule and uncommenting the last line of the rule all; clusters.txt. Must be used with the LRT analysis
# rule find_clusters:
#     input:
#         "deseq_results_LRT.txt"
#     output:
#         "clusters.txt"
#     run:
#         shell(
#             "Rscript scripts/findClusters.R {input} {samples} {output}"
#         )