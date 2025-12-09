
configfile: "config.yaml"

#collect samples from config file
SAMPLES = config["samples"]

# Define the final targets of the workflow
rule all:
    input:
        # All filtered objects
        expand("results/{sample}/filtered.rds", sample=SAMPLES),
        # All QC figures
        expand("figures/{sample}/{version}_qc.png", 
               sample=SAMPLES, 
               version=["unfiltered", "filtered"]),
        # Optional: combined report
        "reports/qc_summary.html"


def get_seq_dir(wildcards):
    """
    Get the sequencing directory for a given sample from the config file.
    """
    return config["sample_dirs"][wildcards.sample]


rule seurat_initial:
    input:
        seq_dir = get_seq_dir  # Calls the function with current sample
    output:
        unfiltered = "results/{sample}/unfiltered.rds",
        filtered = "results/{sample}/filtered.rds"
    script:
        "scripts/seurat_initial.R"

