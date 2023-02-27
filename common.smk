from types import SimpleNamespace
from os import path
pepfile: config["pepfile"]


config["rg_stats_script"] = srcdir(path.join("data_handler/scripts", "gather_rg_stats.py"))
config["sample_stats_script"] = srcdir(path.join("data_handler/scripts", "gather_sample_stats.py"))

def get_fastq(wildcards, pair):
    """Get the input fastq file (R1 or R2), based on the wildcards
    Uses wildcards.sample to determine the saple
    Uses wildcards.read_group to select the correct fastq
    """
    fastq = get_fastq_pep(wildcards.sample, pair)
    # Here, we use the readgroup wildcard to pick the correct fastq file
    read_groups = {
        rg: fastq for rg, fastq in zip(get_readgroups(wildcards.sample), fastq)
    }
    return read_groups[wildcards.read_group]


def get_forward(wildcards):
    """Get the input forward fastq file"""
    return get_fastq(wildcards, "R1")


def get_reverse(wildcards):
    """Get the input reverse fastq file"""
    return get_fastq(wildcards, "R2")


def get_all_trimmed_fastq(wildcards, pair):
    """Get the filenames of all trimmed fastq files for sample"""
    sample = wildcards.sample
    return [
        f"outputs/{sample}/qc-seq/{rg}/{sample}-{rg}-{pair}.fq.gz"
        for rg in get_readgroups(sample)
    ]

def get_all_trimmed_forward(wildcards):
    """Get the filenames for all forward fastq files for sample"""
    return get_all_trimmed_fastq(wildcards, "R1")


def get_all_trimmed_reverse(wildcards):
    """Get the filenames for all revrse fastq files for sample"""
    return get_all_trimmed_fastq(wildcards, "R2")


def get_fastq_pep(sample, pair):
    """Get all fastq files for the specified sample"""
    fastq = pep.sample_table.loc[sample, pair]
    # If a single fastq file is specified, we put it in a list
    if isinstance(fastq, str):
        fastq = [fastq]
    return fastq


def get_all_r1(wildcards):
    """Get all forward fastq files for the specified sample"""
    return get_fastq_pep(wildcards.sample, "R1")


def get_all_r2(wildcards):
    """Get all forward fastq files for the specified sample"""
    return get_fastq_pep(wildcards.sample, "R2")


def get_readgroups(sample):
    """Get the readgroup names per sample"""
    return [f"rg_{i+1}" for i in range(len(get_fastq_pep(sample, "R1")))]
    #return pep.sample_table.loc[sample, "subsample_name"]


def get_readgroup_per_sample():
    """Get all combinations of readgroups and samples"""
    for sample in pep.sample_table["sample_name"]:
        for readgroup in get_readgroups(sample):
            yield readgroup, sample


def get_sample_stats(wildcards):
    """Get the stat files for every readgroup in wildcards.sample"""
    sample = wildcards.sample
    return [f"outputs/{sample}/qc-seq/{rg}/stats.json" for rg in get_readgroups(sample)]


## Functions for module outputs ##


def get_forward_output(wildcards):
    return get_output(wildcards, "R1")


def get_reverse_output(wildcards):
    return get_output(wildcards, "R2")


def get_output(wildcards, pair):
    return f"outputs/{wildcards.sample}/{wildcards.sample}-{pair}.fq.gz"


def get_raw_forward_output(wildcards):
    return get_raw_output(wildcards, "R1")


def get_raw_reverse_output(wildcards):
    return get_raw_output(wildcards, "R2")


def get_raw_output(wildcards, pair):
    return f"outputs/{wildcards.sample}/{wildcards.sample}-{pair}.raw.fq.gz"


def get_sample_json(wildcards):
    return f"outputs/{wildcards.sample}/qc-seq/{wildcards.sample}.seq_stats.json"

def get_vcf(wildcards):
    return f"/exports/me-lcco-aml-hpc/Ramin/genotyping_pipeline/outputs/{wildcards.sample}/snv-indels/{wildcards.sample}.raw.vcf"

def get_bam(wildcards):
    return f"/exports/me-lcco-aml-hpc/Ramin/genotyping_pipeline/outputs/{wildcards.sample}/snv-indels/variants.bam"

module_output = SimpleNamespace(
    forward=get_forward_output,
    reverse=get_reverse_output,
    forward_raw=get_raw_forward_output,
    reverse_raw=get_raw_reverse_output,
    json=get_sample_json,
    vcf=get_vcf,
    bam=get_bam,
)
