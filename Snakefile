include: "common.smk"
from types import SimpleNamespace
pepfile: config["pepfile"]
samples = [SimpleNamespace(sample=sample) for sample in pep.sample_table["sample_name"]]
config["rg_stats_script"] = srcdir(path.join("/exports/archive/me-lcco-aml-archive/Ramin/genotyping_aml_samples/data_handler/scripts", "gather_rg_stats.py"))
config["sample_stats_script"] = srcdir(path.join("/exports/archive/me-lcco-aml-archive/Ramin/genotyping_aml_samples/data_handler/scripts", "gather_sample_stats.py"))


rule all:
    input:
        sample_json=[module_output.json(sample) for sample in samples],
        bam=[module_output.bam(sample) for sample in samples],
        vcf=[module_output.vcf(sample) for sample in samples],


rule fastqc_raw:
    """Runs FastQC for each pair of each read group of each sample given in the config file."""
    input:
        fq1=get_forward,
        fq2=get_reverse,
    output:
        folder=directory("outputs/{sample}/qc-seq/{read_group}/fastqc-{pair}-raw/"),
    threads: 4
    resources:
        mem_mb= "8000",
        runtime= "4h",
        partition= "all" 
    shell:
        """
        mkdir -p {output.folder}
        chmod -R 777 {output.folder}
        fastqc \
            --outdir={output.folder} \
            --dir=/tmp \
            --extract \
            --nogroup \
            --threads={threads} \
            {input.fq1:q} {input.fq2:q}
        """


rule cutadapt:
    """Clip fastq files"""
    input:
        fq1=get_forward,
        fq2=get_reverse,
    output:
        fq1="outputs/{sample}/qc-seq/{read_group}/{sample}-{read_group}-R1.fq.gz",
        fq2="outputs/{sample}/qc-seq/{read_group}/{sample}-{read_group}-R2.fq.gz",
    threads: 8
    resources:
        mem_mb= "8000",
        runtime= "4h",
        partition= "all"
    conda:
        "envs/cutadapt.yml"
    shell:
        """
        cutadapt \
            -a AGATCGGAAGAG \
            -A AGATCGGAAGAG \
            --cores={threads} \
            --compression-level=1 \
            --minimum-length=20 \
            --quality-cutoff=20,20 \
            --output={output.fq1} \
            --paired-output={output.fq2} \
            {input.fq1:q} {input.fq2:q}
        """


rule fastqc_processed:
    """Runs FastQC for each pair of QC-ed inputs."""
    input:
        fq1="outputs/{sample}/qc-seq/{read_group}/{sample}-{read_group}-R1.fq.gz",
        fq2="outputs/{sample}/qc-seq/{read_group}/{sample}-{read_group}-R2.fq.gz",
    output:
        folder=directory("outputs/{sample}/qc-seq/{read_group}/fastqc-{pair}-processed"),
    threads: 4
    resources:
        mem_mb= "8000",
        runtime= "4h",
        partition= "all" 
    shell:
        """
        mkdir -p {output.folder}
        chmod 777 -R {output.folder}
        fastqc  \
            --outdir={output.folder} \
            --dir=/tmp \
            --extract \
            --nogroup \
            --threads={threads} \
            {input.fq1:q} {input.fq2:q}
        """

rule rg_stats:
    """Gathers read statistics on the read group level."""
    input:
        raw1="outputs/{sample}/qc-seq/{read_group}/fastqc-R1-raw",
        raw2="outputs/{sample}/qc-seq/{read_group}/fastqc-R2-raw",
        proc1="outputs/{sample}/qc-seq/{read_group}/fastqc-R1-processed",
        proc2="outputs/{sample}/qc-seq/{read_group}/fastqc-R2-processed",
        rg_stats_script=config["rg_stats_script"],
    output:
        stats="outputs/{sample}/qc-seq/{read_group}/stats.json",
    resources:
        mem_mb= "8000",
        runtime= "4h",
        partition= "all"
    shell:
        """
        python {input.rg_stats_script} \
            --name {wildcards.read_group} \
            {input.raw1} {input.raw2} \
            {input.proc1} {input.proc2} > {output.stats} 
        """


rule sample_stats:
    """Gathers read statistics on the sample level."""
    input:
        rg_stats=get_sample_stats,
        sample_stats_script=config["sample_stats_script"],
    output:
        stats="outputs/{sample}/qc-seq/{sample}.seq_stats.json",
    resources:
        mem_mb= "8000",
        runtime= "4h",
        partition= "all"
    shell:
        """
        python {input.sample_stats_script} \
            {input.rg_stats} > {output.stats} 
        """


rule merge_fastqs_r1:
    """Merges all FASTQ files for a given sample from its read groups."""
    input:
        fqs=get_all_trimmed_forward,
    output:
        merged=temp("outputs/{sample}/{sample}-R1.fq.gz")
    resources:
        mem_mb= "1000",
        runtime= "1h",
        partition= "short"
    shell:
        """
        cp {input.fqs} {output.merged} \
        || \
        cat {input.fqs} > {output.merged} 
        """


rule merge_fastqs_r2:
    """Merges all FASTQ files for a given sample from its read groups."""
    input:
        fqs=get_all_trimmed_reverse,
    output:
        merged=temp("outputs/{sample}/{sample}-R2.fq.gz")
    resources:
        mem_mb= "1000",
        runtime= "1h",
        partition= "short"
    shell:
        """
        cp {input.fqs} {output.merged}  \
        || \
        cat {input.fqs} > {output.merged} 
        """


rule merge_fastqs_raw_r1:
    """Merges all raw FASTQ files for a given sample from its read groups."""
    input:
        fqs=get_all_r1,
    output:
        merged=temp("outputs/{sample}/{sample}-R1.raw.fq.gz"),
    resources:
        mem_mb= "1000",
        runtime= "1h",
        partition= "short"
    shell:
        """
        cp {input.fqs:q} {output.merged}  \
        || \
        cat {input.fqs:q} > {output.merged} 
        """


rule merge_fastqs_raw_r2:
    """Merges all raw FASTQ files for a given sample from its read groups."""
    input:
        fqs=get_all_r2,
    output:
        merged=temp("outputs/{sample}/{sample}-R2.raw.fq.gz"),
    resources:
        mem_mb= "1000",
        runtime= "1h",
        partition= "short"
    shell:
        """
        cp {input.fqs:q} {output.merged}  \
        || \
        cat {input.fqs:q} > {output.merged}
        """
    

rule star_map:
    input:
        ref="/exports/archive/me-lcco-aml-archive/Ramin/alternative_splicing/data/reference-genome/reference.fa",
        genomedir="/exports/archive/me-lcco-aml-archive/Ramin/T-cells/data/reference-genome/genomeIndex2",
        fq1="outputs/{sample}/{sample}-R1.fq.gz",
        fq2="outputs/{sample}/{sample}-R2.fq.gz",
        gtf="/exports/archive/me-lcco-aml-archive/Ramin/T-cells/data/reference-genome/gencode.v42.chr_patch_hapl_scaff.annotation.gtf",
    output:
        bam="outputs/{sample}/snv-indels/Aligned.sortedByCoord.out.bam",
        gene_count="outputs/{sample}/snv-indels/ReadsPerGene.out.tab",
    threads:
        8
    resources:
        mem_mb= "64000",
        runtime= "24h",
        partition= "all"
    params:
        twopassMode="Basic",
        chim_segment=20,
        min_intron_size=50,
    threads:
        8
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.genomedir} \
            --sjdbGTFfile {input.gtf} \
            --readFilesCommand zcat \
            --twopassMode Basic \
            --outFileNamePrefix $(dirname {output.bam})/ \
            --outSAMtype BAM SortedByCoordinate \
            --alignIntronMin {params.min_intron_size} \
            --chimOutType WithinBAM \
            --chimSegmentMin {params.chim_segment} \
            --quantMode GeneCounts \
            --readFilesIn {input.fq1:q} {input.fq2:q}
        """

rule index_bamfile:
    input:
        bam="outputs/{sample}/snv-indels/Aligned.sortedByCoord.out.bam",
    output:
        bai="outputs/{sample}/snv-indels/Aligned.sortedByCoord.out.bam.bai",
    resources:
        mem_mb= "1000",
        runtime= "1h",
        partition= "short"
    shell:
        """
        samtools index {input.bam}
        """

rule add_read_group:
    input:
        bam="outputs/{sample}/snv-indels/Aligned.sortedByCoord.out.bam",
        index="outputs/{sample}/snv-indels/Aligned.sortedByCoord.out.bam.bai"
    output:
        temp("outputs/{sample}/snv-indels/read_grp_fixed.bam")
    resources:
        mem_mb= "4000",
        runtime= "2h",
        partition= "all"
    shell:
        """
        picard AddOrReplaceReadGroups -I {input.bam} -O {output} --RGID 1 --RGPL illumina --RGPU unit1 --RGSM 3 --RGLB lib2
        """

rule picard_index_read_group:
    input:
        "outputs/{sample}/snv-indels/read_grp_fixed.bam"
    output:
        temp("outputs/{sample}/snv-indels/read_grp_fixed.bam.bai")
    resources:
        mem_mb= "1000",
        runtime= "1h",
        partition= "short"
    shell:
        "picard BuildBamIndex -I {input} -O {output}"

rule splice_handling:
    input:
        bam="outputs/{sample}/snv-indels/read_grp_fixed.bam",
        index="outputs/{sample}/snv-indels/read_grp_fixed.bam.bai",
        ref="/exports/archive/me-lcco-aml-archive/Ramin/alternative_splicing/data/reference-genome/reference.fa",
    output:
        temp("/exports/me-lcco-aml-hpc/Ramin/genotyping_pipeline/outputs/{sample}/snv-indels/split_reads.bam")
    threads: 8
    resources:
        mem_mb= "10000",
        runtime= "20h",
        partition= "all"
    shell:
        """
        gatk --java-options '-XX:ConcGCThreads={threads} -Xmx8g' SplitNCigarReads -R {input.ref} -I {input.bam} -O {output} 
        """

rule picard_index_split_reads:
    input:
        "/exports/me-lcco-aml-hpc/Ramin/genotyping_pipeline/outputs/{sample}/snv-indels/split_reads.bam"
    output:
        temp("/exports/me-lcco-aml-hpc/Ramin/genotyping_pipeline/outputs/{sample}/snv-indels/split_reads.bam.bai")
    resources:
        mem_mb= "1000",
        runtime= "1h",
        partition= "short"
    shell:
        "picard BuildBamIndex -I {input} -O {output}"


rule variant_calling:
    input:
        bam="/exports/me-lcco-aml-hpc/Ramin/genotyping_pipeline/outputs/{sample}/snv-indels/split_reads.bam",
        index="/exports/me-lcco-aml-hpc/Ramin/genotyping_pipeline/outputs/{sample}/snv-indels/split_reads.bam.bai",
        ref="/exports/archive/me-lcco-aml-archive/Ramin/alternative_splicing/data/reference-genome/reference.fa",
    output:
        vcf="/exports/me-lcco-aml-hpc/Ramin/genotyping_pipeline/outputs/{sample}/snv-indels/{sample}.raw.vcf",
        bam="/exports/me-lcco-aml-hpc/Ramin/genotyping_pipeline/outputs/{sample}/snv-indels/variants.bam"
    threads: 8
    resources:
        mem_mb= "10000",
        runtime= "20h",
        partition= "all"
    shell:
        """
        gatk --java-options '-XX:ConcGCThreads={threads} -Xmx8g' HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} -bamout {output.bam}
        """

