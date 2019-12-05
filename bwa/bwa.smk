rule map_paired_reads_with_bwa:
    conda:
        pipeline_path + "envs/bwa-samtools.yml"
    params:
        db = lambda w: "references/bwa_db/{}/all.fna".format(w.ref),
        platform = "ILLUMINA"
    input:
        fastq1 = out_dir + "{sample}/clean/{sample}_R1.filHuman.fastq",
        fastq2 = out_dir + "{sample}/clean/{sample}_R2.filHuman.fastq",
    output:
        out_dir + "{sample}/mapping/bwa/{ref}/paired.{ref}.bam",
    benchmark:
        out_dir + "benchmarks/{sample}/map_paired_reads_with_bwa_{ref}.txt"
    log:
        logging_folder+"{sample}/logs/mapping/bwa/{ref}/paired.log.txt"
    shell:
        """
        if ls {output[0]}.tmp* 1> /dev/null 2>&1
        then
            rm {output[0]}.tmp*
        fi
        if ls {out_dir}/{wildcards.sample}/mapping/bwa/{wildcards.ref}/*.bai 1> /dev/null 2>&1
        then
            rm -f {out_dir}/{wildcards.sample}/mapping/bwa/{wildcards.ref}/*.bai
        fi
        (bwa mem -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:{params[platform]}' {params[db]} {input[fastq1]} {input[fastq2]} -v 1 | samtools sort -O BAM -o {output[0]}) 2> {log}
        """

rule map_single_reads_with_bwa:
    conda:
        pipeline_path + "envs/bwa-samtools.yml"
    params:
        db = lambda w: "references/bwa_db/{}/all.fna".format(w.ref),
        platform = "ILLUMINA"
    input:
        fastq = out_dir + "{sample}/clean/single_{sample}.filHuman.fastq",
    output:
        out_dir + "{sample}/mapping/bwa/{ref}/single.{ref}.bam",
    benchmark:
        out_dir + "benchmarks/{sample}/map_single_reads_with_bwa_{ref}.txt"
    log:
        logging_folder+"{sample}/logs/mapping/bwa/{ref}/single.log.txt"
    shell:
        """
        if ls {output[0]}.tmp* 1> /dev/null 2>&1
        then
             rm {output[0]}.tmp*
        fi
        if ls {out_dir}/{wildcards.sample}/mapping/bwa/{wildcards.ref}/*.bai 1> /dev/null 2>&1
        then
            rm -f {out_dir}/{wildcards.sample}/mapping/bwa/{wildcards.ref}/*.bai
        fi
        (bwa mem -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:{params.platform}' {params[db]} {input[fastq]} -v 1 | samtools sort -O BAM -o {output[0]}) 2> {log}
        """

rule filter_reads_on_quality:
    conda:
        pipeline_path + "envs/bwa-samtools.yml"
    input:
        bam = out_dir + "{sample}/mapping/bwa/{ref}/{reads_type}.{ref}.bam",
    params:
        db = lambda w: "references/bwa_db/{}/all.fna".format(w.ref),
    output:
        out_dir + "{sample}/mapping/bwa/{ref}/{sample}.{reads_type}.{ref}_filtered.bam",
    benchmark:
        out_dir + "benchmarks/{sample}/{reads_type}.{ref}.filter_reads_on_quality.txt"
    shell:
        """
        if ls {output[0]}.tmp* 1> /dev/null 2>&1
        then
            rm {output[0]}.tmp*
        fi
        if ls {out_dir}/{wildcards.sample}/mapping/bwa/{wildcards.ref}/*.bai 1> /dev/null 2>&1
        then
            rm -f {out_dir}/{wildcards.sample}/mapping/bwa/{wildcards.ref}/*.bai
        fi
        samtools view -q 60 -b -u -F 3588 -T {params[db]} {input[bam]} | samtools sort -O bam -o {output[0]}
        """

rule bam_stat:
    conda:
        pipeline_path + "envs/pysam.yml"
    input:
        in_bam = out_dir + "{sample}/mapping/bwa/{ref}/{sample}.{reads_type}.{ref}_filtered.bam"
    params:
        logging_level = logging_level
    output:
        out_stat = out_dir + "{sample}/mapping/bwa/{ref}/{reads_type}.stat.tsv"
    benchmark:
        out_dir + "benchmarks/{sample}/{ref}/{reads_type}.bam_stat.txt"
    log:
        logging_folder + "{sample}/{ref}/{reads_type}_bam_stat.log"
    script:
        "scripts/bwa/bam_stat.py"

if os.path.basename(config["link_directory"]) == 'single_links':
    rule sum_stat_on_genome:
        input:
            stat_tsv = out_dir + "{sample}/mapping/bwa/{ref}/{reads_type}.stat.tsv",
            reads_stat = out_dir + "{sample}/clean/single_stat.tsv"
        params:
            refseq_summary_dir = config["refseq_summary_dir"],
            logging_level = logging_level
        output:
            sum_tsv = out_dir + "report/{sample}/{ref}.{reads_type}.stat.sum.tsv"
        benchmark:
            out_dir + "benchmarks/{sample}/{ref}.{reads_type}.sum_stat_on_genome.txt"
        log:
            logging_folder + "{sample}/{ref}.{reads_type}.sum_stat_on_genome.log"
        script:
            "scripts/bwa/sum_genome_coverage.py"
else:
    rule sum_stat_on_genome:
        input:
            stat_tsv = out_dir + "{sample}/mapping/bwa/{ref}/{reads_type}.stat.tsv"
        params:
            refseq_summary_dir = config["refseq_summary_dir"],
            logging_level = logging_level
        output:
            sum_tsv = out_dir + "report/{sample}/{ref}.{reads_type}.stat.sum.tsv"
        benchmark:
            out_dir + "benchmarks/{sample}/{ref}.{reads_type}.sum_stat_on_genome.txt"
        log:
            logging_folder + "{sample}/{ref}.{reads_type}.sum_stat_on_genome.log"
        script:
            "scripts/bwa/sum_genome_coverage.py"
