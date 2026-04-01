nextflow.enable.dsl=2

params.reads = "data/*_{1,2}.fastq.gz"
params.genome = "genome/STAR_index"
params.gtf = "annotation/genes.gtf"

workflow {

    reads_ch = Channel.fromFilePairs(params.reads, flat: true)

    fastqc_out = fastqc(reads_ch)
    multiqc(fastqc_out)

    trimmed = fastp(reads_ch)

    aligned = star(trimmed)

    counts = featureCounts(aligned)

    cleaned = remove_header(counts)

    deseq2(cleaned)
}

process fastqc {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.zip")

    script:
    """
    fastqc ${reads[0]} ${reads[1]}
    """
}

process multiqc {

    input:
    path(qc_files)

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc .
    """
}

process fastp {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1_trimmed.fastq.gz"), path("${sample_id}_2_trimmed.fastq.gz")

    script:
    """
    fastp \
      -i ${reads[0]} -I ${reads[1]} \
      -o ${sample_id}_1_trimmed.fastq.gz \
      -O ${sample_id}_2_trimmed.fastq.gz \
      --detect_adapter_for_pe
    """
}

process star {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    STAR \
      --genomeDir ${params.genome} \
      --readFilesIn $r1 $r2 \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --runThreadN 4 \
      --outFileNamePrefix ${sample_id}_

    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    """
}

process featureCounts {

    input:
    tuple val(sample_id), path(bam)

    output:
    path("counts_${sample_id}.txt")

    script:
    """
    featureCounts \
      -a ${params.gtf} \
      -o counts_${sample_id}.txt \
      -t exon \
      -g gene_id \
      -s 2 \
      -p \
      --countReadPairs \
      ${bam}
    """
}

process remove_header {

    input:
    path(counts)

    output:
    path("clean_counts.txt")

    script:
    """
    tail -n +2 ${counts} > clean_counts.txt
    """
}

process deseq2 {

    input:
    path(counts)

    output:
    path("DESeq2_results_*.csv"),
    path("plots/*")

    script:
    """
    mkdir -p plots
    Rscript bin/script_deseq2.R ${counts}
    """
}
