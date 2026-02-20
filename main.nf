#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// modules
include { FASTQC }             from './modules/fastqc/main.nf'
include { TRIM }               from './modules/trimmomatic/main.nf'
include { BOWTIE2_BUILD }      from './modules/bowtie2_build/main.nf'
include { BOWTIE2_ALIGN }      from './modules/bowtie2_align/main.nf'
include { SAMTOOLS_SORT }      from './modules/samtools_sort/main.nf'
include { SAMTOOLS_IDX }       from './modules/samtools_idx/main.nf'
include { SAMTOOLS_FLAGSTAT }  from './modules/samtools_flagstat/main.nf'
include { BAMCOVERAGE }        from './modules/deeptools_bamcoverage/main.nf'
include { MULTIQC }            from './modules/multiqc/main.nf'
include { MULTIBWSUMMARY }     from './modules/deeptools_multibwsummary'
include { PLOTCORRELATION }    from './modules/deeptools_plotcorrelation'
include { CALLPEAK }           from './modules/macs3_callpeak'
include { INTERSECT }          from './modules/bedtools_intersect'
include { REMOVE }             from './modules/bedtools_remove'
include { ANNOTATE }           from './modules/homer_annotatepeaks'
include { COMPUTEMATRIX }      from './modules/deeptools_computematrix'
include { PLOTPROFILE }        from './modules/deeptools_plotprofile'
include { FIND_MOTIFS }        from './modules/findmotif'

workflow {

    // loading csv file
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .filter { row -> row.name && row.path } // Filter out bad rows
        .map { row -> tuple(row.name, file(row.path)) }
        .set { raw_reads_ch }

    // loading genome
    Channel
        .value(file(params.genome))
        .set { genome_ch }

    // running the modules
    FASTQC(raw_reads_ch)

    TRIM(raw_reads_ch, file(params.adapter_fa))

    BOWTIE2_BUILD(genome_ch)

    BOWTIE2_ALIGN(
        TRIM.out.trimmed,
        BOWTIE2_BUILD.out.index,
        BOWTIE2_BUILD.out.name
    )

    SAMTOOLS_FLAGSTAT(BOWTIE2_ALIGN.out.bam)

    SAMTOOLS_SORT(BOWTIE2_ALIGN.out.bam)

    SAMTOOLS_IDX(SAMTOOLS_SORT.out)

    BAMCOVERAGE(SAMTOOLS_IDX.out)

    // merging for multiQC
    FASTQC.out.zip
        .merge(TRIM.out.log)
        .merge(SAMTOOLS_FLAGSTAT.out.flagstat)
        .set { multiqc_ch }

    MULTIQC(multiqc_ch)

    // Extract BigWig paths and collect into list
    BAMCOVERAGE.out.bigwig
        .map { it[1] }
        .collect()
        .set { bigwig_files_ch }

    // Run multiBigwigSummary
    multibw_out = MULTIBWSUMMARY(bigwig_files_ch)

    // Plot correlation 
    PLOTCORRELATION(multibw_out, Channel.value("spearman"))

    // macs3 peak calling
    // package download issues faced by class -- used provided peaks for further analysis
    //intersect using refs 
    Channel
    .fromPath('/projectnb/bf528/materials/project-2-chipseq/refs/rep1_peaks.narrowPeak')
    .combine(
        Channel.fromPath('/projectnb/bf528/materials/project-2-chipseq/refs/rep2_peaks.narrowPeak')
    )
    .map { bedA, bedB -> tuple("ChIP", bedA, bedB) }
    .set { intersect_input_ch }

    INTERSECT(intersect_input_ch)

    REMOVE(INTERSECT.out.intersect, file(params.blacklist))

    ANNOTATE(REMOVE.out.filtered, file(params.genome), file(params.gtf))
    
    COMPUTEMATRIX(
    BAMCOVERAGE.out.bigwig.filter { meta, bw -> meta.startsWith('IP') },
    file("/projectnb/bf528/materials/project-2-chipseq/refs/hg38_genes.bed"),
    2000
    )
    PLOTPROFILE(COMPUTEMATRIX.out.matrix)

    FIND_MOTIFS(REMOVE.out.filtered, file(params.genome))

}
