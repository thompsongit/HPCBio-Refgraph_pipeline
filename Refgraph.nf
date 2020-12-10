#!/usr/bin/env nextflow

/*parameters that are specified at the command line or via config file*/
params.genome                = false          /*genome fasta file, must specify complete path. Required parameters*/
params.samplePath            = false          /*input folder, must specify complete path. Required parameters*/
params.outputDir             = false          /*output folder, must specify complete path. Required parameters*/
params.singleEnd             = false          /*options: true|false. true = the input type is single end reads; false = the input type is paired reads. Default is false*/
params.assembler             = 'megahit'      /*options: megahit|masurca. Default is megahit*/
params.skipKraken2           = true           /*options: true|false. Default is true which means that kraken2 will be skipped*/

/* parameters for readprep = qctrimming and adaptor removal */
params.skipTrim              = false           /*qc-trimming of reads. options: true|false. Default is false*/     
params.min_read_length       = '20'            /*minimum length of read to be kept after trimming for downstream analysis. Default is 20*/
params.min_base_quality      = '20'            /*minimum base quality. Default is 20*/
params.guess_adapter         = true            /*auto-detect adapter from input file. options: true|false. Default is true*/
params.forward_adaptor       = false           /*adapter sequence to be clipped off (forward). */
params.reverse_adaptor       = false           /*adapter sequence to be clipped off (reverse). Used for paired reads only*.*/


/*output folder paths*/
readPrepPath                 = "${params.outputDir}/read_prep"
trimPath                     = "${params.outputDir}/trimmed"
megahitPath                  = "${params.outputDir}/megahit"
masurcaPath                  = "${params.outputDir}/masurca"
kraken2Path                  = "${params.outputDir}/kraken2"
multiqcPath                  = "${params.outputDir}/multiqc"
metricsPath                  = "${params.outputDir}/assembly_metrics"

/*cluster parameters */
myExecutor                   = 'slurm'
params.myQueue               = 'normal'
defaultCPU                   = '1'
defaultMemory                = '20'
assemblerCPU                 = '12'
assemblerMemory              = '100'
// params.clusterAcct           = " -A h3bionet "

/*software stack*/
params.perlMod               = 'Perl/5.24.1-IGB-gcc-4.9.4'
params.fastpMod              = 'fastp/0.20.0-IGB-gcc-4.9.4'
params.samtoolsMod           = 'SAMtools/1.10-IGB-gcc-8.2.0'
params.megahitMod            = 'MEGAHIT/1.2.9-IGB-gcc-8.2.0'
params.assemblyMetricsMod    = "BBMap/38.36-Java-1.8.0_152"
params.multiqcMod            = "MultiQC/1.7-IGB-gcc-4.9.4-Python-3.6.1"

/*Prepare input*/
genome_file                  = file(params.genome)
genomeStore                  = genome_file.getParent()
if( !genome_file.exists() ) exit 1, "Missing reference genome file: ${genome_file}"
CRAM_Ch1 = Channel.fromPath("${params.samplePath}")

/*

  prepare_genome 
  This process is executed only once

*/


process prepare_genome{
    tag                    { genome }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    storeDir               genomeStore
    validExitStatus        0
    
    input:
    file genome from genome_file

    output:
    file "*.fai" into genome_index_ch
    
    script:
    """
    samtools faidx ${genome}
    """
}

/*

  extract_unmap 
  The input files are in CRAM format which have been previously aligned to genome prepared in previous step

*/

process qc_input {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    publishDir             readPrepPath, mode: "copy"	    
    validExitStatus        0,1
    errorStrategy          'ignore'
    stageOutMode           'copy'
    
    input:
    set val(name), file(CRAM) from CRAM_Ch1	

    output:
    set val(name), file('*_ok.cram') optional true into extract_ch
    
    script:
    """
    samtools quickcheck ${name}
    if [ \$? -eq 0 ]
    then
        cp ${name} ${name.baseName}_ok.cram
    fi
    """

}


process extract_unmap {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    publishDir             readPrepPath, mode: "copy"	    
    validExitStatus        0,1
    errorStrategy          'finish'
    scratch                '/scratch'
    stageOutMode           'copy'
    
    input:
    set val(name), file(okcram) from extract_ch	
    file genome from genome_file
    file index from genome_index_ch

    output:
    set val(name), file('*.fastq') optional true into fq_ch
    file '*'
    
    script:
    if(params.singleEnd){
    """
    samtools view -bt ${index} -S -b -f 4 ${name} > ${name.baseName}.unmap.bam
    if [ -s ${name.baseName}.unmap.bam ]
    then
	      samtools fastq -f 4 ${name.baseName}.unmap.bam > ${name.baseName}_SE_R1.fastq
    fi
    """
    } else {
    """
    samtools view -bt ${index} -S -b -f 4 ${name} > ${name.baseName}.unmap.bam
    if [ -s ${name.baseName}.unmap.bam ]
    then
         samtools fastq -f 12 ${name.baseName}.unmap.bam -1 ${name.baseName}_PE_R1.fastq -2 ${name.baseName}_PE_R2.fastq
    fi
    #samtools fastq -f 68 -F 8 ${name.baseName}.unmap.bam > ${name.baseName}_SE_R1.fastq
    #samtools fastq -f 132 -F 8 ${name.baseName}.unmap.bam > ${name.baseName}_SE_R2.fastq
    """    
    }

}

/*

  trimming 

*/

process trimming {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             trimPath, mode: "copy"
    module                 params.fastpMod
    validExitStatus        0,1
    errorStrategy          'finish'
    scratch                '/scratch'
    stageOutMode           'copy'

    input:
    set val(name), file(reads) from fq_ch
    
    output:
    set val(name), file('*.trimmed.fq') optional true into trim_ch
    set val(name), file('*.json') optional true into multiqc_ch
    file '*'
	    
    script:
    trimOptions      = params.skipTrim ? ' ' :  ' -l 20 -q 20  --cut_right --cut_right_window_size 3 --cut_right_mean_quality 20 '
    adapterOptionsSE = params.guess_adapter ? ' ' : ' --adapter_sequence="${params.forward_adaptor}" '
    adapterOptionsPE = params.guess_adapter ? ' --detect_adapter_for_pe ' : ' --detect_adapter_for_pe --adapter_sequence="${params.forward_adaptor}"  --adapter_sequence_r2="${params.reverse_adaptor}"  '

    if(params.singleEnd){
    """
    fastp --in1 ${reads[0]} --out1 "${name.baseName}.SE.R1.trimmed.fq" ${adapterOptionsSE} ${trimOptions} --thread ${task.cpus} -w ${task.cpus} --html "${name.baseName}"_SE_fastp.html --json "${name.baseName}"_SE_fastp.json
    """
    } else {
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 "${name.baseName}.PE.R1.trimmed.fq"  --out2 "${name.baseName}.PE.R2.trimmed.fq" --unpaired1 "${name.baseName}.unpR1.trimmed.fq" --unpaired2 "${name.baseName}.unpR2.trimmed.fq" ${adapterOptionsPE}  ${trimOptions} --thread ${task.cpus} -w ${task.cpus}  --html "${name.baseName}"_PE_fastp.html --json "${name.baseName}"_PE_fastp.json
    """
    }
}



/*
  *Megahit for different input data types:
  *megahit -1 pe_1.fq -2 pe_2.fq -o out  # 1 paired-end library
  *megahit --12 interleaved.fq -o out # one paired & interleaved paired-end library
  *megahit -1 a1.fq,b1.fq,c1.fq -2 a2.fq,b2.fq,c2.fq -r se1.fq,se2.fq -o out # 3 paired-end libraries + 2 SE libraries
*/


process megahit_assemble {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   assemblerCPU
    queue                  params.myQueue
    memory                 "$assemblerMemory GB"
    module                 params.megahitMod 
    publishDir             megahitPath , mode:'copy'
    validExitStatus        0,1
    errorStrategy          'finish'
    scratch                '/scratch'
    stageOutMode           'copy'
    
    input:
    set name, file(fastqs)  from trim_ch  

    output:
    set val(name), file('*.final.contigs.fa') optional true into assembly_ch
    file '*'

    script:
    if(params.singleEnd){
    """
    megahit -1 ${fastqs[0]} -o ${name.baseName}.megahit_results
    cp ${name.baseName}.megahit_results/final.contigs.fa ${name.baseName}.final.contigs.fa
    """
    } else {
    """
    megahit -1 ${fastqs[0]} -2 ${fastqs[1]}  -o ${name.baseName}.megahit_results
    cp ${name.baseName}.megahit_results/final.contigs.fa ${name.baseName}.final.contigs.fa

    """
    }
}


process calc_assembly_metrics {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.assemblyMetricsMod 
    publishDir             metricsPath , mode:'copy', overwrite: true
    
    input:
    set name, file(contig)  from assembly_ch  

    output:
    set val(name), file('*.stats') optional true into metrics_ch
    file '*'

    script:
    """
    stats.sh ${contig} > ${name.baseName}.stats
    """

}
/*

  masurca
  To be added later

*/

/*

  kraken2
  To be added later

*/

/*

  QC all samples

*/
process MultiQC_readPrep {
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.multiqcMod
    publishDir             multiqcPath, mode: 'copy', overwrite: true
 
    input:
    file('*') from multiqc_ch.collect()

    output:
    file "*"

    """
    multiqc .
    """
} 

process Assembly_metrics {
    executor               myExecutor
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.perlMod
    publishDir             metricsPath, mode: 'copy', overwrite: true
 
    input:
    file('*') from metrics_ch.collect()

    output:
    file "*"

    """
    cat *.stats > assembly_metrics.txt
    """
} 

