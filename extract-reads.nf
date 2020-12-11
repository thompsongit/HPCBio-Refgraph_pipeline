#!/usr/bin/env nextflow

params.genome                = false          /*genome fasta file, must specify complete path. Required parameters*/
params.samplePath            = false          /*input folder, must specify complete path. Required parameters*/
params.outputDir             = false          /*output folder, must specify complete path. Required parameters*/
params.singleEnd             = false

readPrepPath                 = "${params.outputDir}/read_prep"
trimPath                     = "${params.outputDir}/trimmed"

myExecutor                   = 'slurm'
params.myQueue                      = 'normal'
defaultCPU                   = '1'
defaultMemory                = '20'
assemblerCPU                 = '12'
assemblerMemory              = '500'

params.samtoolsMod           = 'SAMtools/1.10-IGB-gcc-8.2.0'

/*Prepare input*/
genome_file                  = file(params.genome)
genomeStore                  = genome_file.getParent()
if( !genome_file.exists() ) exit 1, "Missing reference genome file: ${genome_file}"
CRAM_Ch1 = Channel.fromFilePairs("${params.samplePath}", size: 1)

/*

  prepare_genome 
  This process is executed only once
  NOTE: This process should remain the same as in the original workflow, will be merged in

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

   Read extraction.  This step is tricky.
   
   For single end reads, we only need to extract unmapped as there isn't a mapped mate.  
   However this also means that the additional workflows downstream cannot be run as they 
   paired end read data (location placement requires having mapped mates)
   
   samtools view -f 4 > se_unmapped.bam
   For paired end reads, we need to extract three classes of reads, keeping the mates
   
   1. Both unmapped
   2. R1 mapped, R2 unmapped
   3. R2 mapped, R1 unmapped
   
   This can be accomplished with samtools though the logic is a little tricky with the bit
   flags.  Best way would be to extract pairs where any of the reads are unmapped, then 
   pull out subsets using flags. Otherwise we're running through the 
   
   # both unmapped; -f means unmapped, 
   # bit flag 12  = both reads unmapped (bit flag 4 & 8)
   # bit flag 2304 = not primary alignment (this removes these), not supplemental
   samtools view -hb -f 12 -F 2304 alignments.bam > both_unmapped.bam
    
   # R1 mapped, R2 not
   # bit flag 4   = R1 unmapped
   # bit flag 2312 = Mate unmapped and not primary alignment (removes these), not supplemental
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                set
   samtools view -hb -f 4 -F 2312 alignments.bam  > R1_unmapped.bam
   
   # R2 mapped, R1 not
   # bit flag 8    = R2 (mate) unmapped
   # bit flag 2308 = R1 (read) unmapped and not primary alignment (removes these), not supplemental
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                set
   samtools view -hb -f 8 -F 2308 alignments.bam  > R2_unmapped.bam

*/

/*

  extract_unmap 
  The input files are in CRAM format which have been previously aligned to genome prepared in previous step

*/

process extract_unmapped {
    tag                    { name }
    executor               myExecutor
    cpus                   4
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    publishDir             readPrepPath, mode: "copy"
    validExitStatus        0,1
    errorStrategy          'finish'
    scratch                '/scratch'
    stageOutMode           'copy'
    
    input:
    set val(id), file(cram) from CRAM_Ch1	
    file genome from genome_file
    file index from genome_index_ch

    output:
    set val(id), file("${id}.both-unmapped.R{1,2}.fastq") optional true into fq_pe_ch
    set val(id), file("${id}.orphans.unmapped.fastq") optional true into fq_se_ch
    file "${id}.improper.bam" optional true // all improper pairs; we save this for now, might be useful later
    set val(id), file("${id}.unmapped.bam") optional true into unmapped_bam_ch  // all unaligned + mates
    
    script:
    if(params.singleEnd) {
    
    """
    # TODO: UNTESTED!!!!
    # we only need to extract reads that are unmapped, no worries about pairing
    samtools view -@ ${defaultCPU} -hbt ${index} -f 4 -o ${id}.unmapped.bam ${cram} 
    
    # convert to FASTQ
    samtools fastq -@ ${defaultCPU} ${id}.SE.unmapped.bam > ${id}.orphans.unmapped.fastq
    """
    
    } else {
    
    """
    # two stages; grab any non-properly paired reads (includes unmapped)
    samtools view -@ ${task.cpus} -hbt ${index} -G 2 -o ${id}.improper.bam ${cram}
    
    # now capture the three classes; this is much faster than parsing the full BAM each time
    
    # both unmapped
    samtools view -@ ${task.cpus} -hbt ${index} \\
        -f 12 -F 2304 -o ${id}.both-unmapped.bam ${id}.improper.bam

    samtools fastq -@ ${task.cpus} ${id}.both-unmapped.bam \\
        -1 ${id}.both-unmapped.R1.fastq -2 ${id}.both-unmapped.R2.fastq
        
    # R1 only unmapped
    samtools view -@ ${task.cpus} -hbt ${index} \\
        -f 4 -F 2312 -o ${id}.R1-unmapped.bam ${id}.improper.bam

    samtools fastq -@ ${task.cpus} ${id}.R1-unmapped.bam \\
        -1 ${id}.R1-unmapped.R1.fastq -2 ${id}.R1-unmapped.R2.fastq

    # R2 only unmapped
    samtools view -@ ${task.cpus} -hbt ${index} \\
        -f 8 -F 2308 -o ${id}.R2-unmapped.bam ${id}.improper.bam
        
    samtools fastq -@ ${task.cpus} ${id}.R1-unmapped.bam \\
        -1 ${id}.R2-unmapped.R1.fastq -2 ${id}.R2-unmapped.R2.fastq
    
    # combine unmapped BAM files for later analysis
    samtools merge -@ ${task.cpus} ${id}.unmapped.bam ${id}.both-unmapped.bam ${id}.R1-unmapped.bam ${id}.R2-unmapped.bam
    
    # combined SE unmapped FASTQ into one file for assembly
    cat ${id}.R1-unmapped.R1.fastq ${id}.R2-unmapped.R2.fastq >> ${id}.orphans.unmapped.fastq
    
    # we could combine the mapped reads here, but we'll use the original BAM instead
    ## cat ${id}.R1-unmapped.R2.fastq ${id}.R2-unmapped.R1.fastq >> ${id}.orphans.mapped.fastq
    
    """
    
    }

}

