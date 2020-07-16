#!/usr/bin/env nextflow
params.genome = "data/genome.fa"
params.sample = "data/sample.cram"
params.assembler = 'megahit'

genome_file = file(params.genome)
sample_file = file(params.sample)

process prepare_genome{
    module 'SAMtools'
    input:
    file genome from genome_file
    output:
    file "${genome}.fai" into genome_index_ch
    
    script:
    """
    samtools faidx ${genome}
    """
}

process extract_unmap {
    module 'SAMtools'
    input:
    file sample from sample_file	
    file genome from genome_file
    file index from genome_index_ch
    output:
    file 'unmap.bam' into bam_ch

    script:
    """
    samtools view -bt index -S -b -f 4 sample > unmap.bam
    """
}

process bamtofastq {
    module 'SAMtools'
    input:
    file 'unmap.bam' from bam_ch
    output:
    file '*' into fq_ch

    script:
    """
    samtools fastq -f 12 unmap.bam -1 PEr1.fastq -2 PEr2.fastq
    samtools fastq -f 68 -F 8 unmap.bam > SEr1.fastq
    samtools fastq -f 132 -F 8 unmap.bam > SEr2.fastq
    """
}

if(params.assembler == 'megahit'){
process megahit_assemble {
    module 'MEGAHIT' 
    input:
    file '*' from fq_ch  
    output:
    file 'megahit_results' into assembly_ch

    script:
    """
    megahit -1 PEr1.fastq -2 PEr2.fastq -r SEr1.fastq,SEr2.fastq -o megahit_results
    """
}
}

/*
  *Megahit for different input data types:
  *megahit -1 pe_1.fq -2 pe_2.fq -o out  # 1 paired-end library
  *megahit --12 interleaved.fq -o out # one paired & interleaved paired-end library
  *megahit -1 a1.fq,b1.fq,c1.fq -2 a2.fq,b2.fq,c2.fq -r se1.fq,se2.fq -o out # 3 paired-end libraries + 2 SE libraries
  */

library_ch = Channel.fromPath('library')
process remove_contaminants {
    module 'Kraken2'
    input:
    file 'library' from library_ch
    file 'megahit_results' from assembly_ch

    output:
    file 'contigs_kraken2.txt' into classification_ch
    file 'sample_specific_contigs.fa' into result
          
    script:
    """
    kraken2 --db library megahit_results/final.contigs.fa > contigs_kraken2.txt
    grep U contigs_kraken2.txt > unclassified.txt
    awk '{print $2}' unclassified.txt > unclassified_readID.txt
    grep -w -A 1 -f unclassified_readID.txt megahit_results/final.contigs.fa --no-group-separator > sample_specific_contigs.fa
    """
}

