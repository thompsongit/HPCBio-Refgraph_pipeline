#!/usr/bin/env nextflow
params.genome = "data/genome.fa"
params.sample = "data/sample.cram"
params.assembler = 'megahit'
params.skipKraken2 = false

genome_file = file(params.genome)
sample_file = file(params.sample)

process prepare_genome{
    module 'SAMtools/1.10'
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
    module 'SAMtools/1.10'
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
    module 'SAMtools/1.10'
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

process trimming {
    module 'fastp/0.20.0'
    input:
    file '*' from fq_ch
    output:
    file '*' into trim_ch
    script:
    """
    fastp -i PEr1.fastq -o PEr1_trim.fastq -I PEr2.fastq -O PEr2_trim.fastq --unpaired1 upr1_trim.fastq --unpaired2 upr2_trim.fastq -l 20 -q 20
    fastp -i SEr1.fastq -o SEr1_trim.fastq -l 20 -q 20
    fastp -i SEr2.fastq -o SEr2_trim.fastq -l 20 -q 20
    """
}

params.megahitPath = "./results/megahit/"
if(params.assembler == 'megahit'){
process megahit_assemble {
    module 'MEGAHIT/1.2.9' 
    publishDir params.megahitPath, mode:'link'
    input:
    file '*' from trim_ch  
    output:
    file 'megahit_results' into assembly_ch

    script:
    """
    megahit -1 PEr1_trim.fastq -2 PEr2_trim.fastq -r upr1_trim.fastq,upr2_trim.fastq,SEr1_trim.fastq,SEr2_trim.fastq -o megahit_results
    """
}
}

params.masurcaPath = "./results/masurca/"
else if(params.assembler == 'masurca'){
process masurca_assemble {
    module 'MaSuRCA/3.2.3'
    publishDir params.masurcaPath, mode:'link'
    input:
    file '*' from trim_ch
    output:
    file '*' into assembly_ch
    
    script:
    """
    masurca config.txt
    bash assemble.sh
    """
}
}

else {println "Please choose assembler from megahit and masurca"}
/*
  *Megahit for different input data types:
  *megahit -1 pe_1.fq -2 pe_2.fq -o out  # 1 paired-end library
  *megahit --12 interleaved.fq -o out # one paired & interleaved paired-end library
  *megahit -1 a1.fq,b1.fq,c1.fq -2 a2.fq,b2.fq,c2.fq -r se1.fq,se2.fq -o out # 3 paired-end libraries + 2 SE libraries
  */

params.kraken2Path = "./results/Kraken2/"
library_ch = Channel.fromPath('library')
if(params.skipKraken2 == false){
process remove_contaminants {
    module 'Kraken2/2.0.8'
    publishDir params.kraken2Path, mode:'link'
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
}
