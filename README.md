# Description

Command-line scripts and pipelines to perform de-novo assembly with unmapped reads of human datasets 

# Dependencies

This program expects the following tools/languages to be installed as modules and be available in your path:

- fastp
- SAMtools
- MEGAHIT
- MASURCA
- KRAKEN2

...

# Installation instructions

Clone this repo to your local (linux-based) cluster environment....

# Running the pipeline

Make sure  nextflow is installed and in your path.

The basic command to run the pipeline is: 

<pre>
nextflow run Refgraph.nf
</pre>


There are a number of parameters whose values can be set or reset at the command line or via a config file.

Examples of config files are included in this repo in the conf folder.


To run the pipeline with a config file, please type:

<pre>
nextflow run -c config Refgraph.nf
</pre>


To run the pipeline and specify parameters at the command line, please type all parameters in the same line or by using the backslash to continue in the next line like this:

<pre>

nextflow run Refgraph.nf \
--genome                = "/some/path/data/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
--samplePath            = "/some/path/data/1000genome/LWK/*.cram" \
--outputDir             = "/some/path/results/1000genome_LWK_20Samples" \
--batchid               = "1000genome_LWK_20Samples"


</pre>

Required parameters:

-  --genome                = false          /*genome fasta file, must specify complete path. Required parameters*/
-  --samplePath            = false          /*input folder, must specify complete path. Required parameters*/
-  --outputDir             = false          /*output folder, must specify complete path. Required parameters*/
-  --singleEnd             = false          /*options: true|false. true=the input type is single end reads; false=the input type is paired reads. Default is false*/
-  --assembler             = 'megahit'      /*options: megahit|masurca. Default is megahit*/
-  --skipKraken2           = true           /*options: true|false. Default is true which means that kraken2 will be skipped*/


Optional parameters for readprep --  qc-trimming and adaptor removal

-  --skipTrim              = false           /*qc-trimming of reads. options: true|false. Default is false*/     
-  --min_read_length       = '20'            /*minimum length of read to be kept after trimming for downstream analysis. Default is 20*/
-  --min_base_quality      = '20'            /*minimum base quality. Default is 20*/
-  --guess_adapter         = true            /*auto-detect adapter from input file. options: true|false. Default is true*/
-  --forward_adaptor       = false           /*adapter sequence to be clipped off (forward). */
-  --reverse_adaptor       = false           /*adapter sequence to be clipped off (reverse). Used for paired reads only*.*/

# Results

The results are placed in the path specified with the outputDir parameter. 
This is the folder structure of the results folder:

<pre>
assembly_metrics
megahit
multiqc
read_prep
trimmed  
</pre>

-  The results of the first step, extract unmapped reads, go in the read_prep folder.
-  The results of the second step, qc-trimming reads, go in the trimmed folder.
-  The results of the third step, generate assembly with megahit, go in the megahit folder.
-  The results of QC step on the reads go in the multiqc folder.
-  The results of the QC step on the assembly go in the assembly_metrics folder.


# Troubleshooting

text goes here
