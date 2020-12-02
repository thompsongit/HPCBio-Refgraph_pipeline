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

# Results

text goes here

# Trobleshooting

text goes here
