
Download the transcript sequence file here: https://www.genome.ucsc.edu/cgi-bin/hgTables?hgsid=1912921966_2XNnb1bvwfSn0mWNWurFuaasykXp
the currently downloaded file using "human-hg38-gene and gene predication-GENECODEV41-genome-sequence" then chose CDS Exons and 3 extra downstream and upstream bases and Exons in upper cases, everything else in lower case.

## after downloaded the transcrip file
## run perl generate_all_variant.pl $ftranscript $fout
## it will generate all variants' mutationrate in each transcript.
## because one variant could be involved in several transcripts, the output have repeated variants.
## for the downstream analysis
## need to annotate each variants with the required scores and variant function using the same version of transcript 
## then chose the variants fits to the requirment 
## sum up each types of variants' mutation rate per variant type per transcript.
## the mutation rate of framshift is nonsense_mutation*1.25
 
