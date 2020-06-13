##########################################################################
#		 		pipeline for haploid assembly validation
##########################################################################

#step 1. allele separation vis BLAST
makeblastdb -in apple_phased_diploid_scaffolds.fasta -dbtype nucl
	
blastn -task megablast -query apple_phased_diploid_scaffolds.fasta -db apple_phased_diploid_scaffolds.fasta -evalue 1e-20 -num_threads 60 -perc_identity 90 -word_size 100 -max_target_seqs 5 -outfmt "6 qseqid sseqid qlen slen length qstart qend sstart send pident evalue bitscore" -out apple_self_megablast 

perl scripts/01parse_megablast.pl apple_phased_diploid_scaffolds.fasta apple_self_megablast

cat p1_secondary.fa p3_unmap.fa > apple_phased_diploid_scaffolds.alternative.fa

rm p1_secondary.fa p3_unmap.fa 

mv p2_primary.fa apple_phased_diploid_scaffolds.primary.fa

#step 2. assembly validation
Validation of the two separated assemblies (apple_phased_diploid_scaffolds.primary.fa and apple_phased_diploid_scaffolds.alternative.fa) was implemented using a pipeline same to that for the haploid assembly (see haploid pipeline for details).

