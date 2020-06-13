########################################################################
##	 				 	diploid assembly improvement	
########################################################################

#Note: assembly of PacBio HiFi reads was described in the haploid assembly improvement pipeline. The unitig assembly was used for diploid genome improvement 

#step 1. align hifiasm unitig assembly to the diploid scaffolds using Nucmer
fasta_tool --chunk 20 apple.p_utg.fasta

for i in *.fasta;do echo "nohup nucmer -t 5 -p $i apple_phased_diploid_scaffolds.fa $i &";done > nucmer.sh; sh nucmer.sh

cat *.delta | grep -v '^NUCMER'  > aln.delta # first two lines of aln.delta need to modified for downstream analysis

delta-filter -q -r aln.delta > aln.qr.delta

delta-filter -g aln.qr.delta > aln.qrg.delta

show-coords -c -d -l -o -r -T aln.qrg.delta > aln.qrg.delta.coords

#step 2. gap closing
perl scripts/GapCloser_diploid_assembly_Mummer.pl \
-r apple_phased_diploid_scaffolds.fa \
-d apple.p_utg.fasta \
-a aln.qrg.delta.coords \
-c 0.7 --overlap 3000 \
--maxdist 3000 --seedsize 5000

#step 3. second round improvement 
The output sequence file from step 2 was used as input for second round improvment using HiCanu unitig assembly following the same pipeline described above

#step 4. remove redundancy
In the final improved assembly, a self all-vs-all blast was performed to remove redundant scaffolds

