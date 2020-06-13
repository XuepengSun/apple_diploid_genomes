#############################################################################################
#					pipeline for haploid assembly validation
#############################################################################################

#step 1. map mate-pair (MP) reads against scaffold assembly
Note: the cleaned and deduplicated reads were mapped with "bwa mem" and the properly paired alignments were sorted and deduplicated. Each MP library with different insertion size was mapped separatedly and then combined to a final output file named "mp.bam".

#step 2. calculate fragment coverage 
bedtools genomecov -bga -pc -ibam mp.bam > mp.frag.cov

#step 3. calculate fragment coverage on gap regions as well as 2kb regions flanking the gaps
perl scripts/01depth_span_gaps.pl scaffolds.fa mp.frag.cov > gap_depth.txt

#step 4. get 2kb flanking sequences of gaps whose coverage is lower than 20 or half of the 2kb flanking regions
perl scripts/02filter_depth_span_gaps.pl scaffolds.fa gap_depth.txt > flanking2kb.fa

#step 5. blastn search of the 2kb sequences against GDDH13 genome 
makeblastdb -in GDDH13_1-1_formatted.fasta -dbtype nucl

blastn -task megablast -word_size 30 -db GDDH13_1-1_formatted.fasta -num_threads 40 -evalue 1e-5 -max_target_seqs 20 -outfmt "6 qseqid sseqid qlen slen length qstart qend sstart send pident evalue bitscore" -query flanking2kb.fa -out flanking2kb.bn

#step 6. record suspicious gaps based on fragment coverage and blastn
perl scripts/03filter_megablast.pl gap_depth.txt flanking2kb.bn > bn_susp_gaps.txt

#step 7. split scaffolds at suspicious gaps
perl scripts/04split_genome.pl scaffolds.fa bn_susp_gaps.txt > scaffolds.split.fa

#step 8. validate scaffolds using the high-density genetic map (20k SNPs data)
bwa mem -t 60 scaffolds.split.fa ../resource/20k.snp.fa > 20k.snp.sam 

perl scripts/20k_processing.pl scaffolds.split.fa ../resource/20K_SNP_MARKERS.MAP.dat ../resource/20K_SNP_MARKERS.dat 20k.snp.sam ../resource/20k_map_dist.txt > 20k.map.csv

## (optional) get synteny between the scaffold assembly and GDDH13 genome, requires GeneMark-ES, diamond and ALLMAPS. The output is "synteny.bed". 
gmes_petap.pl --ES --min_contig 2000 --cores 40 --sequence scaffolds.split.fa

perl -ne '$_=~/(^\S+).*(GeneMark.*)/;print $1,"\t",$2,"\n"' genemark.gtf > genemark.m.gtf

perl scripts/GetSeqFromGTF.pl scaffolds.split.fa genemark.m.gtf

perl scripts/gtf2gff3.pl genemark.m.gtf > genemark.gff
	 
diamond makedb --in GDDH13_1-1_prot.fasta --db GDDH13_1-1_prot.fasta

diamond makedb --in output.faa --db output.faa

diamond blastp -d GDDH13_1-1_prot.fasta.dmnd -q output.faa -p 50 -k 1 -e 1e-10 --more-sensitive -f 6 -o scaf_to_ref.bp
	 
diamond blastp -d output.faa.dmnd -q GDDH13_1-1_prot.fasta -p 50 -k 1 -e 1e-10 --more-sensitive -f 6 -o ref_to_scaf.blastp 
	 
cat *.blastp | perl -ne '{chomp;@s=split;@t=($s[0],$s[1]);next if $ck{$s[0]};$ck{$s[0]}=1;@t=sort @t;$j=join("|",@t);push @{$hash{$j}},$s[-1]}END{foreach $i(keys %hash){if(scalar @{$hash{$i}} ==2){($a,$b)=$i=~/(.*)\|(.*)/;print $a,"\t",$b,"\t",$hash{$i}->[0],"\n"}}}' > ortholog.txt
	
perl -ne 'chomp;@s=split;if($s[2] eq "gene"){$s[8]=~/ID=(.*?);/;print $s[0],"\t",$s[3],"\t",$s[4],"\t",$1,"\t100\t",$s[6],"\n"}' genemark.gff > genemark.bed

# gene_models_20170612.gff3 is the annotation file of GDDH13 genome
perl -ne 'chomp;@s=split;if($s[2] eq "gene"){$s[8]=~/ID=gene:(.*?);/;print $s[0],"\t",$s[3],"\t",$s[4],"\t",$1,"\t100\t",$s[6],"\n"}' gene_models_20170612.gff3 > ref.bed 
	 
python -m jcvi.assembly.syntenypath bed ortholog.txt --qbed genemark.bed --sbed ref.bed --scale=100000 -o synteny.bed

#note: if synteny step is skipped, replace "synteny.bed" with '20k.map.csv' in the following script
#parameter marker_thres can be adjused to get the optimal value for breakpoint selection
perl scripts/05check_scaf_with_marker.pl 20k.map.csv synteny.bed > position_to_split.txt

#manually examine the position_to_split.txt file carefully, and create a 'split.local.txt' file in a format like "scaffoldID\tstartPosition\tendPostion" to indicate highly suspicious region to break.
perl scripts/06get_gap_and_fragCov.pl scaffolds.split.fa split.local.txt gap_depth.txt > split.local.gaps.txt
  
#the split.local.gaps.txt includes gaps within the suspicious regions identified based on the genetic map and optionally genomic synteny. It is highly recommened to make some selection on the gaps (e.g. coverage < 20 and others based on manual examination) before going ahead to second round scaffold split.

perl scripts/04split_genome.pl scaffolds.split.fa split.local.gaps.txt > scaffolds.split2.fa

#step 9. anchor scaffolds with genetic maps and genomic synteny using ALLMAPS
bwa mem -t 60 scaffolds.split2.fa ../resource/8K_SNP.fa > 8k.snp.sam 
bwa mem -t 20 scaffolds.split2.fa ../resource/20k.snp.fa > 20k.snp.sam 

perl scripts/20k_processing.pl scaffolds.split2.fa ../resource/20K_SNP_MARKERS.MAP.dat ../resource/20K_SNP_MARKERS.dat 20k.snp.sam ../resource/20k_map_dist.txt > 20k.map.csv

perl scripts/8k_processing.pl scaffolds.split2.fa ../resource/8K_SNP_LineageGroup.dat 8k.snp.sam ../resource/8k_map_dist.txt > 8k.snp.csv
	
python -m jcvi.assembly.allmaps merge 20k.map.csv 8k.snp.csv -o GeneticMap.bed
	
perl -ne '$_=~s/8k-//;$_=~s/20k-//;print $_' GeneticMap.bed > GeneticMap.bed2; mv GeneticMap.bed2 GeneticMap.bed
	
# combined genetic map and syntenic map
python -m jcvi.assembly.allmaps mergebed GeneticMap.bed synteny.bed -o all.maps.bed
	
python -m jcvi.assembly.allmaps path all.maps.bed scaffolds.split2.fa

#step 10. merge split scaffolds if they are anchored close togther
perl scripts/07merge_anchored_scaf.pl all.maps.chr.agp scaffolds.split2.fa > scaffolds.split2.m.fa

#step 11. merge unplaced scaffolds and output the final scaffold sequences
perl scripts/08merge_unplaced_scaf.pl scaffolds.split2.m.fa all.maps.unplaced.fasta > scaffolds.final.fa
	
