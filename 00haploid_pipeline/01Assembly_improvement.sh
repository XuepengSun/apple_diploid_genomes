#step 1. de novo assembly of PacBio HiFi reads with HiCanu and Hifiasm
##assembly with Hifiasm 
hifiasm -o Gala -t 100 ccs_reads.fastq.gz

##assembly with Hicanu
canu genomeSize=750m -p Gala_canu -d data corOutCoverage=100 OvlMerThreshold=300 -pacbio-hifi ccs_reads.fastq.gz

##a primary assembly of each of the Hifiasm and Hicanu contigs was generated using Purge Haplotigs
pbmm2 align --sort -j 100 --preset HiFi -N 1 Gala.p_ctg.gfa.fasta CCS.fasta.fofn Gala.p_ctg.bam
purge_haplotigs hist -b Gala.p_ctg.bam -g Gala.p_ctg.gfa.fasta -t 50
purge_haplotigs contigcov -i aligned.bam.gencov -l 3 -m 26 -h 200
purge_haplotigs purge -g Gala.p_ctg.gfa.fasta -c coverage_stats.csv -a 70 -t 60

#step 2. gap closing based on whole genome alignment 
##anchor the primary assembly into pseudochromosomes using RaGOO
ragoo.py -t 90 hifiasm_primary_contigs.fasta reference.fasta

##align the Hifiasm chromosomes againt DenovoMagic assembly. Each chromosome was aligned separately
fasta_tool --split hifiasm.ragoo.fasta

#align sequences using Nucmer
for i in {1..17};do echo "nohup nucmer -t 5 -p chr$i chr${i}.fasta chr${i}_RaGOO.fasta &"; done > nucmer.sh; sh nucmer.sh

#filter the alignments using delta-filter
for i in {1..17};do echo "nohup delta-filter -g chr${i}.delta > chr${i}.g.delta &"; done > filter.sh; sh filter.sh

#output coordinates of filtered alignments
for i in {1..17};do show-coords -c -d -l -o -r -T chr${i}.g.delta > chr${i}.g.delta.coords; done

cat *.coords > aln.coords 

#prepare Hifiasm contig orders from the RaGOO output directory
perl -ne '$ARGV=~/(chr.*?)_orderings/;print $1,"\t",$_' orderings/chr* > hifiasm.ordering

#prepare DenovoMagic scaffold orders from the AGP file
perl -ne 'next if /^#/;@s=split;next unless $s[0]=~/chr/i && $s[5] ne "map";$s[8]="+" if $s[8] eq "?";print $s[0],"\t",$s[5],"\t",$s[8],"\n"' Gala.agp >  Gala.ordering

#gap closing 
perl -i -pe 's/_RaGOO//' hifiasm.ragoo.fasta

perl scripts/GapCloser_with_mummer.pl \
--dscaf hifiasm_primary_contigs.fasta \
--dchr hifiasm.ragoo.fasta \
--dord hifiasm.ordering \
--rscaf Gala.scaffold.fasta \
--rchr Gala.chr.fasta \
--rord Gala.ordering \
--minalnlen 10000 \
--coords aln.coords 

#the output files include improved scaffold assembly ("new_scaffold.fa") and chromosome assembly ("new_chromosome.fa"), these can be used as input for second round improvment with Hicanu assembly using the same pipeline

#in the final improved assembly, a self all-vs-all blast was performed to remove redundant scaffolds

