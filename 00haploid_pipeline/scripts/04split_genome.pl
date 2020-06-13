#!/usr/bin/perl -w
use Bio::SeqIO;

my ($genome,$gapfile) = @ARGV;

my %hash;
open GAP,$gapfile || die;
while(<GAP>){
	chomp;
	my @s = split /\t/;
	$hash{$s[1]}{$s[2]} = $s[3];
}
close GAP;

my $io = new Bio::SeqIO(-file=>$genome,-format=>'fasta');
while($seq=$io->next_seq){
	if(!$hash{$seq->id}){print '>',$seq->id,"\n",$seq->seq,"\n"}
	else{
		my $st = 1;
		my $n = 1;
		foreach $i(sort{$a<=>$b} keys %{$hash{$seq->id}}){
		  my $end = $i;
		  my $sq = $seq->subseq($st,$end);
		  print '>',$seq->id,"|",$n,"\t",$hash{$seq->id}{$i},"\n",$sq,"\n";
		  $st = $i + $hash{$seq->id}{$i}+1 ;
		  $n++;
		}
		if($st <= $seq->length){
			$sq = $seq->subseq($st,$seq->length);
		 	print '>',$seq->id,"|",$n,"\t",0,"\n",$sq,"\n";
		}
	}
}
$io->close;

