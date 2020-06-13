#!/usr/bin/perl -w
use Bio::SeqIO;

my ($genome,$range,$fragcov) = @ARGV;

my (%cov);

open IN,$range || die;
while(<IN>){
	chomp;
	my @s = split;
	$hash{$s[0]}{$s[1]} = $s[2];
}
close IN;

open IN,$fragcov || die;
while(<IN>){
	chomp;
	my @s = split;
	$cov{$s[0]}{$s[1]} = $_;
}
close IN;


my $io=new Bio::SeqIO(-file=>$genome,-format=>'fasta');
while($seq=$io->next_seq){
	next if !$hash{$seq->id};
	my $xl = $seq->seq;
	while($xl=~/[Nn]+/g){
		my $st = length($`) + 1;
		my $end = length($`) + length ($&);
		foreach $i(keys %{$hash{$seq->id}}){
			if($st >= $i && $end <= $hash{$seq->id}{$i}){
				my $x = length($`);
				my $c = "na";
				if($cov{$seq->id}{$x}){
					$c=$cov{$seq->id}{$x};
				}
				print $seq->id,"\t",$seq->id,"\t",length($`),"\t",length($&),"\t",$c,"\n";
			}
		}
	}
}


