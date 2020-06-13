#!/usr/bin/perl -w 
use Bio::SeqIO;
use Statistics::Basic qw(:all);

my ($fasta,$cov) = @ARGV;

my $flank = 2000;

open COV,$cov || die;
while(<COV>){
	chomp;
	my @s = split;
	for($s[1]..$s[2]-1){push @{$hash{$s[0]}},$s[3]}
}
close COV;

print STDERR "\ndone!\n";

my $io=new Bio::SeqIO(-file=>$fasta,-format=>'fasta');
while($seq=$io->next_seq){
	my $xl = $seq->seq;
	while($xl=~/[Nn]+/g){
		my $ltmp = length($`) - $flank;
		my $lst = $ltmp < 0 ? 0 : $ltmp;
		my $lend = length($`) - 1;
		my $gst = length($`) ;
		my $gend = length($`) + length($&) - 1 ;
		my $rst = length($`) + length($&);
		my $rtmp = length($`) + length($&) + $flank - 1 ;
		my $rend = $rtmp > $seq->length -1 ? $seq->length -1 : $rtmp ;
		my $lmedian = median(@{$hash{$seq->id}}[$lst..$lend]);
		my $gmedian = median(@{$hash{$seq->id}}[$gst..$gend]);
		my $rmedian = median(@{$hash{$seq->id}}[$rst..$rend]);
		print $seq->id,"\t",length($`),"\t",length($&),"\t",$lmedian,"\t",$gmedian,"\t",$rmedian,"\n";
	}
}


