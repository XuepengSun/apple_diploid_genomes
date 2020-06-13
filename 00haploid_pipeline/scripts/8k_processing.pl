#!/usr/bin/perl -w
use Bio::SeqIO;

my ($genome,$data,$sam,$dist) = @ARGV;

my $io=new Bio::SeqIO(-file=>$genome,-format=>'fasta');
while(my $seq=$io->next_seq){
	$seqs{$seq->id}=$seq->length;
}
$io->close;

open IN,$dist;
while(<IN>){
	next if /^\s+/;
	chomp;
	my @s = split;
	$distance{$s[0]}=$s[1];
}
close IN;

open IN,$data;
while(<IN>){
	next if $.==1;
	my @s = split /\t/;
	$hash{$s[0]} = {'LG'=>"LG$s[1]","dist"=>$s[2]};
}
close IN;

print "Scaffold ID,scaffold position,LG,genetic position\n";

open SAM,$sam;
while(<SAM>){
	next if /^\@/;
	my @s = split;
	next if $s[4] < 20;
	if($s[1] == 0){
		$pos = $s[3] + length($s[9]);
	}
	elsif($s[1]==16){
		$pos = $s[3] - 1;
	}
	#print $s[2],"\t",$seqs{$s[2]},"\t",$pos,"\t",$hash{$s[0]}{'LG'},"\t",$distance{$hash{$s[0]}{'LG'}},"\t",$hash{$s[0]}{'dist'},"\n";
	print $s[2],",",$pos,",",$hash{$s[0]}{'LG'},",",$hash{$s[0]}{'dist'},"\n";
}
close SAM;









