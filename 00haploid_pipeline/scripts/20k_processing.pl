#!/usr/bin/perl -w
use Bio::SeqIO;

my ($genome,$infor,$data,$sam,$dist) = @ARGV;

my $io=new Bio::SeqIO(-file=>$genome,-format=>'fasta');
while(my $seq=$io->next_seq){
	$seqs{$seq->id}=$seq->length;
}
$io->close;

open IN,$infor || die;
while(<IN>){
	next if $. == 1;
	chomp;
	my @s = split /\t/;
	next if $s[1] ne "";
	next if $s[-1] !~/\d+/;
	$valid{$s[0]}={'LG'=>"LG$s[2]",'dist'=>$s[-1]}
}
close IN;

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
	if($s[5] eq $s[6]){$strand{$s[0]} = 1}
	else{$strand{$s[0]} = -1}
}
close IN;

print "Scaffold ID,scaffold position,LG,genetic position\n";

open SAM,$sam;
while(<SAM>){
	next if /^\@/;
	my @s = split;
	next if $s[4] < 20;
	next if !$valid{$s[0]};
	my $strd = $s[1] == 0 ? 1 : -1;
	if($strand{$s[0]} * $strd > 0){
		$pos = $s[3] + length($s[9]);
	}
	else{
		$pos = $s[3] - 1;
	}
	#print $s[2],"\t",$seqs{$s[2]},"\t",$pos,"\t",$valid{$s[0]}{'LG'},"\t",$distance{$valid{$s[0]}{'LG'}},"\t",$valid{$s[0]}{'dist'},"\n";
	print $s[2],",",$pos,",",$valid{$s[0]}{'LG'},",",$valid{$s[0]}{'dist'},"\n";
}
close SAM;









