#!/usr/bin/perl -w
use Bio::SeqIO;

my ($fasta,$dep) = @ARGV;

my $flank = 2000;
my $minDepth = 20;

my $io = new Bio::SeqIO(-file=>$fasta,-format=>'fasta');
while($seq=$io->next_seq){$hash{$seq->id}=$seq->seq}
$io->close;


open IN,$dep || die;
while(<IN>){
	chomp;
	my @s = split;
	$s[3]=~s/,//g;
	$s[4]=~s/,//g;
	$s[5]=~s/,//g;
	my $mean = ($s[3]+$s[5]) / 2;
	if($s[4] < $mean / 2 || $s[4] < $minDepth){
		my $ltmp = $s[1] - $flank;
		my $lst = $ltmp < 0 ? 0 : $ltmp;
		my $llen = $s[1] - $lst;
		my $rst = $s[1] + $s[2];
		my $rtmp = $rst + $flank ;
		my $rend = $rtmp > length($hash{$s[0]}) ? length($hash{$s[0]}) : $rtmp ;
		my $rlen = $rend - $rst;
		my $lseq = substr($hash{$s[0]},$lst,$llen);
		my $rseq = substr($hash{$s[0]},$rst,$rlen);
		my @ltp = split(/[Nn]+/,$lseq);	
		my @rtp = split(/[Nn]+/,$rseq);
	
		print '>',$s[0].'_'.$s[1].'_'.$s[2].'_L',"\n",$ltp[-1],"\n";
		print '>',$s[0].'_'.$s[1].'_'.$s[2].'_R',"\n",$rtp[0],"\n";
	}
}
close IN;

