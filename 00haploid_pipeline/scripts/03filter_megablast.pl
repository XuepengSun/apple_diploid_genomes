#!/usr/bin/perl -w

my ($gap_cov,$blast) = @ARGV;

my $thres = 0.5 ;
my $cov_thres = 20;

open BLAST, $blast || die;
while(<BLAST>){
	chomp;
	my @s = split;
	next if $ck{$s[0]};
	$ck{$s[0]} = 1;
	next if $s[4] / $s[2] < $thres;
	$s[0]=~/(.*)_(.)/;
	$data{$1}{$2} = \@s;		
}
close BLAST;

open IN,$gap_cov || die;
while(<IN>){
	chomp;
	my @s = split;
    $s[3]=~s/,//g;
    $s[4]=~s/,//g;
    $s[5]=~s/,//g;
    my $mean = ($s[3]+$s[5]) / 2;
    if($s[4] < $mean / 2 || $s[4] < $cov_thres){
		my $id = $s[0].'_'.$s[1].'_'.$s[2];
		if(!$data{$id}){print "none","\t",$_,"\n"}
		else{
			if(!$data{$id}{'L'} || !$data{$id}{'R'}){print "single\t",$_,"\n"}
			else{
				if($data{$id}{'L'}->[1] ne $data{$id}{'R'}->[1]){print "translocation\t",$_,"\n"}
				elsif(($data{$id}{'L'}->[7]-$data{$id}{'L'}->[8])*($data{$id}{'R'}->[7]-$data{$id}{'R'}->[8]) < 0){print "inversion\t",$_,"\n"}	
				elsif(abs($data{$id}{'L'}->[7] - $data{$id}{'R'}->[8]) > 100000){
					print "toolong\t",$_,"\n";
				}
				else{next}	
			}	
		}
	}
}
close IN;


