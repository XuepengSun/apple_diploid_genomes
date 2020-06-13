#!/usr/bin/perl -w 
use Bio::SeqIO;

my ($scaffolds,$unplaced) = @ARGV;

open SPLIT,$scaffolds || die;
while(<SPLIT>){
	if(/^>/){
		chomp;
		if(/\|/){
			$_=~/^>(\S+)\s+(\d+)/;
			push @{$gap{$1}},$2;
		}
	}
}
close SPLIT;

my $io=new Bio::SeqIO(-file=>$unplaced,-format=>'fasta');

while($seq=$io->next_seq){
	$data{$seq->id} = $seq->seq;
	my $sname = $seq->id;
	if($sname=~/\|/){
		$sname=~/(.*)\|(.*)/;
		push @{$join{$1}},$2;
	}
}

open OUT, ">unplaced_combined_list.txt";

foreach $s(keys %join){
	my @t = sort{$a<=>$b} @{$join{$s}};
	next if scalar @t < 2;
	push @t, 100000000;
	my %p = ();
	while(@t){
		my $x = shift @t;
		last if $x == 100000000;
		if($x + 1 == $t[0]){
			$p{$x}++;
			$p{$t[0]}++;
		}
		else{
			if(%p && scalar keys %p >=2){
				my @m = sort{$a<=>$b} keys %p;
				my $new_id = $s.'|'.$m[0].'-'.$m[-1];
				print OUT $new_id,"\n";
				my $tmp_seq = "";
				for $i(0..$#m-1){
					my $new_name = $s.'|'.$m[$i];
					my $sq = $data{$new_name}. 'N' x $gap{$new_name}[0];
					$tmp_seq.=$sq;
					delete($data{$new_name});
					delete($gap{$new_name});
				}
				my $last = $s.'|'.$m[-1];
				$tmp_seq.=$data{$last};
				$data{$new_id} = $tmp_seq;
				push @{$gap{$new_id}},$gap{$last}[0];
				delete($data{$last});
				delete($gap{$last});
				%p=();
			}
		}
	}
}

close OUT;

foreach(keys %data){
	print ">",$_;
	if($_=~/\|/){print "\t",$gap{$_}[0],"\n"}
	else{print "\n"}
	print $data{$_},"\n";
}









