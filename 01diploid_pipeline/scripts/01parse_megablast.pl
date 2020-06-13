#!/usr/bin/perl -w
use Bio::SeqIO;
use List::Util qw(max min);

my $identity = 95;
my $coverage = 0.5;
my $genome = $ARGV[0];
my $blast = $ARGV[1];

my $io = new Bio::SeqIO(-file=>$genome,-format=>'fasta');
while($seq = $io->next_seq){$res{$seq->id}=$seq->seq}

open BLAST,$blast || die;
while(<BLAST>){
			my @s = split;
			next if ($s[0] eq $s[1]) ;
			next if ($s[9] < $identity);
			push @{$hash{$s[0]}{$s[1]}{$s[5]}},$s[6];
}
close BLAST;

print STDERR "\nreading blast done\n";

open A,'>p1_secondary.fa';
open B,'>p2_primary.fa';
open C,'>p3_unmap.fa';

foreach $i(keys %hash){
		my $name = "";
		my $longest = 0;
		foreach $j(keys %{$hash{$i}}){
			my @pos = ();
			foreach $m(sort{$a<=>$b} keys %{$hash{$i}{$j}}){
				my $n = max @{$hash{$i}{$j}{$m}};
				if(@pos && $m <= $pos[-1]){
					$pos[-1] = max($pos[-1],$n);
				}
				else{push @pos,($m,$n)}
			}
			my $len = 0;
			while(@pos){
				my $x = shift @pos;
				my $y = shift @pos;
				$len += $y - $x + 1;
			}
			if($len > $longest){
				$longest = $len;
				$name = $j;
			}
		}
		if($longest / length($res{$i}) >= $coverage){
			print A '>',$i,"\n",$res{$i},"\n";
			$final{$i}++;
		}
		else{print B '>',$i,"\n",$res{$i},"\n";$final{$i}++}
}

foreach(keys %res){
	print C '>',$_,"\n",$res{$_},"\n" if !$final{$_};
}

close A;
close B;
close C;
