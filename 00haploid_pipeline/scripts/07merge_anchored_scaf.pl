#!/usr/bin/perl -w 
use Bio::SeqIO;

my ($chr_agp,$genome) = @ARGV;

open AGP,$chr_agp || die;
while(<AGP>){
	next if /^#/;
	chomp;
	my @s = split;
	if($s[4] eq 'W'){
		push @{$hash{$s[0]}},$s[5];
		push @{$strand{$s[0]}},$s[8];
		$anchor{$s[5]} = 1;
	}
}
close AGP;

open IN,$genome || die;
while(<IN>){
	if(/^>/){
		chomp;
		$_=~s/^>//;
		@t = split(/\s+/,$_);
		if(scalar @t == 2){
			($id,$gap) = @t;
			push @{$gap_size{$id}},$gap;
		}
		else{$id = $t[0]}
	}
	else{
		chomp;
		$seq{$id}=$_;
	}
}
close IN;

open OUT,">join_list.txt";

foreach $lg(keys %hash){
	my @t = @{$hash{$lg}};
	my @strd = @{$strand{$lg}};
	my %s = ();
	my %strand_tag = ();
	my $name = "";
	my $nid = 0;
	while(@t){
		my $v = shift @t;
		my $orient = shift @strd;
		$strand_tag{$nid}{$orient} = 1;
		if(my ($scaf,$num)=$v=~/(.*)\|(.*)/){
			if($scaf ne $name || ($strand_tag{$nid}{'+'} && $strand_tag{$nid}{'-'})){$nid++}
			$s{$nid}{'scaf'} = $scaf;
			$s{$nid}{'split'}{$num} = 1;
			$strand_tag{$nid}{$orient} = 1;
			$name = $scaf;
		}
		else{$nid++}
	}
	foreach $i(keys %s){
		next if scalar keys %{$s{$i}{'split'}} < 2;
		my $sid = $s{$i}{'scaf'};
		my @range = sort{$a<=>$b} keys %{$s{$i}{'split'}};
		my %u = ();
		my $tpid = 0;
		for $j($range[0]..$range[-1]){
			my $tmp_scaf = $sid.'|'.$j;
			if(!$s{$i}{'split'}{$j} && $anchor{$tmp_scaf}){$tpid++}
			else{
				$u{$tpid}{$j} = 1;
			}
		}
		foreach $m(keys %u){
			my @nn = sort{$a<=>$b} keys %{$u{$m}};
			next if scalar @nn < 2;
			my $new_id = $sid."|".$nn[0].'-'.$nn[-1];
			print OUT $new_id,"\n";
			my $new_seq = "";
			for $dx(0..$#nn-1){
				$old_id = $sid.'|'.$nn[$dx];
				my $tp_seq = $seq{$old_id} . 'N' x $gap_size{$old_id}[0];
				$new_seq.=$tp_seq;
				delete($seq{$old_id});
				delete($gap_size{$old_id});
			}
			my $last_id = $sid.'|'.$nn[-1];
			$seq{$new_id}=$new_seq.$seq{$last_id};
			push @{$gap_size{$new_id}},$gap_size{$last_id}[0];
			delete($seq{$last_id});
			delete($gap_size{$last_id});
		}
	}
}
close OUT;

foreach(keys %seq){
	print ">",$_;
	if($gap_size{$_}){print "\t",$gap_size{$_}[0],"\n"}
	else{print "\n"}
	print $seq{$_},"\n";
}










