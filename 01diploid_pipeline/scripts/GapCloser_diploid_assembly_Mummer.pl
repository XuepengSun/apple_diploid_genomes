#!/usr/bin/perl
use warnings;
use Term::ANSIColor;
use Bio::SeqIO;
use Getopt::Long::Descriptive;
use Bio::SeqIO;
use threads;
use File::Find;
use List::Util qw(sum min max);

my ($opt, $usage) = describe_options(
  'Usage: %c %o <args>',
	[],
  ['rscaf|r=s',"reference scaffold sequence in fasta format", {required => 1}],
	['dscaf|d=s',"query scaffold sequence in fasta format",{required => 1}],
	['alignment|a=s',"alignment between reference and query scaffolds in coords format",{required => 1}],
	['mincov|c=f',"minimal coverage of unitigs or scaffolds, default=0.7",{default => 0.7}],
	['overlap|o=i',"maximal overlap length between alignments, default=1000",{default=>1000}],
	['maxdist|m=i',"maximal distance to scaffold between two seqeunces, default=1000",{default=>1000}],
	['break|b=i',"stop scaffolding in inconsecutive region, yes (1,default) or no (0)",{default=>1}],
	['segratio|s=f',"overlap alignment ratio to stop scaffolding, default=0.8",{default => 0.8}],
	['identity|i=f',"minimal identity for anchors, default=99.5",{default => 99.5}],
	['maxchunk|k=i',"maximal chunk size to replace, default=1000000",{default=>1000000}],
	['seedsize|e=i',"minimal anchor alignment size, default=10000",{default=>10000}],
	['prefix|p=s',"prefix of output, default=tmp_dip",{default => "tmp_dip"}],
  ['help|h', "print usage message and exit", {shortcircuit => 1}],
	['version|v',"version 0.2",{shortcircuit => 1}],
  []
);
print($usage->text), exit if $opt->help;
print STDERR "\n";

my $subFile = $opt->prefix.'.substitutions';
my $scafFile = $opt->prefix.'.scaffolding';
my $delFile = $opt->prefix.'.deletions';
my $table = $opt->prefix.'.nameTable.txt';
my $fasta = $opt->prefix.'.sequences.fasta';

LOG("reading ".$opt->rscaf);
my $rscaf = &read_seq($opt->rscaf);

LOG("reading ".$opt->dscaf);
my $dscaf = &read_seq($opt->dscaf);

LOG("mapping gaps in ".$opt->rscaf);
my $gaps = &locate_gaps($opt->rscaf);

LOG("starting to process ".$opt->alignment);
my ($coords,$strand) = &read_coords($opt->alignment,$opt->identity);

LOG("finding substitutions");
my ($sub,$redundant,$bound) = &substitution($coords,$strand,$subFile,$delFile);

LOG("finding connections");
&connection($sub,$coords,$strand,$bound,$scafFile); 

LOG("Finalizing");
combine($subFile,$scafFile,$fasta,$table);



sub combine{
	my ($subs,$conns,$fasta,$table) = @_;
	my (%bridge,%dup,%info,$add,$del,%exclude);
	open CON,$conns || die;
	while(<CON>){
		next if $. == 1;
		my @s = split;
		push @{$info{$s[0]}{$s[1]}},$_;
		$dup{$s[2]}++;
	}
	close CON;
	foreach $u(keys %info){
		foreach $i(keys %{$info{$u}}){
			my $ck = 0;
			my @pos = ();
			my %list = ();
			$delTMP = 0;
			foreach(@{$info{$u}{$i}}){
				my @k = split;
				$ck++ if $dup{$k[2]} > 1;
				$ck++ if ($k[6]-$k[5]+1)/length($rscaf->{$k[2]}) < 0.95;
				push @pos,($k[3],$k[4]);
				$list{$k[2]} =1;
				$delTMP += length($rscaf->{$k[2]});
			}
			next if $ck > 0;
			$del+=$delTMP;
			@pos = sort{$a<=>$b} @pos;
			$add += $pos[-1] - $pos[0] + 1;
			my $id = $u.'_'.$i;
			$bridge{$id} = substr($dscaf->{$u},$pos[0]-1,$pos[-1]-$pos[0]+1);
			foreach(keys %list){$exclude{$_}=$id}
		}
	}
	LOG("scaffolding: add $add bp; delete $del bp");
	
	my %data = ();
	my %repSeq = ();
	open SUB,"sort -k1,1 -k2n,2 $subs|" || die;
	while(<SUB>){
		next if $. == 1;
		my @s = split;
		next if $exclude{$s[0]};
		next if($s[2]-$s[1]+1<$opt->seedsize || $s[4]-$s[3]+1<$opt->seedsize);
		push @{$data{$s[0]}},$_;
	}
	close SUB;
	foreach $s(keys %data){
		$exclude{$s} = $s;
		my @aln = @{$data{$s}};
		for($i=0;$i<=$#aln;$i++){
			my @k = split(/\s+/,$aln[$i]);
			my $g = 0;
			foreach $m(sort{$a<=>$b} keys %{$gaps->{$s}}){
				last if $m > $k[2];
				next if $gaps->{$s}{$m} < $k[1];
				$g++;
			}
			if($g>0){
				$tmpseq=substr($dscaf->{$k[5]},$k[3]-1,$k[4]-$k[3]+1);
				if($k[6] eq '-1'){
					$tmpseq = join("",reverse split("",$tmpseq));
					$tmpseq=~tr/ATGCatgc/TACGtacg/;
				}
				$repSeq{$s} .= $tmpseq;
				$add+=$k[4]-$k[3]+1;
				$del+=$k[2]-$k[1]+1;
			}
			else{$repSeq{$s}.=substr($rscaf->{$s},$k[1]-1,$k[2]-$k[1]+1)}
			if($i < $#aln){
				my @p = split(/\s+/,$aln[$i+1]);
				next unless $p[1] > $k[2] + 1;
				$rst = $k[2] + 1;
				$rend = $p[1] - 1;
				my $rfseq = substr($rscaf->{$s},$k[2],$rend-$rst+1);
				if($k[5] ne $p[5]){
					$repSeq{$s} .= $rfseq;
				}
				else{
					$g = 0;
					foreach $m(sort{$a<=>$b} keys %{$gaps->{$s}}){
						last if $m > $rend;
						next if $gaps->{$s}{$m} < $rst;
						$g++;
					}
					if($g > 0){
						my @nn = sort{$a<=>$b}($k[3],$k[4],$p[3],$p[4]);
						my $dist = ($k[4]-$k[3]+1)+($p[4]-$p[3]+1)-($nn[-1]-$nn[0]+1);
						if($dist >= 0){
							$del += $rend - $rst + 1;
							next;
						}
						elsif($dist + $opt->maxchunk < 0){
							$repSeq{$s} .= $rfseq;
						}
						else{
							$dst = $nn[1]+1;
							$dend = $nn[2]-1;
							$tmpseq = substr($dscaf->{$k[5]},$dst-1,$dend-$dst+1);
							if($k[6] eq '-1'){
								$tmpseq = join("",reverse split("",$tmpseq));
								$tmpseq=~tr/ATGCatgc/TACGtacg/;
							}
							$repSeq{$s} .= $tmpseq;
							$add += $dend - $dst + 1;
							$del += $rend - $rst + 1;
						}
					}
					else{
						$repSeq{$s} .= $rfseq;
					}
				}
			}
		}
		my @h = split(/\s+/,$aln[0]);
		if($h[1]>1){
			$repSeq{$s} = substr($rscaf->{$s},0,$h[1]-1).$repSeq{$s};
		}
		my @b = split(/\s+/,$aln[-1]);
		if($b[2] < length($rscaf->{$s})){
			$repSeq{$s} = $repSeq{$s} . substr($rscaf->{$s},$b[2]);
		}
	}
	open OUT,">$fasta";
	open TMP,">comparisonTMP.txt";
	foreach $i(keys %{$rscaf}){
		next if $exclude{$i};
		print OUT '>',$i,"\n",$rscaf->{$i},"\n";
	}
	foreach(keys %repSeq){
		print OUT '>',$_,"\n",$repSeq{$_},"\n";
		print TMP $_,"\t",length($rscaf->{$_}),"\t",length($repSeq{$_}),"\n";
	}
	foreach(keys %bridge){print OUT '>',$_,"\n",$bridge{$_},"\n"}
	close OUT;
	close TMP;

	open EX,">$table" ;
	foreach(keys %exclude){print EX $_,"\t",$exclude{$_},"\n" if $_ ne $exclude{$_}	}
	close EX;

	LOG("final added sequences: $add bp");
	LOG("final deleted sequences: $del bp");
}


sub connection{
	my $replace = shift;
	my $data = shift ;
	my $strand = shift;
	my $bound = shift;
	my $output = shift;
	open OUT,">$output";
	print OUT "unitig\tgroup\tscaffold\tUnitigSt\tUnitigEnd\tScafSt\tScafEnd\tStrand\n";
	my (%valid,%large,%range,%margin);
	foreach $u(keys %$replace){
		my @kk = @{$replace->{$u}};
		while(@tt=splice(@kk,0,7)){
			$valid{$tt[0]}{$tt[-2]}=1;
			push @{$range{$u}},@tt[3..4];
		}
	}
	my %lostStrand = ();
	foreach $s(keys %$data){
		foreach $u(keys %{$data->{$s}}){
			my @pos = ();
			foreach(@{$data->{$s}{$u}}){
				my @k = split;
				my ($st,$end) = $k[2] < $k[3]?($k[2],$k[3]):($k[3],$k[2]);
				push @pos,($k[0],$k[1],$st,$end);
				$lostStrand{$u}{$s}{$k[12]} = 1;
			}	
			my ($tt,$rs,$qs,$left,$right) = &get_range(\@pos);
			my @k = @$tt;
			my ($rleft,$rright) = (100000000000000,0);
			while(@pp = splice(@k,0,4)){
				$rleft = $rleft < $pp[0] ? $rleft : $pp[0];
				$rright = $rright > $pp[1] ? $rright : $pp[1];
			}									  
			$large{$u}{$s}=[$left,$right,$rleft,$rright];
		}
	}
	foreach $u(keys %large){
		next unless $replace->{$u};
		my %new = ();
		my %order = ();
		my $num =0;
		my @occupy = sort{$a<=>$b} @{$range{$u}};
		foreach $s(keys %{$large{$u}}){
			$new{$large{$u}{$s}->[0]} = $s;
		}
		my @scafs = map{$new{$_}} sort{$a<=>$b} keys %new;
		next if scalar @scafs < 2;
		foreach $s(@scafs){
			if($valid{$s}{$u}){
				## gap size control - two scaffolds should locate within 1kb
				if($order{$num}){
					my $dist = $large{$u}{$s}->[0] - $large{$u}{$order{$num}[-1]}->[1];
					my $op = $large{$u}{$order{$num}[-1]}->[1] - $large{$u}{$s}->[0];
					my $sp1 = $large{$u}{$order{$num}[-1]}->[1] - $large{$u}{$order{$num}[-1]}->[0] + 1;
					my $sp2 = $large{$u}{$s}->[1] - $large{$u}{$s}->[0] + 1;
					my $rd1=$strand->{$order{$num}[-1]}{$u} eq '1' ? length($rscaf->{$order{$num}[-1]})-$large{$u}{$order{$num}[-1]}->[1] : $large{$u}{$order{$num}[-1]}->[0];
					my $rd2 = $strand->{$s}{$u} eq '1' ? $large{$u}{$s}->[0] : length($rscaf->{$s}) - $large{$u}{$s}->[1]; 
					if($dist>$opt->maxdist || $op/$sp1 >$opt->segratio || $op/$sp2>$opt->segratio || ($rd1 + $rd2)>$opt->maxdist * 2){
						$num++;
					}
				}
				push @{$order{$num}},$s;
			}
			else{
				($qst,$qend,$rst,$rend) = @{$large{$u}{$s}};
				if((($qst<100 && $qend-$occupy[0]<$opt->overlap && $occupy[0]-$qend <=$opt->maxdist) || ($occupy[-1]-$qst<$opt->overlap && length($dscaf->{$u})-$qend<100 && $qst-$occupy[-1]<=$opt->maxdist)) && ($rst<100 || length($rscaf->{$s})-$rend<100)){
					my @str = keys %{$lostStrand{$u}{$s}};
					if(scalar @str == 1){
						push @{$order{$num}},$s;
					}
					else{if($opt->break == 1){$num++}}
				}
				else{
					if($opt->break == 1){$num++}
				}
			}
		}					
		foreach $i(keys %order){
			next if scalar @{$order{$i}} < 2;
			foreach(@{$order{$i}}){
				if($strand->{$_}{$u}){$ss = $strand->{$_}{$u}}
				else{
					@sp = keys %{$lostStrand{$u}{$_}};  ### order
					$ss = $sp[0];
				}
				print OUT $u,"\t",$i,"\t",$_,"\t",join("\t",@{$large{$u}{$_}}),"\t",$ss,"\n";
			}
		}
	}
	close OUT;
}

sub substitution{
	my $data = shift ;
	my $strand = shift;
	my $output = shift;
	my $delFile = shift;
	my (%clean,%bound,%size,%array,%remove);
	foreach $s(keys %$data){
		foreach $u(keys %{$data->{$s}}){
			my @tmp = ();
			foreach(@{$data->{$s}{$u}}){
				my @t = split;
				next if $t[12] ne $strand->{$s}{$u};
				my ($st,$end) = $t[2]<$t[3]?($t[2],$t[3]):($t[3],$t[2]);
				push @tmp,($t[0],$t[1],$st,$end);
			}
			my ($tt,$rs,$qs,$left,$right) = &get_range(\@tmp);
			next unless $rs/length($rscaf->{$s})>=$opt->mincov || $qs/length($dscaf->{$u})>=$opt->mincov;
			$array{$u}{$s} = $tt;
			$size{$u}{$s} = $qs;
			$bound{$u}{$s} = [$left,$right];
		}
	}
	foreach $u(keys %array){
			my @range = ();
			my %candidate = ();
			foreach $s(sort{$size{$u}{$b}<=>$size{$u}{$a}} keys %{$size{$u}}){
					my @tt = @{$array{$u}{$s}};
					if(!@range){
						$candidate{$s} = 1;
						while(@w=splice(@tt,0,4)){
							push @range,($s,@w,$u,$strand->{$s}{$u});
						}
					}
					else{
						my $skip = 0;
						my $delete = 0;
						foreach $t(keys %candidate){
							my @w2 = sort{$a<=>$b}(@{$bound{$u}{$t}},@{$bound{$u}{$s}});
							my $len1 = $bound{$u}{$t}->[1] - $bound{$u}{$t}->[0] + 1;
							my $len2 = $bound{$u}{$s}->[1] - $bound{$u}{$s}->[0] + 1;
							my $overlap = $len1 + $len2 - ($w2[-1]-$w2[0]+1);
							$skip = 1 if $overlap > $opt->overlap;
							$delete = 1 if $overlap / $len2 >= 0.9;
						}
						if($skip == 1){
							$remove{$s} = 1 if $delete == 1;
							next;
						}
						$candidate{$s} = 1;
						while(@w=splice(@tt,0,4)){
							push @range,($s,@w,$u,$strand->{$s}{$u});
						}
					}
				}
			$clean{$u} = \@range;
	}
	open OUT,">$output";
	print OUT "Scaf\tScafSt\tScafEnd\tUnitigSt\tUnitigEnd\tUnitig\tStrand\n";
	foreach $i(keys %clean){
		my @k = @{$clean{$i}};
		while(@tt=splice(@k,0,7)){print OUT join("\t",@tt),"\n"}
	}
	close OUT;
	open DEL,">$delFile";
	foreach(keys %remove){print DEL $_,"\n"};
	close DEL;

	return(\%clean,\%remove,\%bound);
}

sub read_coords{
	my $file = shift;
	my $minID = shift;
	my (%data, %strand,%match);
	open COORDS,$file || die;
	while(<COORDS>){
		next unless /^\d+/;
		my @s = split;
		next if $s[6] < $minID;
		push @{$data{$s[13]}{$s[14]}},$_;
		$match{$s[13]}{$s[14]} //= 0;
		if($s[5]>$match{$s[13]}{$s[14]}){
			$strand{$s[13]}{$s[14]} = $s[12];
			$match{$s[13]}{$s[14]} = $s[5];
		}
	}
	close COORDS;
	return(\%data,\%strand);
}

sub get_range{
	my $d = shift;
	my (%ref,%qry,@tmp,@final,$qsize,$rsize); 
	while(@t=splice(@$d,0,4)){
		$ref{$t[0]}{$t[1]} = [$t[2],$t[3]];
	}
	foreach $i(sort{$a<=>$b} keys %ref){
		my @j = sort{$b<=>$a} keys %{$ref{$i}};
		if(!@tmp || $tmp[-3]<$i){push @tmp,($i,$j[0],@{$ref{$i}{$j[0]}})}
		else{
			$tmp[-3] = $tmp[-3] < $j[0] ? $j[0] : $tmp[-3];
			my @order = sort{$a<=>$b} ($tmp[-2],$tmp[-1],@{$ref{$i}{$j[0]}});
			($tmp[-2],$tmp[-1]) = ($order[0],$order[-1]);
		}
	}
	while(@k=splice(@tmp,0,4)){
		$qry{$k[2]}{$k[3]} = [$k[0],$k[1]];
	}
	foreach $i(sort{$a<=>$b} keys %qry){
	    my @j = sort{$b<=>$a} keys %{$qry{$i}};
			if(!@final || $final[-1]<$i){push @final,(@{$qry{$i}{$j[0]}}),$i,$j[0]}
			else{
				$final[-1] = $final[-1] < $j[0] ? $j[0] : $final[-1];
				@order = sort{$a<=>$b} ($final[-4],$final[-3],@{$qry{$i}{$j[0]}});
				($final[-4],$final[-3]) = ($order[0],$order[-1]);
			}
	}
	my @tt = @final;
	my $left = 100000000000;
	my $right = 0;
	while(@pp = splice(@tt,0,4)){
		$rsize += $pp[1]-$pp[0]+1;
		$qsize += $pp[3]-$pp[2]+1;
		$left = $left < $pp[2] ? $left : $pp[2];
		$right = $right > $pp[3] ? $right : $pp[3];
	}
	return(\@final,$rsize,$qsize,$left,$right);
}

sub read_seq{
	my $file = shift ;
	my %data = ();
	my $io=new Bio::SeqIO(-file=>$file,-format=>'fasta');
	while($seq=$io->next_seq){$data{$seq->id}=$seq->seq}
	$io->close;
	return (\%data);
}

sub locate_gaps{
	my %data = ();
	my ($num,$size) = (0,0); 
	my $io=new Bio::SeqIO(-file=>$_[0],-format=>'fasta');
	while($seq=$io->next_seq){
		my $xl = $seq->seq;
		while($xl=~/[Nn]+/g){
			my $st = length($`)+1;
			my $end = length($`) + length($&);
			$data{$seq->id}{$st} = $end;
			$num++;
			$size += length($&);
		}
	}
	$io->close;
	LOG("Gap number: ".$num."; Gap size: ".$size);
	return (\%data);
}

sub LOG {
	$| = 0;
	my $print = shift;
	my $time = localtime;
	print STDERR color("green"),'[',$time,'] ',color("red"),$print,"\n",color("reset");
}

sub ERROR {
	my $print = shift;
	print STDERR "\n====================================================\n";
	print STDERR "Error:  ",color("green"),$print,color("reset"),"\n";
	print STDERR "====================================================\n";
}


