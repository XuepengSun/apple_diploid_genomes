#!/usr/bin/perl -w
use Bio::SeqIO;
use List::Util qw(max min);

my ($gentic_map,$synteny_map) = @ARGV;

$marker_thres = 5;

my (%hash,%dist,%pos,%comp);

open GMAP,$gentic_map || die;
while(<GMAP>){
	next if $.==1;
	chomp;
	my @s = split /,/;
	$hash{$s[0]}{$s[2]}++;
	push @{$pos{$s[0]}{$s[2]}},$s[1];
}
close GMAP;

open SY,$synteny_map || die;
while(<SY>){
	chomp;
	my @s = split;
	my ($chr)=$s[-1]=~/(.*):/;
	my $max = max($s[1],$s[2]);
	$comp{$s[0]}{$chr}++;
	push @{$dist{$s[0]}{$chr}},$max;
}
close SY;

foreach $scaf(keys %hash){
	my %tmp = ();
	my @keys = sort{$hash{$scaf}{$b}<=>$hash{$scaf}{$a}} keys %{$hash{$scaf}};
	if(scalar @keys > 1){
		if($hash{$scaf}{$keys[1]} >= $marker_thres){
			print STDERR $scaf,"\t",$keys[0],"|",$hash{$scaf}{$keys[0]},"\t",$keys[1],"|",$hash{$scaf}{$keys[1]},"\n";
			foreach(@{$pos{$scaf}{$keys[0]}}){$tmp{$_} = 1}
			foreach(@{$pos{$scaf}{$keys[1]}}){$tmp{$_} = 2}
			
			my @t = sort{$a<=>$b} keys %tmp;
			my $class = "";
			my ($st,$end,%pp);

			while(@t){
				my $v = shift @t;
				if($class eq ""){
					($class,$st,$end) = ($tmp{$v},$v,$v);
				}
				else{
					if($tmp{$v} ne $class){
						my $j=$st.'-'.$end;
						push @{$pp{$class}},$j;
						($class,$st,$end) = ($tmp{$v},$v,$v);
					}
					else{$end = $v}
				}
			}

			$j=$st.'-'.$end;
			push @{$pp{$class}},$j;
			
			foreach $i(keys %pp){
				print "genetic\t$scaf";
				foreach(@{$pp{$i}}){print "\t",$_}
				print "\n";
			}

			%tmp=();
			my @k = sort{$comp{$scaf}{$b}<=>$comp{$scaf}{$a}} keys %{$comp{$scaf}};	
			foreach(@{$dist{$scaf}{$k[0]}}){$tmp{$_} = 1}
			foreach(@{$dist{$scaf}{$k[1]}}){$tmp{$_} = 2}
			
			@t = sort{$a<=>$b} keys %tmp;
                        $class = "";
			%pp=();
			$st = "";
			$end = "";
                        while(@t){
                                my $v = shift @t;
                                if($class eq ""){
                                        ($class,$st,$end) = ($tmp{$v},$v,$v);
                                }
                                else{
                                        if($tmp{$v} ne $class){
                                                my $j=$st.'-'.$end;
                                                push @{$pp{$class}},$j;
                                                ($class,$st,$end) = ($tmp{$v},$v,$v);
                                        }
                                        else{$end = $v}
                                }
                        }

                        $j=$st.'-'.$end;
                        push @{$pp{$class}},$j;

                        foreach $i(keys %pp){
                                print "syntenic\t$scaf";
                                foreach(@{$pp{$i}}){print "\t",$_}
                                print "\n";
                        }

		}
	}
}


