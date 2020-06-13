#!/usr/bin/perl 
use warnings;
use Term::ANSIColor;
use Bio::SeqIO;
use Getopt::Long::Descriptive;
use Bio::SeqIO;
use threads;
use File::Find;

# version 03-11-2020


my ($opt, $usage) = describe_options(
  'Usage: %c %o <args>',
  [],
  ['dscaf|k=s',"donor or query scaffold sequence in fasta format", {required => 1}],
  ['dchr|d=s',"donor or query chromosome sequence in fasta format", {required => 1}],
  ['dord|n=s',"donor scaffold order on chromosome, [chr\tscaffold\tstrand]", {required => 1}],
  ['rscaf|p=s',"receptor or reference scaffold sequence in fasta format", {required => 1}],
  ['rchr|r=s',"receptor or reference chromosome sequence in fasta format", {required => 1}],
  ['rord|j=s',"receptor scaffold order on chromosome, [chr\tscaffold\tstrand]", {required => 1}],
	['coords|c=s',"coords file from mummer", {required => 1}],
  ['minalnlen|l=i',"minimal alignment length to use, default: 10000",{default  => 10000}],
	['maxinterval|m=i',"maximal interval between anchors to allow, default: 1000000",{default=>1000000}],
  ['tmpdir|f=s',"temp fold, default: temp_gapCloser",{default  => "temp_gapCloser"}],
  ['help|h', "print usage message and exit", {shortcircuit => 1}],
	['version|v',"version 0.2",{shortcircuit => 1}],
  []
);
print($usage->text), exit if $opt->help;
print STDERR "\n";
LOG("checking gaps in reference sequences");

unless (-d $opt->tmpdir){mkdir $opt->tmpdir}
my $info = $opt->tmpdir.'/seq_substitution_info.txt';
my $raw = $opt->tmpdir.'/seq_substitution_raw.txt';
my $chrFile = $opt->tmpdir.'/new_chromosome.fa';
my $scfFile = $opt->tmpdir.'/new_scaffold.fa';
my $summary = $opt->tmpdir.'/summary.txt';
my $fagp = $opt->tmpdir.'/new_agp.txt';

open INFO,">$info";
open CHR,">$chrFile";
open SCF,">$scfFile";
open SUM,">$summary";
open FAGP,">$fagp";
open RAW,">$raw";

print INFO "chr\trefSTART\trefEND\tqrySTART\tqryEND\n";
print RAW "chr\trefSTART\trefEND\tqrySTART\tqryEND\ttype\tGaps\n";


my $rGaps = &find_gap($opt->rchr);
my $dscafLen = &read_seq($opt->dscaf);
my $dchrSeq = &read_seq($opt->dchr);
my $rchrSeq = &read_seq($opt->rchr);
my $dscafOrder = &find_order($opt->dord,$dscafLen);
my $fake = &fake_chr($opt->rord,$opt->rscaf);

LOG("checking reference scaffold orders");
check_fake_chr($fake,$rchrSeq);

LOG("processing gaps");
my %coords = ();
my %crdsToScaf = ();
open COORDS,$opt->coords || die;
while(<COORDS>){
	my $line = $_;
	$line=~s/^\s+|\s+$//g;
	next unless $line=~/^\d+/;
	my @s = split(/\s+/,$line);	
	my ($qchr) = $s[14]=~/(chr\d+)/i;
	my ($st,$end) = $s[2]<$s[3]?($s[2],$s[3]):($s[3],$s[2]);
	push @{$coords{$s[13]}},$line;
	foreach $i(keys %{$dscafOrder->{$qchr}}){
		if($i<=$st && $dscafOrder->{$qchr}{$i}->[0]>=$end){
			push @{$crdsToScaf{$s[13]}{$dscafOrder->{$qchr}{$i}->[1]}},$line;
			#print STDERR $dscafOrder->{$qchr}{$i}->[1],"\t",$line,"\n";
		}
	}
}
close COORDS;


my %new = ();
my $add = 0;
my $remove = 0;
foreach my $c(keys %crdsToScaf){
	my %tmp = map{$_=>length($dscafLen->{$_})} keys %{$crdsToScaf{$c}};
	my @list = sort{$tmp{$b}<=>$tmp{$a}} keys %tmp;
	my %ref = ();
	foreach my $i(@list){
		my @tmpqry = ();
		my @tmpref = ();
		my @range = sort{$a<=>$b} keys %ref;
		foreach(@{$crdsToScaf{$c}{$i}}){
			my @t = split;
			next if $t[4] < $opt->minalnlen && $t[5] < $opt->minalnlen;
			my ($st,$end) = $t[2]<$t[3] ? ($t[2],$t[3]):($t[3],$t[2]);
			next unless (!@range || ($t[0]>=$range[-1] || $t[1]<=$range[0]));

			## or change below to $tmpref[-1]<$t[0] and $tmpqry[-1]<$st to minimize substitution
			if((!@tmpref || !$tmpref[-1]<$t[0]-1) && (!@tmpqry || $tmpqry[-1]<$st-1)){ 
				push @tmpqry,($st,$end);
				push @tmpref,($t[0],$t[1]);
			}
			else{
				$tmpref[-1] = $tmpref[-1] < $t[1] ? $t[1] : $tmpref[-1];
				$tmpqry[-1] = $tmpqry[-1] < $end ? $end : $tmpqry[-1];
			}
		}
		while(@ss = splice(@tmpref,0,2)){
			my @qq = splice(@tmpqry,0,2);
			$ref{$ss[0]}=[$ss[1],$i,$qq[0],$qq[1]];
		}
	}
	my @replace = ();
	my @sites = sort{$a<=>$b} keys %ref;
	while(@p = splice(@sites,0,1)){
		my ($end,$scf,$qst,$qend) = @{$ref{$p[0]}};
		my $opGap = 0;
		foreach(keys %{$rGaps->{$c}}){
			next if $p[0] > $rGaps->{$c}{$_} || $end < $_;
			$opGap++;
		}
		if($opGap > 0){
			push @replace,($p[0],$end,$qst,$qend);
			print RAW join("\t",($c,$p[0],$end,$qst,$qend,"within_anchor",$opGap)),"\n";
		}
		if(@sites){
			my $next = $sites[0];
			my @nn = @{$ref{$next}};
			if($scf eq $nn[1]){
				my $left = $end + 1;
				my $right = $next - 1;
				$opGap = 0;
				foreach(keys %{$rGaps->{$c}}){
					next if $left > $rGaps->{$c}{$_} || $right < $_;
					$opGap++;
				}
				if($opGap > 0){
					my $qqleft = $qend + 1;
					my $qqright = $ref{$next}->[2]-1;
					if(($qqright-$qqleft+1)-($right-$left+1)<=$opt->maxinterval){    ## v0.2
						push @replace,($left,$right,$qqleft,$qqright);
						print RAW join("\t",($c,$left,$right,$qqleft,$qqright,"between_anchor",$opGap)),"\n";
					}
				}
			}
		}
	}
	$offset = 0;
	while(@pp = splice(@replace,0,4)){
		print INFO $c,"\t",join("\t",@pp),"\n";
		if($pp[0] -1 > $offset){
			$pre = substr($fake->{$c},$offset,$pp[0]-1-$offset);
		}
		else{$pre = ""}
		$subseq = substr($dchrSeq->{$c},$pp[2]-1,$pp[3]-$pp[2]+1);
		$new{$c}.=$pre.$subseq;
		$offset = $pp[1];
		$remove += $pp[1]-$pp[0]+1;
		$add += $pp[3]-$pp[2]+1
	}
	if($offset < length($rchrSeq->{$c})){
		$new{$c}.= substr($fake->{$c},$offset);
	}
	## add long ends of query chromosomes
	my @head = split(/\s+/,$coords{$c}[0]);
	my @tail = split(/\s+/,$coords{$c}[-1]);
	if($head[4]>=300 && $head[2] > $head[0] * 2 && $head[2] > 1000){
		$new{$c} = substr($dchrSeq->{$c},0,$head[2]-1).substr($new{$c},$head[0]-1);
		$add += $head[2] - 1;
	}
	if($tail[4]>=300 && $tail[8]-$tail[3]>$tail[7]-$tail[1]*2 && $tail[8]-$tail[3]>1000){
		my $plen = length($new{$c}) - ($tail[7]-$tail[1]);
		$new{$c} = substr($new{$c},0,$plen).substr($dchrSeq->{$c},$tail[3]);
		$add += $tail[8] - $tail[3];
	}
}

my $scfNum = 0;
foreach my $j(keys %new){
	print CHR '>',$j,"\n",$new{$j},"\n";
	my @split = split(/\++/,$new{$j});
	foreach(@split){
		$scfNum ++;
		print SCF '>scaffold',$scfNum,'n',"\n",$_,"\n";
		print FAGP $j,"\t",'scaffold',$scfNum,'n',"\t+\n";
	}
}

print SUM "Nucleotides deleted from reference: ",$remove,"\n";
print SUM "Nucleotides added to the reference:",$add,"\n";

close CHR;
close SUM;
close SCF;
close FAGP;
close INFO;
close RAW;

sub check_fake_chr{
	my ($fake,$origin) = @_;
	my $check = 0;
	foreach $r(keys %$fake){
			while($fake->{$r}=~/\++/g){
				my $pre = length($`);
				my $match = length($&) + length($`);
				my $preseq1 = substr($fake->{$r},$pre-50,50);
				my $flwseq1 = substr($fake->{$r},$match,50);
				my $preseq2 = substr($origin->{$r},$pre-50,50);
				my $flwseq2 = substr($origin->{$r},$match,50);
				$check++ if $preseq1 ne $preseq1 || $flwseq1 ne $flwseq2;
			}
	}
	if($check == 0){LOG("chromosome check done")}
	else{ERROR("Error in reference AGP file") & exit}
}

sub fake_chr{
  my ($agp,$scf) = @_;
  my %rscf = ();
  my %fakechr =();
  my %tmpAGP = ();
  my $io=new Bio::SeqIO(-file=>$scf,-format=>'fasta');
  while(my $seq=$io->next_seq){$rscf{$seq->id}=$seq->seq}
  $io->close;

  open AGP,$agp || die;
  while(<AGP>){chomp;my @s =split;push @{$tmpAGP{$s[0]}},($s[1],$s[2])}
  close AGP;

  foreach my $i(keys %tmpAGP){
    while(my @t = splice(@{$tmpAGP{$i}},0,2)){
      my $xl = $rscf{$t[0]};
      if($t[1] eq '-'){
        $xl=join("",reverse split("",$xl));
        $xl=~tr/ATGCatgc/TACGtacg/;
      }
      $fakechr{$i}.=$xl;
      if(@{$tmpAGP{$i}}){$fakechr{$i} .= '+' x 100}
    }
  }
  return(\%fakechr);
}

sub find_order{
	my ($agp,$len) = @_;
	my %dList = ();
	my %done = ();
	my $offset = 0;
	open IN,$agp || die;
	while(<IN>){
		my @s = split;
		if(!$done{$s[0]}){
			$done{$s[0]} = 0;
			$offset = 0;
		}
		my $st = $offset + $done{$s[0]} * 100 + 1;
		my $end = length($len->{$s[1]}) + $st - 1;
		$dList{$s[0]}{$st} = [$end,$s[1]];
		$done{$s[0]} ++ ;
		$offset += length($len->{$s[1]});
		#print $s[0],"\t",$s[1],"\t",length($len->{$s[1]}),"\t",$st,"\t",$end,"\n";
	}
	close IN;
	return(\%dList);
}

sub read_seq{
	my $file = shift;
	my %tmp = ();
	my $io=new Bio::SeqIO(-file=>$file,-format=>'fasta');
	while($seq=$io->next_seq){
		$tmp{$seq->id} = $seq->seq;
	}
	$io->close;
	return(\%tmp);
}

sub find_gap{
	my $file = shift;
	my %tmp = ();
	my ($st,$end) = ();
	my $io=new Bio::SeqIO(-file=>$file,-format=>'fasta');
	while($seq=$io->next_seq){
		my $xl = $seq->seq;
		while($xl=~/[Nn]+/g){
			$st = length($`) + 1;
			$end = length($`) + length($&);
			$tmp{$seq->id}{$st} = $end;
		}
	}
	$io->close;
	return(\%tmp);
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

