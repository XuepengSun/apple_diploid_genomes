#!/usr/bin/perl -w
use Bio::SeqIO;
use Bio::Seq;

unless(scalar @ARGV == 2){
	print "\nUseage: perl $0 <genome> <gtf file>\n\n";
	exit;
}


my ($genome,$gtf)=@ARGV;

my (%sequence,%strand,%gene);

my $io=new Bio::SeqIO(-file=>$genome,-format=>'fasta');
while(my $seq=$io->next_seq){
	$sequence{$seq->id}=$seq->seq;
}

open GTF,$gtf || die "cannot find gtf file";
while(<GTF>){
	chomp;
	my @s=split /\t+/;
	if($s[2] eq 'CDS'){
		my ($id)=$s[8]=~/gene_id\s+\"(.*?)\"/;
		$gene{$id}.=substr($sequence{$s[0]},$s[3]-1,$s[4]-$s[3]+1);
		$strand{$id}=$s[6];
	}
}
close GTF;

open NUCL,'>output.fna';
open PROT,'>output.faa';

foreach $gn(keys %gene){
	my $seqobj=Bio::Seq -> new(-seq=>$gene{$gn},-alphabet=>'dna');
	if($strand{$gn} eq '+'){
		print NUCL '>',$gn,"\n",$gene{$gn},"\n";
		my $trans=$seqobj->translate;
		print PROT '>',$gn,"\n",$trans->seq,"\n";
	}
	else{
		my $nucl=$seqobj->revcom;
		print print NUCL '>',$gn,"\n",$nucl->seq,"\n";
		my $trans=$seqobj->revcom->translate;
                print PROT '>',$gn,"\n",$trans->seq,"\n";
	}
}
close NUCL;
close PROT;

