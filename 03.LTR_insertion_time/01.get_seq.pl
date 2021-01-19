#! /usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

#system("rm -rf ltrfinder.five_prime_LTR.fa ltrfinder.three_prime_LTR.fa");

my $genome_fasta = shift; #基因组序列
my @seqfile = <./ltr_finder_out/*.out>;
my %h;

for my $seq_file (@seqfile){
    chomp $seq_file;
    $seq_file =~ /.*\/(.*)\.fa\.out/;
    my $seqid = $1;
    open (SEQ, "< $seq_file") || die "$!";
    while (<SEQ>){
	chomp;
	next unless (/^\[/);
	my @a = split /\t/;
        my $index = $a[0];
        $index =~ /\[\s*(\d+)\]/;
	my $name = "def_".$1;
	my ($start1, $end) = (split /-/, $a[2])[0,1];
	my ($len1, $len2) = (split /,/, $a[3])[0,1];
	my $strand = $a[-4];
	my $start2 = $end - $len2 + 1;
	my $len0 = $end - $start1 + 1;
	$h{$seqid}{$name}{start1} = $start1;
	$h{$seqid}{$name}{start2} = $start2;
	$h{$seqid}{$name}{strand} = $strand;
	$h{$seqid}{$name}{len1} = $len1;
	$h{$seqid}{$name}{len2} = $len2;
	$h{$seqid}{$name}{len0} = $len0;
    }
}

open (O,">ltrfinder.seq.fa");
open (O1,">ltrfinder.five_prime_LTR.fa");
open (O2,">ltrfinder.three_prime_LTR.fa");

my $fa = Bio::SeqIO ->new(-file=>"$genome_fasta",-format=>"fasta");
while (my $seq = $fa ->next_seq){
    my $id = $seq ->id;
    my $seq = $seq ->seq;
    if (exists $h{$id}){
	for my $key (sort keys %{$h{$id}}){
	    my ($start1, $len1, $strand) = ($h{$id}{$key}{start1}, $h{$id}{$key}{len1}, $h{$id}{$key}{strand});
	    my ($start2, $len2) = ($h{$id}{$key}{start2}, $h{$id}{$key}{len2});
	    my $len0 = $h{$id}{$key}{len0};
	    my $newseq0 = substr($seq, $start1, $len0);
	    my $newseq1 = substr($seq, $start1, $len1);
	    my $newseq2 = substr($seq, $start2, $len2);
	    if ($strand eq "-"){
		$newseq1 = reverse($newseq1);
		$newseq1 =~ tr/[ATCG]/[TAGC]/;
		$newseq2 = reverse($newseq2);
		$newseq2 =~ tr/[ATCG]/[TAGC]/;
		$newseq0 = reverse($newseq0);
		$newseq0 =~ tr/[ATCG]/[TAGC]/;
	    }
	    print O1 ">$id.$key\n$newseq1\n";
	    print O2 ">$id.$key\n$newseq2\n";
	    print O ">$id.$key\n$newseq0\n";
	}
    }
}
    
close O1;
close O2;


