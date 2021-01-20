#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $file = shift;

my $sum_GC_num=0;
my $sum_len=0;
my $sum_GC_content;
my %count;

open (R,">Whole_genome_GC_content.txt") || die ($!);
open (O,">per_Scaffold_GC_content.txt") || die($!);

my $fa = Bio::SeqIO->new(-file=>"$file",-format=>'fasta');
while(my $seq = $fa->next_seq){
    my $id = $seq->id;
    my $seq = $seq->seq;                    

    my @a = split//,$seq;
    my $num=0;
    foreach my $code(@a){
        if ($code =~ m/G|C/){
            $num++;                         
        }else{
            next;
        }
    }
    my $len = length $seq;                  
    my $GC_per = $num/$len;                 
    $sum_len = $len + $sum_len;             
    $sum_GC_num += $num;                   
    print O "$id\t$len\t$num\t$GC_per\n";   
}
$sum_GC_content = $sum_GC_num/$sum_len;     
print R "Whole genome lenth: $sum_len\nWhole genome GC number: $sum_GC_num\nWhole genome GC content: $sum_GC_content\n";
close O;
close R;
