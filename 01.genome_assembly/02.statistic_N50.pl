#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $fasta;
my $file = shift;

$file =~ /\/.*\/(.*)/g;
my $filename =$1;
my %len;

open (O,">N50_of_Scaffolds.txt") || die ($!);
$fasta = Bio::SeqIO ->new (-file =>$file,-format=>'fasta');
while (my $seq_obj = $fasta -> next_seq ){
      my $name = $seq_obj -> id;
      my $seq = $seq_obj -> seq;
      my $length = $seq_obj -> length;
      $len{$name} = $length;            
  }

my @sort_hash = sort {$len{$a}<=>$len{$b}}(keys %len); 
my $seq_num = @sort_hash;                              
my $whole_length = 0;

foreach my $hash_key(@sort_hash){
    chomp $hash_key;
    $whole_length += $len{$hash_key};
}

my $half_len =  $whole_length/2;
my $number = 0;
foreach my $hash(@sort_hash){
    $number += $len{$hash};
    if($number >= $half_len){
        print O " the N50 of contigs/scaffolds is $len{$hash} ";
        exit;
    }
}
close O;
