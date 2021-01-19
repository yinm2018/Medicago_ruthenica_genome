#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $file = shift;


my $alength = 0;
my $allen = "$file.total.len";
open S, ">$file.sca.len";   
open L, ">$allen";          

my $fa = Bio::SeqIO->new(-file=>$file, -format=>'fasta');
while(my $seq_obj = $fa->next_seq){
    my $id = $seq_obj->id;
    my $seq = $seq_obj->seq;   
    my $len = length($seq);
    print S "$id\t$len\n";
    $alength += $len;
}
print L "$alength\n";

close L;
close S;

