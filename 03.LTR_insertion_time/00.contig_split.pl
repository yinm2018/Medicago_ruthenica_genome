use strict;
use warnings;
use Bio::SeqIO;

my $file=shift;
my $fa=Bio::SeqIO->new(-file=>$file,-format=>'fasta');
my(%gene);
while(my $seq_obj=$fa->next_seq){
    my $id=$seq_obj->id;
    my $seq=$seq_obj->seq;
    open (O,">contig_split\/$id.fa");
    print O ">$id\n$seq\n";
    close O;
}
