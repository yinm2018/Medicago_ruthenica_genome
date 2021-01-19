use strict;
use warnings;
use Bio::SeqIO;

my @in=<align_soloLTRs/*/*3.5.ltr.fa.align.dnadist>;
open (O,">$0.out");
print O "ID\tdnadist\n";
for my $in (@in){
    $in=~/align_soloLTRs\/(\S+)\//;
    my $name=$1;
    my $fasta="align_soloLTRs/$name/$name.3.5.ltr.fa.align";
    open (F,"$in");
    <F>;
    my $line=<F>;
    chomp $line;
    close F;
    my @line=split(/\s+/,$line);
    my $solofa=Bio::SeqIO->new(-format=>"fasta",-file=>"$fasta");
    while (my $seq=$solofa->next_seq) {
        my $id=$seq->id;
        $id=~/(.*)\.def/ or die "$id\n";
    }
    print O "$name\t$line[2]\n";
}
close O;


