# r = 1.3e-8
# t = K / 2r
use strict;
use warnings;

my $file = shift;
my $species = shift;

my $r = 1.3E-8;    
my $t;
my $T;
open F,$file;

open O,">$species\.insert_time\.tab";
print O"ID\tdnadist\tinsert_time\n";

open T,">$species\.Mya\.tab";
print T"Mya\n";

while (<F>){
    chomp;
    next if /ID/;
    my @a = split /\t/;
    $t = $a[1] / (2 * $r);
    $T = $t / 1000000;
    print O"$[0]\t$a[1]\t$r\t$t\t$T\n";
    print T"$T\n";
}
close F;
close O;
close T;
