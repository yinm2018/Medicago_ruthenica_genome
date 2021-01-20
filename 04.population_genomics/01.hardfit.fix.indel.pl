use strict;
use warnings;

$a = shift;
$b = shift;

open IN, "zcat $a |";
open OUT, "| gzip - > $b";
while(<IN>){
chomp;
if(/^#/){
	print OUT "$_\n";
 	next;
}
@c = split/\s+/;
if($c[6] =~ /my_snp_filter/){
	next;
}else{
	print OUT "$_\n";
	}	
}




