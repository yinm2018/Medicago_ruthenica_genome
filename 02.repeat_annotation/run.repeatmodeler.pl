#!/usr/bin/perl

use strict;
use warnings;

my $ref=shift or die "perl $0 ref.fa\n";
die "no such file: $ref\n" if (! -e "$ref");
my $repeatmask="your_RepeatMasker_sorftware_path";
my @in;
@in=<RM*/consensi.fa.classified>;
die"no classified files\n"  if scalar(@in) != 1;
print "$repeatmask -lib $in[0] -pa 30 $ref 2>&1 | tee 03.RepeatMasker.log\n";
#system ("$repeatmask -lib $in[0] -pa 30 $ref 2>&1 | tee 03.RepeatMasker.log");
