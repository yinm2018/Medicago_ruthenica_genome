#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use Getopt::Long;
use MCE::Map;
use v5.10;
use zzIO;

#my $GCstats_file='./GCstat_merge.txt';

my $window = 1000;
my($in_vcf, $in_indel, $opt_help);
my $type='snp';
my $thread=10;
GetOptions (
	'help|h!' => \$opt_help,
	'vcf=s' => \$in_vcf,
	'indel=s' => \$in_indel,
	'p=i' => \$thread,
	't=s' => \$type,
);

die "usage: perl xxx.pl -vcf raw_snp.vcf.gz -indel raw_indel.vcf.gz [-p 10] (| bgzip -c ) > out.vcf\n" if ($opt_help || !-e $in_vcf || !-e $in_indel );

print STDERR "\n** Now in node: ";
print STDERR`hostname`;

foreach ($in_indel, $in_vcf){
	die "file not exist: $_\n" unless -e $_;
}

print STDERR "* Processing indel\n";
my %indel = &read_indel($in_indel);

my $INVCF = open_in_fh($in_vcf);
my @idline;
print STDERR "Processing snp_vcf: $in_vcf\n";
while ( <$INVCF> ){
	if ( /^##/ ){
		print $_;
		next;
	}
	if (/^#/){
		chomp;
		print $_."\n";
		print STDERR $_."\n";
		if (/^#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+\w+/){
			@idline=split(/\s+/,$_);
			last;
		}else{
			print STDERR;
			last;
		}
	}
}
die "NO idline: \n" if ! @idline;

MCE::Map::init {
   chunk_size => 'auto',
   max_workers => $thread,
};

if ($thread == 1) {
	while (<$INVCF>) {
		my $result = &flt($_);
		say $result if defined $result;
	}
} else {
	my @t = grep{defined$_} mce_map_f { &flt } $INVCF;
	say STDERR 'run ok';
	MCE::Map::finish;
	say $_ foreach @t;
}

close $INVCF;
print STDERR "/home/share/software/samtools/htslib/tabix -p vcf xxxOUT.vcf.gz\n";
# system("/home/share/software/samtools/htslib/tabix -p vcf $outvcf");

THEEND:

print STDERR "\n** ALL Done\n";
print STDERR`date`;

exit 0;
sub readdp{
	my ($in)=@_;
	my %r;
	my $D;
	open ($D,"$in")||die"no $in\n";
	while (<$D>) {
		chomp;
		next if /^#/;
		my @a=split(/\s+/,$_);
		$r{$a[0]}=$a[5];
		# $r{$a[0]}=$a[6];
	}
	close $D;
	return %r;
}

sub read_rep {
	my $F;
	my ($repfile)=@_;
	my %rep;
	# print STDERR `wc -l $repfile`;
	open ($F,"$repfile") or die "can't open rep_gff3_file, $repfile, $!"; #gff3
	while ( <$F> ) {
		chomp;
		## gff3:
		##   0     1      2    3    4    5      6     7       8
		## seqid source type start end score strand phase attributes
		next if /^#/;
		next if /^$/;
		my @a=split(/\s+/,$_);
		next if $a[1] eq "TRF";
		#$rep{$a[0]}{$a[3]}=$a[4];
		my $start_1w=int($a[3]/$window);
		my $end_1w=int($a[4]/$window);
		my $sub=$end_1w-$start_1w;
		if ($sub==0) {
			$rep{$a[0]}{$start_1w}{$a[3]}=$a[4];
		}elsif ($sub==1){
			$rep{$a[0]}{$start_1w}{$a[3]}=$end_1w * $window;
			$rep{$a[0]}{$end_1w}{$start_1w * $window} = $a[4];
		}else{
			$rep{$a[0]}{$start_1w}{$a[3]}=($start_1w+1) * $window;
			$rep{$a[0]}{$end_1w}{($end_1w-1)*$window}=$a[4];
			foreach (1..$sub-1){
				$rep{$a[0]}{$start_1w+$_}{ALL}++;
			}
		}
	}
	close $F;
	return %rep;
}



sub read_indel {
	my ($inindel)=@_;
	my %indel;
	my $F;
	 ($inindel=~/\.gz$/) ? open ($F,"zcat $inindel|") : open ($F,"< $inindel") ;
	while ( <$F> ) {
		##  VCF:
		##   0    1   2  3    4   5      6     7    8      9...
		## CHROM POS ID REF  ALT QUAL FILTER INFO FORMAT   ...
		chomp;
		next if $_=~/^#/;
		my @a=split(/\s+/,$_);
		my $info=get_info($a[7]);
		## 舍去 QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0
		next if ( exists($$info{'QD'}) && $$info{'QD'} < 2 );
		next if ( exists($$info{'FS'}) && $$info{'FS'} > 200 );
		next if ( exists($$info{'ReadPosRankSum'}) && $$info{'ReadPosRankSum'} < -20 );
		next if ( exists($$info{'SOR'}) && $$info{'SOR'} > 10.5 );
		next if ( exists($$info{'InbreedingCoeff'}) && $$info{'InbreedingCoeff'} < -0.8 );

		# next if $a[6] ne "PASS";
		$indel{$a[0]}{$a[1]+$_}=0 foreach (-10..10);  #### indel +- 10bp 筛选
	}
	close $F;
	return %indel;
}



sub flt{
	my $line=$_;
	chomp $line;
	return $line if $line=~/^#/;
	my @a=split(/\s+/,$line);

	##  VCF:
	##   0    1   2  3    4   5      6     7    8      9...
	## CHROM POS ID REF  ALT QUAL FILTER INFO FORMAT   ...
	my $info=get_info($a[7]);

	return undef if ( %indel and exists $indel{$a[0]}{$a[1]} );

	$a[6]='PASS';
	return join("\t",@a);
}

sub get_info() {
	my $in = $_[0] or die;
	my %ret;
	my @a = split ';',$in;
	foreach my $a (@a) {
		my @b = split '=',$a;
		$ret{$b[0]} = defined $b[1] ? $b[1] : 0;
	}
	return \%ret;
}





