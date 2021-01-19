use strict;
use warnings;
use Bio::SeqIO;

my $mucle="~user/software/muscle3.8.31/muscle3.8.31_i86linux64"; # path-to-mucle;
my $clustalw2="~user/software/ClustalW/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2"; # path-to-clustalw2;
my $dnadist="~user/software/phylip/phylip-3.696/exe/dnadist"; # path-to-dnadist;
die "please check the path of the software used in this script\n" if ((! -e "$mucle") || (! -e "$clustalw2") || (! -e "$dnadist"));

my $dir="align_soloLTRs";
`mkdir $dir` if (! -e "$dir");

my %seq;
&readfasta("ltrfinder.three_prime_LTR.fa");
&readfasta("ltrfinder.five_prime_LTR.fa");

my $id_name;
open IN, "< ltrfinder.three_prime_LTR.fa";
# system("rm $0.sh");
open (SH, "> $0.sh");
while (<IN>){
    chomp;
    if (/>(.*)/){
	$id_name = $1;
	next;
    }
    my $dir2="$dir/$id_name";
    # print STDERR "$dir2";
    `mkdir -p $dir2` if (! -e "$dir2");
    open (O,">$dir2/$id_name.3.5.ltr.fa");
    print O ">$id_name#five_prime_LTR\n$seq{$id_name}{five_prime_LTR}\n";
    print O ">$id_name#three_prime_LTR\n$seq{$id_name}{three_prime_LTR}\n";
    close O;
    open (O,">$dir2/$id_name.3.5.ltr.fa.align.phylip");
    print O "$id_name.3.5.ltr.fa.align.phy\nD\nY\n";
    close O;
    # open SH, ">>$0.sh";
    print SH "cd $dir2; $mucle -in $id_name.3.5.ltr.fa -out $id_name.3.5.ltr.fa.align; $clustalw2 -INFILE=$id_name.3.5.ltr.fa.align  -CONVERT -TYPE=DNA -OUTFILE=$id_name.3.5.ltr.fa.align.phy -OUTPUT=PHYLIP; $dnadist < $id_name.3.5.ltr.fa.align.phylip; mv outfile $id_name.3.5.ltr.fa.align.dnadist; cd ../../\n";
    # close SH;
}
# close SH;

sub readfasta{
    my ($infile_fa)=@_;
    die "no $infile_fa\n" if (! -e "$infile_fa");
    $infile_fa=~/(five_prime_LTR|three_prime_LTR)/;
    my $type=$1;
    my $infile_fa_read=Bio::SeqIO->new(-file=>"$infile_fa",-format=>"fasta");
    while (my $seq=$infile_fa_read->next_seq){
	my $id=$seq->id;
	my $seq=$seq->seq;
	$seq{$id}{$type}=$seq;
    }
}
