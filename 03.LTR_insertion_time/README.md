# Genome size evolution of Medicago ruthenica

## LTR insertion time

[00.contig_split.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/03.LTR_insertion_time/00.contig_split.pl)

```
mkdir contig_split ltr_finder_out

perl 00.contig_split.pl Genome.fa
```

The prediction of full-length LTR by LTR_Finder.

```
LTR_Finder/ltr_finder contig_split/ctg1.fa -s LTR_Finder/source/tRNA/Athal-tRNAs.fa -w 2 >./ltr_finder_out/ctg1.fa.out
```

[01.get_seq.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/03.LTR_insertion_time/01.get_seq.pl)

[02.align_3.5.LTR.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/03.LTR_insertion_time/02.align_3.5_LTR.pl)

[03.merge_distmat_result.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/03.LTR_insertion_time/03.merge_distmat_result.pl)

[04.alculate_insertion_time.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/03.LTR_insertion_time/04.alculate_insertion_time.pl)

```
perl 01.get_seq.pl Gemome.fa

perl 02.align_3.5_LTR.pl | sh

perl 03.merge_distmat_result.pl

perl 04.alculate_insertion_time.pl 03.merge_distmat_result.pl.out
```
