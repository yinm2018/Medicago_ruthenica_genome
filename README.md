# Workflow of Medicago ruthenica genomic analysis

## Genome Assembly

### Canu (version 1.8)
```
 canu -correct -p GS180050-03 -d Canu –pacbio merged.subreads.fasta genomeSize=902970000 corOutCoverage=80 saveOverlaps=false
```

### [SmartDenovo](https://github.com/ruanjue/smartdenovo)

```
perl smartdenovo.pl –k 16 –J 5000 –p GS180050-03
```

### FALCON(version 3.1)

```
fc_run	fc_run.cfg

fc_run.cfg 

[General]
input_fofn=input.fofn 
input_type=raw pa_DBdust_option= pa_fasta_filter_option=pass 
target=assembly 
skip_checks=False 
LA4Falcon_preload=false

#### Data Partitioning 
pa_DBsplit_option=-x500 -s300
ovlp_DBsplit_option=-x500 -s300

#### Repeat Masking 
pa_HPCTANmask_option= 
pa_REPmask_code=0,300;0,300;0,300

#### Pre-assembly 
genome_size= 919260000 seed_coverage=30 length_cutoff=-1
pa_HPCdaligner_option=-v -B128 -M24
pa_daligner_option=-e.8 -l2000 -k14 -h480 -w8 -s100
falcon_sense_option=--output-multi --min-idt 0.70 --min-cov 2 --max-n-read 1800 falcon_sense_greedy=False

####Pread overlapping
ovlp_daligner_option=-e.9 -l1000 -k24 -h1024 -w6 -s100
ovlp_HPCdaligner_option=-v -B128 -M24 -l500


####Final Assembly
overlap_filtering_setting=--max-diff 100 --max-cov 100 --min-cov 2
fc_ovlp_to_graph_option=
length_cutoff_pr=5000
```

### [Quickmerge](https://github.com/mahulchak/quickmerge)

```
quickmerge -hco 5.0 -c 1.5 -lm 5000 -l 300000
nucmer -l 100 -p out self_oneline.fa hybrid_oneline.fa
delta-filter -i 95 -r -q out.delta >out.rq.delta
quickmerge -d out.rq.delta -q hybrid_oneline.fa -r self_oneline.fa -hco 5.0 -c 1.5 -lm 5000 -l 10000
```

### blasr(version5.3.2)

```
blasr merged.subreads.fasta genome.fa --out blasr.bam --bam --bestn 5 --minMatch 18 --nproc 4 --minSubreadLength 1000 --minAlnLength 500 --minPctSimila rity 70 --minPctAccuracy70 --hitPolicy randombest --randomSeed 1
```

### arrow(version2.2.2)

```
samtools faidx genome.fa
arrow --algorithm=arrow -v -j8 blasr.bam -r merged.fasta -o genome.fa
```

### bwa(version 0.7.9a) samtools(version 0.1.19) 

```
bwa index -p genome genome.fa
bwa mem mem -M -k 30 genome read_1.fq.gz read_2.fq.gz | samtools sort -m 1G -o align.bam
```

### Pilon (version 1.22)

```
java -Xmx150G -jar pilon-1.22.jar --genome draft.fa --frags agn.bam
```

## Genome Statistics

- Genome Size
- N50
- Scaffold Length
- GC content

## Identification of repetitive elements

### TRF

```
trf407b.linux64 Genome.final.fa 2 7 7 80 10 50 2000 -d -h
```

### RepeatProteinMask
```
/RepeatMasker-lastest/RepeatMasker/RepeatProteinMask -engine abblast -noLowSimple -pvalue 1e-04 Genome.final.fa
```

### RepeatModeler
```
RepeatMasker/RepeatModeler/BuildDatabase -name Medru Genome.final.fa 2>&1 | tee 01.BuildDatabase.log

RepeatMasker/RepeatModeler/RepeatModeler -pa 25 -database Medru 2>&1 | tee 02.RepeatModeler.log; perl run.repeatmodeler.pl Genome.fa
```

## Genome Annotation

### De novo predection

#### Geneid

```
geneid/bin/geneid -3 -P geneid/param/arabidopsis.param.Aug_4_2004 ctg1.fa > geneid.ctg1.gff
```

#### Augustus
```
augustus.2.5.5/bin/augustus --species=arabidopsis ctg1.fa > augustus.ctg1.gff
```

#### GeneMarke-ET

```
genemark_hmm_euk_linux_64/ehmm/gmhmme3 -f gff3 -m ~user/sorftware ~user/software/genemark_hmm_euk_linux_64/ehmm/m_truncatula.mod -o  genemark.ctg1.gff ctg1.fa
```

#### SNAP

```
snap ~user/software/snap/HMM/A.thaliana.hmm ctg1.fa -gff > snap.ctg1.gff
```

### Homologous-based prediction

GeneWise was used to predict gene models
[http://www.ebi.ac.uk/~birney/wise2](http://www.ebi.ac.uk/~birney/wise2)

### transcriptome-based prediction

[https://github.com/PASApipeline/PASApipeline/wiki](https://github.com/PASApipeline/PASApipeline/wiki)

```
PASApipeline-v2.3.3/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -r -g Genome.final.fa -t polished.transcripts.rename.fa.clean --ALIGNERS gmap,blat
```

## LTR insertion time

```
mkdir contig_split ltr_finder_out

perl 00.contig_split.pl Genome.fa

~user/software/LTR_FINDER_parallel/bin/LTR_FINDER.x86_64-1.0.7/ltr_finder contig_split/ctg1.fa -s ~user/software/LTR_Finder/source/tRNA/Athal-tRNAs.fa -w 2 >./ltr_finder_out/ctg1.fa.out

perl 01.get_seq.pl gemome.fa 

perl 02.align_3.5_LTR.pl | sh

perl 03.merge_distmat_result.pl

perl 04.alculate_insertion_time.pl 03.merge_distmat_result.pl.out
```

## SNP calling

### Mapping reads to the reference genome

#### Build index for reference genome

```
bwa index ref.fasta
samtools faidx ref.fasta
java -jar ~user/software//picard/picard-tools-1.129/picard.jar CreateSequenceDictionary REFERENCE=ref.fasta OUTPUT=ref.dict
```

#### First mapping: Use the BWA MEM; default parameters (Pairwise sequenced file)

```
bwa mem -t 30 -R '@RG\tID:$samplename\tPL:illumina\tPU:illumina\tLB:$samplename\tSM:$samplename' ref.fasta sample.1.fq.gz sample.2.fq.gz | samtools sort -O bam -T /tmp/sample -o 01.bwa/sample.sort.bam
```

#### Merge the BAM files from different sequencing cell belongs to same sample

```
samtools merge sample.sort.bam sample.L1.bam sample.L2.bam ...
```

#### Remove the duplicated Reads

```
java -Xmx10g -jar picard.jar MarkDuplicates INPUT=01.bwa/sample.sort.bam OUTPUT=02.rmdup/sample.rmdup.bam METRICS_FILE=02.rmdup/sample.dup.txt REMOVE_DUPLICATES=true ; samtools index 02.rmdup/sample.rmdup.bam
```

#### Get intervals of INDEL with GATK

```
java -jar GenomeAnalysisTK.jar -nt 30 -R ref.fasta -T RealignerTargetCreator -o 03.realign/sample.realn.intervals -I 02.rehead/sample.rmdup.bam 2>03.realign/sample.realn.intervals.log
```

#### Realign on the INDEL intervals

```
java -jar GenomeAnalysisTK.jar -R ref.fa -T IndelRealigner -targetIntervals 03.realign/sample.realn.intervals -o 03.realign/sample.realn.bam -I 2.rmdup/sample.rmdup.bam 2>03.realign/sample.realn.bam.log
```

### Variantion detection

#### bam2gvcf

```
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R reference.fa -I 03.realign/sample.realn.bam -nct 15 -ERC GVCF -o 01.gvcf/sample.gvcf.gz -variant_index_type LINEAR -variant_index_parameter 128000
```

#### gvcf2vcf
```
java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R reference.fa -V sample1.gvcf.gz -V sample2.gvcf.gz -V (...) -o Pop.vcf.gz
```

#### vcf2snp,vcf2indel

````
java -jar GenomeAnalysisTK.jar -T SelectVariants -R reference.fa -V Pop.vcf.gz -selectType SNP -o Pop.SNP.vcf.gz

java -jar GenomeAnalysisTK.jar -T SelectVariants -R reference.fa -V Pop.vcf.gz -selectType INDEL -o Pop.INDEL.vcf.gz
````
### SNP FILTER

#### SNP

```
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference.fa -V Pop.SNP.vcf.gz --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o Pop.HDflt.SNP.vcf.gz
```

#### InDel

```
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference.fa -V Pop.INDEL.vcf.gz --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter" -o Pop.HDflt.INDEL.vcf.gz
```

#### remove hdfilter

```
perl 01.hardfit.fix.indel.pl Pop.HDflt.INDEL.vcf.gz Pop.HDflted.INDEL.vcf.gz

perl 01.hardfit.fix.snp.pl Pop.HDflt.SNP.vcf.gz Pop.HDflted.SNP.vcf.gz
```

#### filter indel within 5bp upstream & downstream of SNPs

```
perl 02.flt.indel.5bp.pl -vcf Pop.HDflted.SNP.vcf.gz -indel Pop.HDflted.INDEL.vcf.gz -p 10 | gzip -c > Pop.HDflted.flt.5bp.SNP.vcf.gz
```

#### vcftools filter

```
vcftools --gzvcf Pop.HDflted.flt.5bp.SNP.vcf.gz --maf 0.05 --min-meanDP 3 --max-meanDP 24 --max-missing 0.7 --hwe 0.001 --minGQ 20 --recode --recode-INFO-all --out Pop.final.SNP.vcf
```