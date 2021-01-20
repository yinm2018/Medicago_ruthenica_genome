# Genome assembly of Medicago ruthenica

The genome assembly of Medicago ruthenica 

## (1) Contig Assembly

### Canu (version 1.8)

The cleaned PacBio long reads were corrected for subsequent contig assembly by Canu.  

```
canu -correct -p GS180050-03 -d Canu –pacbio merged.subreads.fasta genomeSize=902970000 corOutCoverage=80 saveOverlaps=false
```

### [SmartDenovo](https://github.com/ruanjue/smartdenovo)

```
perl smartdenovo.pl –k 16 –J 5000 –p GS180050-03
```

### FALCON (version 3.1)

```
fc_run  fc_run.cfg

fc_run.cfg

[General]
input_fofn=input.fofn
input_type=raw
pa_DBdust_option=
pa_fasta_filter_option=pass
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
genome_size= 919260000
seed_coverage=30
length_cutoff=-1
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

The two version assembly results from FALCON and SmartDenove were integrated by Quickmerge.

```
quickmerge -hco 5.0 -c 1.5 -lm 5000 -l 300000
nucmer -l 100 -p out self_oneline.fa hybrid_oneline.fa
delta-filter -i 95 -r -q out.delta >out.rq.delta
quickmerge -d out.rq.delta -q hybrid_oneline.fa -r self_oneline.fa -hco 5.0 -c 1.5 -lm 5000 -l 10000
```


## (2) The Polish of Contig-scale Assembly

The contig-scale pre-assembly was polished by Illumina and PacBio clean subreads.

### blasr (version5.3.2)

```
blasr merged.subreads.fasta genome.fa --out blasr.bam --bam --bestn 5 --minMatch 18 --nproc 4 --minSubreadLength 1000 --minAlnLength 500 --minPctSimila rity 70 --minPctAccu
racy70 --hitPolicy randombest --randomSeed 1
```

### arrow (version2.2.2)

```
samtools faidx genome.fa
arrow --algorithm=arrow -v -j8 blasr.bam -r merged.fasta -o genome.fa
```

### bwa(version 0.7.9a), samtools (version 0.1.19)

```
bwa index -p genome genome.fa
bwa mem mem -M -k 30 genome read_1.fq.gz read_2.fq.gz | samtools sort -m 1G -o align.bam
```

### Pilon (version 1.22)

```
java -Xmx150G -jar pilon-1.22.jar --genome draft.fa --frags align.bam
```


## (3) Chromosome-scale assembly based on Hi-C data

### BWA (version 0.7.9a)

```
bwa index draft.genome
```

### Juicer

```
juicer/misc/generate_site_positions.py MboI your_species_name draft.genome
juicer/scripts/juicer.sh -g your_species_name -s MboI -z draft.genome -p assembly -y restriction_sites_path/Mru_MboI.txt -D juicer -d your_work_path -t 30 -S early
```

### LANCHESIS

The pseudochromsomes was constructed by LANCHESIS with default parameters.
[https://github.com/shendurelab/LACHESIS](https://github.com/shendurelab/LACHESIS)

### Juicebox

The visual error correction for chromsome-scale pre-assembly was accomplished by Juicebox.
[https://github.com/aidenlab/Juicebox](https://github.com/aidenlab/Juicebox)


## (4) Genome Statistics

The general imformation statistics for the genome.

- Genome Size  [00.statistics.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/00.genome.statistic/00.statistic.pl)
- Scaffold Length [00.statistics.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/00.genome.statistic/00.statistic.pl)
- N50 [02.statistic_N50.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/00.genome.statistic/02.statistic_N50.pl)
- GC content [01.statistics_GC_content.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/00.genome.statistic/01.statistic_GC_content.pl)
