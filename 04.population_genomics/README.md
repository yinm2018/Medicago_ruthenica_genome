# Population Genetics Analysis of Medicago ruthenica

## (1) SNP Calling

### 1.Mapping Reads To The Reference Genome

#### Build index for reference genome

```
bwa index ref.fasta
samtools faidx ref.fasta
java -jar picard-tools-1.129/picard.jar CreateSequenceDictionary REFERENCE=ref.fasta OUTPUT=ref.dict
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

#### Get intervals of InDel with GATK

```
java -jar GenomeAnalysisTK.jar -nt 30 -R ref.fasta -T RealignerTargetCreator -o 03.realign/sample.realn.intervals -I 02.rehead/sample.rmdup.bam 2>03.realign/sample.realn.in
tervals.log
```

#### Realign on the InDel intervals

```
java -jar GenomeAnalysisTK.jar -R ref.fa -T IndelRealigner -targetIntervals 03.realign/sample.realn.intervals -o 03.realign/sample.realn.bam -I 2.rmdup/sample.rmdup.bam 2>0
3.realign/sample.realn.bam.log
```

### 2.Variation Detection

#### bam2gvcf

```
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R reference.fa -I 03.realign/sample.realn.bam -nct 15 -ERC GVCF -o 01.gvcf/sample.gvcf.gz -variant_index_type LINEAR -var
iant_index_parameter 128000
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

### 3.SNP Filter

#### SNP

```
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference.fa -V Pop.SNP.vcf.gz --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o Pop.HDflt.SNP.vcf.gz
```

#### InDel

```
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference.fa -V Pop.INDEL.vcf.gz --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter" -o Pop.HDflt.INDEL.vcf.gz
```

#### remove hdfilter

[01.hardfit.fix.indel.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/tree/main/04.population_genomics/01.hardfit.fix.indel.pl)
```
perl 01.hardfit.fix.indel.pl Pop.HDflt.INDEL.vcf.gz Pop.HDflted.INDEL.vcf.gz
```

[01.hardfit.fix.snp.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/tree/main/04.population_genomics/01.hardfit.fix.snp.pl)

```
perl 01.hardfit.fix.snp.pl Pop.HDflt.SNP.vcf.gz Pop.HDflted.SNP.vcf.gz
```

#### filtering SNPs within 5bp upstream & downstream of InDels

[02.flt.indel.5bp.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/tree/main/04.population_genomics/02.flt.indel.5bp.pl)

```
perl 02.flt.indel.5bp.pl -vcf Pop.HDflted.SNP.vcf.gz -indel Pop.HDflted.INDEL.vcf.gz -p 10 | gzip -c > Pop.HDflted.flt.5bp.SNP.vcf.gz
```

#### vcftools filter

```
vcftools --gzvcf Pop.HDflted.flt.5bp.SNP.vcf.gz --maf 0.05 --min-meanDP 3 --max-meanDP 24 --max-missing 0.7 --hwe 0.001 --minGQ 20 --recode --recode-INFO-all --out Pop.final.SNP.vcf
```


## (2) Admixture

```
vcftools --vcf Pop.final.SNP.vcf --plink --out Pop.final.SNP.vcf
plink --noweb --ped Pop.final.SNP.vcf.ped --map Pop.final.SNP.vcf.map --recode 12 --out admixture/Pop.final.SNP.vcf.extract
cd admixture
admixture --cv -j30 -B[100] Pop.final.SNP.vcf.extract.ped Pop_number > Pop.final.SNP.vcf.extract.log.out
```

## (3) Pop_history

### [psmc](https://github.com/lh3/psmc)

```
samtools mpileup -C50 -uf ../Medru.chr.fa lh2.realn.bam | bcftools call -c | vcfutils.pl vcf2fq -d 3 -D 20 | gzip > lh2.diploid.fq.gz
fq2psmcfa -q20 lh2.diploid.fq.gz >lh2.diploid.psmcfa
splitfa lh2.diploid.psmcfa >lh2.split.psmcfa
psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o ./01.result/lh2.round.1.psmc lh2.split.psmcfa
...
psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o ./01.result/lh2.round.100.psmc lh2.split.psmcfa
```

### [smcpp](https://github.com/popgenmethods/smcpp)

```
smc++ vcf2smc Pop.final.SNP.vcf.gz Pop1/chr1.smc.gz chr1 Pop1:S1,S2,...
smc++ vcf2smc Pop.final.SNP.vcf.gz Pop2/chr1.smc.gz chr1 Pop2:S1,S2,...
smc++ estimate -o Pop1.analysis/ 1.05e-8 Pop1/chr*.smc.gz
smc++ estimate -o Pop2.analysis/ 1.05e-8 Pop2/chr*.smc.gz
smc++ split -o split Pop1.analysis/model.final.json Pop2.analysis/model.final.json data/*.smc.gz
```
