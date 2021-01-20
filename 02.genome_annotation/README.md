# Genome annotaion of Medicago ruthenica

## 1.Identification of repetitive elements

### TRF

```
trf407b.linux64 Genome.fa 2 7 7 80 10 50 2000 -d -h
```

### RepeatMasker

```
RepeatMasker/RepeatMasker -pa 30 -species "All" -nolow -norna -no_is -gff Genome.fa 2>&1 | tee 02.RepeatMasker.log
```

### RepeatProteinMask

```
RepeatMasker/RepeatProteinMask -engine abblast -noLowSimple -pvalue 1e-04 Genome.fa 2>&1 | tee 02.RepeatProteinMask.log
```

### RepeatModeler

[run.repeatmodeler.pl](https://github.com/yinm2018/Medicago_ruthenica_genome/blob/main/02.genome_annotation/run.repeatmodeler.pl)

```
RepeatMasker/RepeatModeler/BuildDatabase -name Medru Genome.fa 2>&1 | tee 01.BuildDatabase.log

RepeatMasker/RepeatModeler/RepeatModeler -pa 25 -database Medru 2>&1 | tee 02.RepeatModeler.log; perl run.repeatmodeler.pl Genome.fa
```


## 2.Gene Prediction

### De novo prediction

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
genemark_hmm_euk_linux_64/ehmm/gmhmme3 -f gff3 -m genemark_hmm_euk_linux_64/ehmm/m_truncatula.mod -o  genemark.ctg1.gff ctg1.fa
```

#### SNAP

```
snap snap/HMM/A.thaliana.hmm ctg1.fa -gff > snap.ctg1.gff
```

### Homologous-based prediction

GeneWise was used to predict gene models
[http://www.ebi.ac.uk/~birney/wise2](http://www.ebi.ac.uk/~birney/wise2)

### Transcriptome-based prediction

[https://github.com/PASApipeline/PASApipeline/wiki](https://github.com/PASApipeline/PASApipeline/wiki)

```
PASApipeline-v2.3.3/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -r -g Genome.final.fa -t polished.transcripts.rename.fa.clean --ALIGNERS gmap,blat
```
