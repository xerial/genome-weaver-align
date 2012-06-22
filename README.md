
# Genome Weaver - A Toolkit for Genome Sciences

* align - FM-index based Aligner
* lens	- Data format reader (FASTA,FASTQ,BED,WIG,GFF,SAM/BAM,etc.)

## Install
```
$ make install
```

Default installation folder is ```$HOME/local/bin/genome-weaver```

### Building BWT index 
```
$ JVM_OPT="-Xmx=32g" genome-weaver bwt hg19.fa 
```

### Single-end alignment
```
$ genome-weaver align -r hg19.fa (fastq file)  > (sam file)
```

### Paired-end alignment
 (soon)
 
### Configure your Git
```
$ git confgi core.eol lf
$ git config core.autocrlf input
```

