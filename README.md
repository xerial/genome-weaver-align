
## Install
```
$ make install
```

Default installation folder is ```$HOME/local/bin/genome-weaver```

## Building BWT index 
```
$ JVM_OPT="-Xmx=32g" genome-weaver bwt hg19.fa 
```

## Single-end alignment
```
$ genome-weaver align -r hg19.fa (fastq file)  > (sam file)
```

## Paired-end alignment
 (soon)
 
## Configure yoru git
```
$ git confgi core.eol lf
$ git config core.autocrlf input
```

