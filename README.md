# SRA2R
SRA2R, a package to import SRA data directly into R

http://ec2-52-90-90-120.compute-1.amazonaws.com/

```
chmod 400 ncbi_sra2r.pem
ssh -i 'ncbi_sra2r.pem' ubuntu@ec2-52-90-90-120.compute-1.amazonaws.com
su - userXXX
```

```
download.file('http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.7/sratoolkit.2.5.7-ubuntu64.tar.gz',destfile='~/sratoolkit.2.5.7-ubuntu64.tar.gz')
untar('~/sratoolkit.2.5.7-ubuntu64.tar.gz',compressed=TRUE)
```

http://r-pkgs.had.co.nz/

## sra-stat example (http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=sra-stat)

```
(known good but large)
sra-stat --quick --xml SRR390728
(smaller file)
sra-stat --quick --xml SRR2971307
```

## sra-pileup (http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=sra-pileup)

Command-line

```
sra-pileup -r chr20:1530960-1540960 SRR2971307
```

In R

```
x = read.delim(text = system('~/sratoolkit.2.5.7-ubuntu64/bin/sra-pileup -r chr20:1530960-1540960 SRR2971307',intern=TRUE))
```

