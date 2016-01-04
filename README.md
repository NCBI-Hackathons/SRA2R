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
