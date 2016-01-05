# SRA2R
SRA2R, a package to import SRA data directly into R

## The AWS Instance Details

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

### git setup for github

Before using git with github, you'll want to do this on the AWS instance.

```
ssh-keygen
cat ~/.ssh/id_rsa.pub
```

copy the contents of that string to the github settings/ssh keys location

```
git config --global user.name 'My Name'
git config --global user.email 'myemail@email.com'
```

### Developing on the package

Checkout the package using git (or Rstudio) and change the working directory to the 
SRA2R directory (with the DESCRIPTION file in it).

```
install.packages('devtools')
devtools::document()
devtools::load_all()
```

## SRA ToolKit examples

### sra-stat example 

http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=sra-stat

```
(known good but large)
sra-stat --quick --xml SRR390728
(smaller file)
sra-stat --quick --xml SRR2971307
(small and no alignment)
sra-stat --quick --xml ERR1162649
```

### sra-pileup

http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=sra-pileup

Command-line

```
sra-pileup -r chr20:1530960-1540960 SRR2971307
```

In R

```
x = read.delim(text = system('~/sratoolkit.2.5.7-ubuntu64/bin/sra-pileup -r chr20:1530960-1540960 SRR2971307',intern=TRUE))
```


## ncbi ngs SDK details

- /usr/include/ngs (interfaces for C++ ngs)
- /usr/include/ncbi-vdb (NGS.hpp)
- /usr/local/share/doc/ngs (javadoc)
- LD_LIBRARY_PATH = /usr/local/ngs/ngs-sdk/lib64:/usr/local/ncbi/ncbi-vdb/lib64:


## R and Rcpp documentation of interest

- http://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf
- http://adv-r.had.co.nz/Rcpp.html
- http://r-pkgs.had.co.nz/
- http://statr.me/rcpp-note/index.html
- https://cran.r-project.org/web/packages/Rcpp/vignettes/


