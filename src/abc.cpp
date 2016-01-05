#include <Rcpp.h>
#include <ncbi-vdb/NGS.hpp>
#include <ngs-bam/ngs-bam.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>


#include <math.h>
#include <iostream>

using namespace ngs;
using namespace std;
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//s

//' The readCount in the read collection.
//'
//' This simply returns the full read count.
//'
//' @param acc An accession or a path to an actual SRA file (with .sra suffix)
//' @return the number of reads in the collection
//' @export
//' @examples
//' readCount('SRR000123')
// [[Rcpp::export]]
long readCount(Rcpp::String acc) {
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );
  return run.getReadCount();
}

//' The reads in the read collection.
//'
//' This simply returns the full read count.
//'
//' @param acc An accession or a path to an actual SRA file (with .sra suffix)
//' @param n The number of reads to return
//' @return the reads in the collection
//' @export
//' @examples
//' reads('SRR000123')
// [[Rcpp::export]]
CharacterVector reads(Rcpp::String acc, int n) {
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );
  ReadIterator rgi = run.getReads( Read::all );
  CharacterVector out(n);
  for(int i = 0; rgi.nextRead() & ( i < n ) ; i++) {
    while ( rgi.nextFragment() ) {
      out[i] = rgi.getFragmentBases().toString();
    }
  }
  return out;
}

