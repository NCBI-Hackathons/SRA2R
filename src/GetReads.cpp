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



//' The readCount in the read collection.
//'
//' This simply returns the full read count.
//'
//' @param acc An accession or a path to an actual SRA file (with .sra suffix)
//' @return the number of reads in the collection
//' @export
//' @examples
//' getFastqCount('SRR000123')
// [[Rcpp::export]]
long getFastqCount(Rcpp::String acc) {
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );
  long MAX_ROW = run.getReadCount ();
  return MAX_ROW;
}


//' The reads in the read collection.
//'
//' This returns the all reads.
//'
//' @param acc An accession or a path to an actual SRA file (with .sra suffix)
//' @param n The number of reads to return
//' @return the reads in the collection
//' @export
//' @examples
//' getFastqReads('SRR000123')
// [[Rcpp::export]]
CharacterVector getFastqReads(Rcpp::String acc) {
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );
  long MAX_ROW = run.getReadCount (); 
  
  ReadIterator rgi = run.getReads( Read::all );
  CharacterVector out(MAX_ROW);
  for(int i = 0; rgi.nextRead() & ( i < MAX_ROW ) ; i++) {
    while ( rgi.nextFragment() ) {
      out[i] = rgi.getFragmentBases().toString();
    }
  }
  return out;
}

//' The reads in the read collection.
//'
//' This returns the all reads.
//'
//' @param acc An accession or a path to an actual SRA file (with .sra suffix)
//' @param n The number of reads to return
//' @return the reads in the collection
//' @export
//' @examples
//' getFastqReadsWithQuality('SRR000123')
// [[Rcpp::export]]
Rcpp::List getFastqReadsWithQuality(Rcpp::String acc) {
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );
  long MAX_ROW = run.getReadCount (); 
  
  ReadIterator rgi = run.getReads( Read::all );
  CharacterVector reads(MAX_ROW);
  CharacterVector qualities(MAX_ROW);

  
  for(int i = 0; rgi.nextRead() & ( i < MAX_ROW ) ; i++) {
    while ( rgi.nextFragment() ) {
      reads[i] = rgi.getFragmentBases().toString();
      qualities[i] = rgi.getFragmentQualities().toString();
    }
  }
  
  List fastqRead = Rcpp::List::create(Rcpp::Named("read") = reads,
                                      Rcpp::Named("quality") = qualities);
  return fastqRead;
}


