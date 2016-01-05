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

//' Get Fastq in the read collection
//' 
//' Fastq strings (without the quality) are returned as a List of strings
//' @param acc An accession or a path to an actual SRA file (with .sra suffix)
//' @return List of reads
//' @export
//' @examples
//' getFastqCount('SRR000123')
// [[Rcpp::export]]
vector<string> getFastq(Rcpp::String acc, int n) {
  // open requested accession using SRA implementation of the API
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );
  //Rcpp::String run_name = run.getName ();
  
  // compute window to iterate through
  long MAX_ROW = run.getReadCount (); 
  
  //start iterator on reads
  ReadIterator it = run.getReadRange ( 1, MAX_ROW, Read::all );
  
  long i;
  for ( i = 0; it.nextRead (); ++ i )
  {
    cout << it.getReadId();
    
    //iterate through fragments
    while ( it.nextFragment () ){
      cout << '\t' <<  it.getFragmentBases ();
      cout << '\t' <<  it.getFragmentQualities ();
    }
    
    cout << '\n';
  }
  
  cerr << "Read " << i << " spots for " << run_name << '\n';
}


