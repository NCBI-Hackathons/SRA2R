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
//' 
/* getReferenceBases
 *  return sub-sequence bases for Reference
 *  "offset" is zero-based
 */

/* getLength
 *  returns the length of the reference sequence
 */
// [[Rcpp::export]]
List refBases(Rcpp::String acc) {
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );
  
  ReferenceIterator ri = run.getReferences() ;
  
  vector<std::string> rlist;
  
  while ( ri.nextReference() ) {
    rlist.push_back(ri.getReferenceBases(0));
  }
  return Rcpp::List::create(Rcpp::Named("referencebases") = rlist);
  
}


