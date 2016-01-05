
#include <Rcpp.h>
#include <ncbi-vdb/NGS.hpp>
#include <ngs-bam/ngs-bam.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>

#include <ngs/Reference.hpp>
#include <ngs/Alignment.hpp>
#include <ngs/PileupIterator.hpp>

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
//' @param refname Reference name for pile up
//' @param start An in for position of start of pileup
//' @param stop An in for position of stop of pileup
//' @param MinPileUpDepth Coverage required
//' @return the number of reads in the collection
//' @export
//' @examples
//' getPileUp('SRR390728')
// [[Rcpp::export]]
DataFrame getPileUp(Rcpp::String acc, Rcpp::String refname, int start = 1, int stop = 0, int MinPileUpDepth = 0 ) {
  
  // open requested accession using SRA implementation of the API
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );


  // get requested reference
  ngs::Reference ref = run.getReference ( refname );
  if ( start == 1 && stop == 0 ){
    stop =  ref.getLength();
  }
  
  
  // start iterator on requested range
  long count = stop - start + 1;
  PileupIterator it = ref.getPileupSlice ( start-1 /*0-based*/, count);
  
  vector<std::string> RefSpec;
  vector<long> RefPos;
  vector<char> RefBase;
  vector<long> PileDepth;

  while ( it.nextPileup ())
  {
 if ( it.getPileupDepth () >= MinPileUpDepth ) { 
        RefSpec.push_back(it.getReferenceSpec ());
        RefPos.push_back( it.getReferencePosition () + 1 );
        RefBase.push_back( it.getReferenceBase () );
        PileDepth.push_back(it.getPileupDepth ( ) );

 } 
  }
  
  return DataFrame::create (
      _["ReferenceSpec"] = RefSpec, _["ReferencePosition"] = RefPos, _["ReferenceBase"] = RefBase, _["PileupDepth"] = PileDepth
    
  );
  
}

