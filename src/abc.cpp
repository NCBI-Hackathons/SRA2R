#include <math.h>
#include <iostream>

#include <R.h>
#include <Rcpp.h>
#include <Rinternals.h>

#include <ncbi-vdb/NGS.hpp>
#include <ngs-bam/ngs-bam.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>

extern "C" {
  #include "S4Vectors_interface.h"
  #include "IRanges_interface.h"
  #include "XVector_interface.h"
  #include "Biostrings_interface.h"
}

using namespace ngs;
using namespace std;
using namespace Rcpp;

// The function below just returns a single value from the API.
// I just need to match up the return type of the function
// to the return type of the API call.

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

// The function below demonstrates using a 
// std::vector template to accumulate reads (strings)
// and return them to R as an R list.

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
List reads(Rcpp::String acc, int n, SEXP lkup) {
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );
  ReadIterator rgi = run.getReads( Read::all );
  
  SEXP r_ans_width, r_ans;
  XVectorList_holder r_ans_holder;
  
  PROTECT(r_ans_width = IntegerVector(1000));
  PROTECT(r_ans = alloc_XRawList("DNAStringSet", "DNAString", r_ans_width));
  
  vector<std::string> out;
  for(int i = 0; rgi.nextRead() & (i < 990) ; i++) {
    cout << i;
    while ( rgi.nextFragment() ) {
      out.push_back(rgi.getFragmentBases().toString());
      Chars_holder r_ans_elt_holder = get_elt_from_XRawList_holder(&r_ans_holder, i);
      Ocopy_bytes_to_i1i2_with_lkup(0, r_ans_elt_holder.length - 1,
                                    (char *)r_ans_elt_holder.ptr, r_ans_elt_holder.length,
                                    "ACTG", 4,
                                    INTEGER(lkup), LENGTH(lkup));
    }
  }
  return List::create (
     _["reads"] = out
  );
}

