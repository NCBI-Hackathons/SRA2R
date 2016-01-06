#include <math.h>
#include <iostream>

#include <Rcpp.h>


extern "C" {
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>
}
  
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
SEXP reads(Rcpp::String acc, int n, SEXP lkup) {
  ReadCollection run = ncbi::NGS::openReadCollection ( acc );
  ReadIterator rgi = run.getReads( Read::all );
  
  SEXP r_ans_width, r_ans;
  XVectorList_holder r_ans_holder;
  
  PROTECT(r_ans_width = IntegerVector(1,1000));
  PROTECT(r_ans = alloc_XRawList("DNAStringSet", "DNAString", r_ans_width));
  
  vector<std::string> out;
  
  for(int i = 0; rgi.nextRead() & (i < 3) ; i++) {
    cout << i;
    while ( rgi.nextFragment() ) {
      std::string str1 = rgi.getFragmentBases().toString();
      const char * abc = str1.c_str();
      cout << abc; 
      Chars_holder r_ans_elt_holder = get_elt_from_XRawList_holder(&r_ans_holder, i);
      Ocopy_bytes_to_i1i2_with_lkup(0, r_ans_elt_holder.length - 1,
                                    (char *)r_ans_elt_holder.ptr, r_ans_elt_holder.length,
                                    abc, str1.length(),
                                    INTEGER(lkup), LENGTH(lkup));
    }
  }
  UNPROTECT(2);
  return r_ans;
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
//' require(Biostrings)
//' lkup = get_seqtype_conversion_lookup('B','DNA')
//' x = SRA2R:::read1('SRR000123',1000,lkup)
//' x
// [[Rcpp::export]]
SEXP read1(SEXP classname, SEXP acc, SEXP lkup) {
  IntegerVector width = IntegerVector::create(4,4,4,4);
  int vv[4] = { 4,4,4,4 };
  std::vector<int> v(&vv[0], &vv[0]+4);
  IntegerVector width2 = wrap(v);
  SEXP r_ans_width, r_ans;
  XVectorList_holder r_ans_holder;
  
  PROTECT(r_ans = alloc_XRawList("BStringSet", "BString", width2));
  r_ans_holder = hold_XVectorList(r_ans);

  for( int i = 0; i < 4; i++) {
    Chars_holder r_ans_elt_holder =
      get_elt_from_XRawList_holder(&r_ans_holder, i);
    memcpy((char *) r_ans_elt_holder.ptr, "ACTG", INTEGER(width)[i] * 
      sizeof(char));
  }
  
  UNPROTECT(1);
  return r_ans;
}

