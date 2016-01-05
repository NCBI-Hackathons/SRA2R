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
//' @return the number of reads in the collection
//' @export
//' @examples
//' getReference('SRR390728')
// [[Rcpp::export]]

DataFrame getReference(Rcpp::String acc) {
  // open requested accession using SRA implementation of the API
   ReadCollection run = ncbi::NGS::openReadCollection ( acc );

   Rcpp::String run_name ( run.getName () );
  
   // get all references
   ReferenceIterator it ( run.getReferences () );
   vector<std::string> out;
   vector<std::string> out2;
   vector<long> out3;
   vector<long> out4;
   
   
       while( it.nextReference () ) {
         out.push_back(it.getCanonicalName()) ;
         out2.push_back(it.getCommonName()) ;
         out3.push_back(it.getLength()) ;
         out4.push_back(it.getAlignmentCount());
       }
       
       return DataFrame::create (
           _["CanonicalNames"] = out, _["CommonNames"] = out2, _["Length"] = out3, _["AlignmentCount"] = out4
       );
}
  


   