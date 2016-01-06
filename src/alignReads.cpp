#include <Rcpp.h>
#include <ncbi-vdb/NGS.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>

#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace ngs;
using namespace std;


//' Read alignment
//'
//' This returns the aligned reads.
//'
//' @param acc An accession or a path to an actual SRA file (with .sra suffix)
//' @param ref The reference name 
//' @param start Start position (inclusive)
//' @param stop End position (inclusive)
//' @return aligned reads
//' @export
//' @examples
//' alignReadsWithRegion('SRR789392','NC_000020.10', 62926240, 62958722)
//' alignReadsWithRegion('SRR789334','NC_000020.10', 62926240, 62958722)
// [[Rcpp::export]]
Rcpp::List alignReadsWithRegion(Rcpp::String acc, Rcpp::String refname, int start, int stop ) {
  try{
    ReadCollection   run = ncbi::NGS::openReadCollection ( acc );
    
    //testing whether there is alignment
    try {
      long alignmentCount = run.getAlignmentCount();
      if (alignmentCount==0) {
        std::printf("no aligned reads availabe"); 
      }
    } catch (ngs::ErrorMsg ngsErr){
      forward_exception_to_r(ngsErr);
      return -1;
    }
    
    try {
      if (!run.hasReference ( refname )) {
        std::string errorAndRefNames = "The accession id "+ string(acc) +" does not have the reference " +  
          string(refname) +  ". The options are:";
        ReferenceIterator refIter = run.getReferences();
        while( refIter.nextReference() ) {
          errorAndRefNames += " " + refIter.getCanonicalName();
        }
        throw std::range_error(errorAndRefNames); 
      }
    } catch (ngs::ErrorMsg ngsErr){
      forward_exception_to_r(ngsErr);
      return -1;
    }
    
    try {
      // get requested reference
      ngs::Reference ref = run.getReference ( refname );
      
      long referenceLength = ref.getLength();
      if (stop<start || stop>referenceLength || start<1) {
        throw std::range_error("wrong reference range, reference length = " + toString(referenceLength)); 
        return -1;
      } 
      
      long count = stop - start + 1;
      AlignmentIterator alignit = ref.getAlignmentSlice ( start, count, Alignment::primaryAlignment );
      vector<std::string> readID;
      vector<std::string> reference;
      vector<int> position;
      vector<std::string> longCigar;
      vector<int> quality;
      vector<int> alignLng;
      vector<std::string> bases; 
      
      while( alignit.nextAlignment() ) {
        readID.push_back(alignit.getReadId().toString());
        reference.push_back(alignit.getReferenceSpec());
        position.push_back(alignit.getAlignmentPosition());
        longCigar.push_back(alignit.getLongCigar( false ).toString());
        quality.push_back(alignit.getMappingQuality());
        alignLng.push_back( alignit.getAlignmentLength());
        bases.push_back(alignit.getAlignedFragmentBases().toString());
      }
      return List::create (
          _["readID"] = readID,
          _["referenceSpec"] = reference,
          _["position"] = position,
          _["longCigar"] = longCigar,
          _["mappingQuality"] = quality,
          _["alignementLength"] = alignLng,
          _["sequence"] = bases
      );
    } catch (ngs::ErrorMsg ngsErr){
      forward_exception_to_r(ngsErr);
      return -1;
    }
  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
    return -1;
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
    return -1;
  } //try ReadCollection run
}


//' Read alignment
//'
//' This returns the aligned reads.
//'
//' @param acc An accession or a path to an actual SRA file (with .sra suffix)
//' @return aligned reads
//' @export
//' @examples
//' alignReads('SRR789392')
//' alignReads('SRR789334')
// [[Rcpp::export]]
Rcpp::List alignReads(Rcpp::String acc) {
  Rcpp::String refname = "";
  
  try{
    ReadCollection run = ncbi::NGS::openReadCollection ( acc );
    //testing whether there is alignment
    try {
      long alignmentCount = run.getAlignmentCount();
      if (alignmentCount==0) {
        std::printf("no aligned reads availabe"); 
      }
    } catch (ngs::ErrorMsg ngsErr){
      forward_exception_to_r(ngsErr);
      return -1;
    }
    
    ReferenceIterator refIter = run.getReferences();
    
    vector<std::string> readID;
    vector<std::string> reference;
    vector<int> position;
    vector<std::string> longCigar;
    vector<int> quality;
    vector<int> alignLng;
    vector<std::string> bases; 
    
    while( refIter.nextReference() ) {
      refname = refIter.getCanonicalName();
      std::printf("%s is aligning... \n", string(refname).c_str());
      
      try {
        // get requested reference
        ngs::Reference ref = run.getReference ( refname );
        long referenceLength = ref.getLength();
        long start = 1;
        long stop = referenceLength;
        long count = stop - start + 1;
        AlignmentIterator alignit = ref.getAlignmentSlice ( start, count, Alignment::primaryAlignment );
        
        while( alignit.nextAlignment() ) {
          readID.push_back(alignit.getReadId().toString());
          reference.push_back(alignit.getReferenceSpec());
          position.push_back(alignit.getAlignmentPosition());
          longCigar.push_back(alignit.getLongCigar( false ).toString());
          quality.push_back(alignit.getMappingQuality());
          alignLng.push_back( alignit.getAlignmentLength());
          bases.push_back(alignit.getAlignedFragmentBases().toString());
        }

      } catch (ngs::ErrorMsg ngsErr){
        forward_exception_to_r(ngsErr);
        return -1;
      }
    } //  while( refIter.nextReference() )

    return List::create (
        _["readID"] = readID,
        _["referenceSpec"] = reference,
        _["position"] = position,
        _["longCigar"] = longCigar,
        _["mappingQuality"] = quality,
        _["alignementLength"] = alignLng,
        _["sequence"] = bases
    );
  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
    return -1;
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
    return -1;
  } //try ReadCollection run
  
}


