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
#include <string>
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
//
//' The readCount in the read collection.
//'
//' This simply returns the full read count.
//' @author Nick Bernstein
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
DataFrame getPileUp(Rcpp::String acc, Rcpp::String refname, int start = 1, int stop = 0, int MinPileUpDepth = 0, bool Quality = false ) {
  
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
  vector<char> AlignedQuality;
  vector<std::string>  AllAlignedQuality;
  vector<char> AlignedBases;
  vector<std::string>  AllAlignedBases;
  vector<std::string> AllAcc;
  while ( it.nextPileup ())
  {
 if ( it.getPileupDepth () >= MinPileUpDepth ) { 
        RefSpec.push_back(it.getReferenceSpec ());
        RefPos.push_back( it.getReferencePosition () + 1 );
        RefBase.push_back( it.getReferenceBase () );
        PileDepth.push_back(it.getPileupDepth ( ) );
        AllAcc.push_back(acc);
        
        AlignedQuality.clear();
        AlignedBases.clear();
        std::string base;          
          
        while ( it.nextPileupEvent() ){
          
        if (Quality){
          AlignedQuality.push_back( it.getAlignmentQuality() );
         }
          
          PileupEvent::PileupEventType e = it.getEventType ();
            
            if(e & PileupEvent::alignment_start)
            {
              base += '^';
              base += (char) (it.getMappingQuality() + 33 );
            }
            if(e & PileupEvent::insertion)
            {
              base += '+';
              StringRef ibases= it.getInsertionBases();
              int c = ibases.size();
              char buf[64];
              if(e & PileupEvent::alignment_minus_strand)
              {
                char *b = buf + sprintf(buf,"%d",c);
                const char *s = ibases.data();
                for(int i=0; i<c;i++,b++,s++)
                {
                  *b=tolower(*s);
                }
                *b='\0';
              }
              else 
                sprintf(buf,"%d%.*s",c,c,ibases.data());
              base += buf;
            }
            if ( ( e & PileupEvent::alignment_minus_strand ) != 0 )
            {
              switch ( e & 7 )
              {
              case PileupEvent::deletion:
                base += '<';
                break;
              case PileupEvent::match:
                base += ',';
                break;
              case PileupEvent::mismatch:
                base += tolower(it.getAlignmentBase ());
                break;
              }
            }
            else
            {
              switch ( e & 7 )
              {
              case PileupEvent::deletion:
                base += '>';
                break;
              case PileupEvent::match:
                base += '.';
                break;
              case PileupEvent::mismatch:
                base += toupper(it.getAlignmentBase ());
                break;
              }
            }
            if(e & PileupEvent::alignment_stop)
            {
              base += '$';
            }
            
        }
        AllAlignedBases.push_back( base);  
        
      if ( Quality ) {
        std::string str(AlignedQuality.begin(),AlignedQuality.end());
        AllAlignedQuality.push_back( str );  
       }
 }
  }
if ( Quality ) {
  return DataFrame::create (
      _["AccensionNumber"] = AllAcc, _["ReferenceSpec"] = RefSpec, _["ReferencePosition"] = RefPos, _["ReferenceBase"] = RefBase, _["PileupDepth"] = PileDepth , _["AllAlignedBases"] = AllAlignedBases,  _["AllAlignedQuality"] = AllAlignedQuality
  );
  }
  else{
    return DataFrame::create (
        _["ReferenceSpec"] = RefSpec, _["ReferencePosition"] = RefPos, _["ReferenceBase"] = RefBase, _["PileupDepth"] = PileDepth , _["AllAlignedBases"] = AllAlignedBases, _["AccensionNumber"] = AllAcc
    
   );
  }

  
}

//Ran 1/8/2016
//> system.time(getPileUp("SRR1596669",'21',1,20000000, Quality =  F))
//user  system elapsed 
//43.599   2.495  46.664 
//> system.time(getPileUp("SRR1596669",'21',1,20000000, Quality =  T))
//user  system elapsed 
//60.158   1.885  70.376 

//Ran 1/13/2016
//> system.time(getPileUp("SRR1596669",'21',1,20000000, Quality =  F))
//  user  system elapsed 
//  48.756   1.363  50.284 
//> system.time(getPileUp("SRR1596669",'21',1,20000000, Quality =  T))
//  user  system elapsed 
//  48.212   1.320  49.536





