
#include <ncbi-vdb/NGS.hpp>
#include <ngs-bam/ngs-bam.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/Reference.hpp>
#include <ngs/AlignmentIterator.hpp>
#include <ngs/Alignment.hpp>

#include <math.h>
#include <iostream>

using namespace ngs;
using namespace std;

class AlignSliceTest
{
public:
  
  static void run_common ( ReadCollection & run, String refname, int start, int stop )
  {
    
    String run_name = run.getName ();
    
    // get requested reference
    Reference ref = run.getReference ( refname );
    
    // start iterator on requested range
    long count = stop - start + 1;
    AlignmentIterator it = ref.getAlignmentSlice ( start, count, Alignment::primaryAlignment );
    
    long i;
    for ( i = 0; it . nextAlignment (); ++ i )
    {
      cout         << it.getReadId ()
                   << '\t' << it.getReferenceSpec ()
                   << '\t' << it.getAlignmentPosition ()
                   << '\t' << it.getLongCigar ( false )         // unclipped
                   << '\t' << it.getAlignedFragmentBases ()
                   << '\n';
    }
    
    cerr << "Read " <<  i <<  " alignments for " <<  run_name << '\n';
  }
  
  static void run_csra ( String acc, String refname, int start, int stop )
  {
    // open requested accession using SRA implementation of the API
    ReadCollection run = ncbi::NGS::openReadCollection ( acc );
    run_common ( run, refname, start, stop );
  }
  
  static void run_bam ( String acc, String refname, int start, int stop )
  {
    // open requested accession using example BAM implementation of the API
    ReadCollection run = NGS_BAM::openReadCollection ( acc );
    run_common ( run, refname, start, stop );
  }
  
  static void run ( String acc, String refname, int start, int stop )
  {
    size_t dot = acc . find_last_of ( '.' );
    if ( dot != string :: npos )
    {
      String extension = acc . substr ( dot );
      if ( extension == ".bam" || extension == ".BAM" )
      {
        run_bam ( acc, refname, start, stop );
        return;
      }
    }
    
    run_csra ( acc, refname, start, stop );
  }
  
  
};

int main ( int argc, char const *argv[] )
{
  if ( argc != 5 )
  {
    cerr << "Usage: AlignSliceTest accession reference start stop\n";
  }
  else try
  {
    ncbi::NGS::setAppVersionString ( "AlignSliceTest.1.1.0" );
    AlignSliceTest::run ( argv[1], argv[2], atoi ( argv[3] ), atoi ( argv[4] ) );
    return 0;
  }
  catch ( ErrorMsg & x )
  {
    cerr <<  x.toString () << '\n';
  }
  catch ( exception & x )
  {
    cerr <<  x.what () << '\n';
  }
  catch ( ... )
  {
    cerr <<  "unknown exception\n";
  }
  
  return 10;
}