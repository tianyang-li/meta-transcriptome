//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif

#include <string.h>

#include "base/CommandLineParser.h"

#include "aligns/KmerAlignCore.h"
#include "analysis/DNAVector.h"
#include <math.h>
#include "base/FileParser.h"
#include "analysis/KmerTable.h"

#include "analysis/CompMgr.h"

void SortPrint(FILE * pReads, svec<IDS> & ids, const vecDNAVector & seq) 
{
  long long i;

  Sort(ids);
  int lastID = -1;
  int id = -1;
  int start = -1;
  int edge = -1;
  int lastStart = -1;
  int lastEdge = -1;
  int ori;

  string line;
  char tmp[1024];

  int lastStartTemp = -1;
  int lastOri = 1;

  for (i=0; i<ids.lsize(); i++) {
    id = ids[i].ID();
    ori = ids[i].Ori();
    start = ids[i].Start();
    edge = ids[i].Edge();
    //cout << id << "\t" << start << "\t" << edge << endl;
    if (id != lastID
#ifndef NO_REVERSE_OUT	  
	|| ori != lastOri
#endif
	) {
      if (lastID != -1) {

	sprintf(tmp, "%d\t%d\t", lastStart, lastEdge);	
	line += tmp;
	if (lastStart > lastStartTemp) {
	  fprintf(pReads, "%s\t", line.c_str());
	  //const DNAVector &d = seq[lastID];

#ifndef NO_REVERSE_OUT	  
	  DNAVector d = seq[lastID];
	  if (lastOri == -1) {
	    d.ReverseComplement();
	    //cout << "Reversing" << endl;
	  } else {
	    //cout << "Forward" << endl;
	  }
#endif 

	  for (int j=0; j<d.isize(); j++) {
	    tmp[1] = 0;
	    tmp[0] = d[j];
	    fprintf(pReads, "%s", tmp);
	  }
	  
	  if (lastOri == -1)
	    fprintf(pReads, "\t-");
	  else
	    fprintf(pReads, "\t+");
	  fprintf(pReads, "\n");
	}
	//fprintf(pReads, "%d\t%d\n", lastStart, lastEdge);
	line = "";
      }
      //fprintf(pReads, "%s\t%d\t%d\t", seq.Name(id).c_str(), start, edge);
      sprintf(tmp, "%s\t%d\t%d\t", seq.Name(id).c_str(), start, edge);
      line = tmp;

      lastStartTemp = start;
    }
    lastID = id;
    lastStart = start;
    lastEdge = edge;
    lastOri = ori;

  }

  if (id != -1) {
    sprintf(tmp, "%d\t%d\t", start, edge);	
    line += tmp;
    if (lastStart > lastStartTemp) {
      fprintf(pReads, "%s\t", line.c_str());
      DNAVector d = seq[id];
      if (ori == -1)
	d.ReverseComplement();
      for (int j=0; j<d.isize(); j++) {
	tmp[1] = 0;
	tmp[0] = d[j];
	fprintf(pReads, "%s", tmp);
      }
      if (lastOri == 1)
	fprintf(pReads, "\t+");
      else
	fprintf(pReads, "\t-");
	
	
      fprintf(pReads, "\n");
    }
    //fprintf(pReads, "%d\t%d\n", start, edge);
  }
}

//========================================================================
//========================================================================
//========================================================================


bool Irregular(char l)
{
  if (l == 'A' || l == 'C' || l == 'G' || l == 'T')
    return false;
  //cout << "Irregular char: " << l << endl;
  return true;
}



string ReadsExt(const string & in) 
{
  char tmp[1024];
  strcpy(tmp, in.c_str());
  int n = strlen(tmp);

  
  for (int i=n-1; i>=0; i--) {
    if (n-i > 6) {
      break;
    }
    if (tmp[i] == '.') {
      tmp[i] = 0;
      string out = tmp;
      out += ".reads";
      return out;
    }
      
  }
  string out = in + ".reads";
  return out;
}

int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-i","read fasta file");
  commandArg<string> gStringCmmd("-g","graph file");
  commandArg<string> oStringCmmd("-o","graph output");
  commandArg<int> kCmmd("-k","kmer size", 24);
  commandArg<bool> strandCmmd("-strand","strand specific", false);
  //commandArg<bool> cCmmd("-nc","do not fully connected graph", false);

  commandLineParser P(argc,argv);
  P.SetDescription("Assembles k-mer sequences.");
  P.registerArg(aStringCmmd);
  P.registerArg(gStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(kCmmd);
  P.registerArg(strandCmmd);
  //P.registerArg(cCmmd);


  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string gString = P.GetStringValueFor(gStringCmmd);
  string oString = P.GetStringValueFor(oStringCmmd);
  bool bStrand = P.GetBoolValueFor(strandCmmd);
  int k = P.GetIntValueFor(kCmmd)+1;


  int i, j;
  //vecbasevector contigBases;
 
  vecDNAVector seq;
  seq.Read(aString);
  
  KmerSequence kmers(k, &seq);
  kmers.Add(seq);
  
  long long m = kmers.GetBoundValue();

  FlatFileParser parser;
  
  parser.Open(gString);

  FILE * pOut = fopen(oString.c_str(), "w");
  string reads = ReadsExt(oString);
  FILE * pReads = fopen(reads.c_str(), "w");

  svec<IDS> ids;
  ids.reserve(seq.isize());

  //string last;
  //int lastNode = -1;

  svec<char> first;

  while (parser.ParseLine()) {
    
    if (parser.GetItemCount() < 4) {
      fprintf(pOut, "%s\n", parser.Line().c_str());
      if (ids.lsize() > 0) {
	SortPrint(pReads, ids, seq);
	//UniqueSort(ids);
	//for (i=0; i<ids.isize(); i++) {
	//fprintf(pReads, "%s\n", seq.Name(ids[i]).c_str());
	//}
      }
      fprintf(pReads, "%s\n", parser.Line().c_str());
      ids.clear();
     
      first.clear();
      //first.reserve(50000);
      continue;
    }

    const string & s = parser.AsString(3);
    int node = parser.AsInt(0);
    int prevNode = parser.AsInt(1);

    const char * p2 = s.c_str();
    if (node >= first.isize())
      first.resize(node + 10000, 'N');
    first[node] = p2[0];
 
    //if (prevNode >= 0 && last == "" ) {
    //  cout << "Potential ERROR" << endl;
    //}
    //if(prevNode >= 0 && prevNode != lastNode) {
    //  cout << "Potential ERROR (2): prevNode = " << prevNode;
    //  cout << " lastNode=" << lastNode << endl;      
    //}
    //int edge = parser.AsInt(0);
    long long edge = prevNode;
    long long n1 = 0;
    long long n2 = 0;
    if (prevNode >= 0) {

      DNAVector sub;
      sub.resize(strlen(s.c_str())+1);
      const char * p = s.c_str();
      for (i=0; i<sub.isize()-1; i++)
	sub[i+1] = p[i];

      //const char * p2 = last.c_str();
      if (first[prevNode] == 'N')
	cout << "ERROR!! " << endl;
      sub[0] = first[prevNode];

      kmers.BasesToNumberCountPlus(ids, n1, sub, edge);

      if (!bStrand) {
	sub.ReverseComplement();
	long long from = ids.lsize();
	kmers.BasesToNumberCountPlus(ids, n2, sub, edge);

	if (n1 + n2 < 0x7FFFFFFF) {

	  for (long long x=from; x<ids.lsize(); x++) {
	    ids[x].SetOri(-1);
	    int len = seq[ids[x].ID()].isize();
	    int pos = ids[x].Start()+1;
	    //cout << "len=" << len << " pos=" << pos;
#ifndef NO_REVERSE_OUT	  
	    ids[x].SetStart(len-pos-k+1);
	    //cout << " new=" << len-pos-k+1 << endl;
#else
	    ids[x].SetStart(pos+1);
#endif
	  }
	} else {
	  cout << "WARNING: k-mer overflow, n=" << n1 + n2 << ". Discarding." << endl;
	  n1 = n2 = 0;
	  //ids.resize(0);	
	}
      }

      /*  if (s == "ATATCACAAAACAATCTTCATTCG") {
	for (int x = 0; x<sub.isize(); x++)
	  cout << sub[x];
	cout << endl;
		
	cout << "Count=" << n1 + n2 << endl;
	}*/
    }

    for (i=0; i<parser.GetItemCount(); i++) {
      if (i>0)
	fprintf(pOut, "\t");
      if (i == 2) {
	fprintf(pOut, "%d", (int)(n1+n2));
      } else {
	fprintf(pOut, "%s", parser.AsString(i).c_str());
      }
    }
    fprintf(pOut, "\n");
    //lastNode = node;
    //last = s;
  }
  if (ids.lsize() > 0) {
    SortPrint(pReads, ids, seq);
  }

  fclose(pOut);
  fclose(pReads);
  return 0;

}
