#include <string>

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "aligns/KmerAlignCore.h"

#include <math.h>
#include "analysis/KmerTable.h"
#include "analysis/NonRedKmerTable.h"
#include "util/mutil.h"

void PrintSeq(const DNAVector & d) {
  int i;
  for (i=0; i<d.isize(); i++) {
    cout << d[i];
    if ((i+1) % 80 == 0)
      cout << endl;
  }
  cout << endl;
}


class Pool
{
public:
  Pool() {}

  void add(int i) {m_index.push_back(i);}

  int size() const {return m_index.isize();}
  int get(int i) const {return m_index[i];}

private:
  svec<int> m_index;
};

float compute_entropy(const DNAVector & kmer) {

  map<char,int> char_map;

  for (unsigned int i = 0; i < (unsigned int)kmer.isize(); i++) {
	
	char c = kmer[i];
	char_map[c]++;
  }

  float entropy = 0;
  
  char nucs[] = { 'G', 'A', 'T', 'C' };

  for (unsigned int i = 0; i < 4; i++) {
	
	char nuc = nucs[i];
	
	int count = char_map[nuc];
	
	float prob = (float)count / kmer.isize();
	
	if (prob > 0) {
	  float val = prob * log(1/prob)/log(2);
	  entropy += val;
	}
  }

  return(entropy);
}

bool IsSimple(const DNAVector & d) 
{
  int i;
  int k = 0;
  for (i=2; i<d.isize(); i++) {
    if (d[i] == d[i-2])
      k++;
  }
  double ratio = (double)k/(double)d.isize();
  if (ratio > 0.6)
    return true;
  else
    return false;


  /*
  double ent = compute_entropy(d);
  //cout << "Entropy: " << ent << endl;
  if (ent < 1.5)
    return true;
  else
  return false;*/
  /*  int i;
  int a=0;
  int c=0;
  int g=0;
  int t=0;

  for (i=0; i<d.isize(); i++) {
    if (d[i] == 'A')
      a++;
    if (d[i] == 'C')
      c++;
    if (d[i] == 'G')
      g++;
    if (d[i] == 'T')
      t++;
  }

  int max = (int)(0.75*(double)d.isize());
  if (a > max || c > max || g > max || t > max)
    return true;

  return false;
  */
}


bool IsShadow(const DNAVector & a, const DNAVector & b, int startA, int startB, int k) 
{
  //return false;


  int i;
  int n = 0;
  int nn = 0;
  int last = -1;
  int len = 0;
  for (i=startA; i<a.isize(); i++) {
    int x = i-startA + startB;
    
    if (x >= b.isize())
      break;
    len++;
    if (a[i] != b[x]) {     
      //cout << "Mismatch @ " << i << " and " << x << endl;
      if (last >= 0) {
	int dist = i-last;
	if (x > 3 && i > 3 && a[i-1] != b[x-1] && a[i-2] != b[x-2])
	  break;
	if (dist == k+1) {
	  n++;
	} else {
	  nn++;
	}
      }
      last = i;
    } else {
      //cout << i << " and " << x << endl;
    }
  }

  int expect = (int)(0.9*(double(len/(k+1)-1)));
  //cout << "Len: " << len << " Expect: " << expect << " Observed: " << n << " diss: " << nn << endl;
  if (n >= expect && n > 4 && nn < n/5) {
    return true;
  }
  return false;

}


double Coverage(const string &s) 
{
  CMTokenizer tk;
  tk.AddDelimiter(";");
  tk.AddDelimiter(" ");
  tk.AddDelimiter("_");
  CMPtrStringList result;
  tk.Tokenize(result, s.c_str());
  
  //for (int i=0; i<result.length(); i++) {
  //  cout << "Token " << i << " " << *result(i) << endl;
  //}

  
  //cout << "Name: " << s << endl;

  if (result.length() < 3)
    return 1;
  double ret = atof(*result(1));

  //cout << "Value: " << ret << endl;
  if (ret < 1.)
    ret = 1.;

  return ret;
}

bool IsGoodCoverage(double a, double b) 
{
  //return true;

  double dev_a = sqrt(a);
  double dev_b = sqrt(b);

  double mean = (a+b)/2;
  
  if (dev_b > dev_a) {
    dev_a = dev_b;
    a = b;
  }

  //cout << "Coverage check: " << a << " mean: " << mean << " dev: " << dev_a << endl;

  if ((a - mean) < 10*dev_a)
    return true;
 
  double ratio = a/b;
  if (ratio < 1.)
    ratio = 1./ratio;

  if (ratio < 100.)
    return true;
  else
    return false;
}



class Welder
{
public:
  Welder(int k, int kk) {
    m_k = k;
    m_kk = kk;
    m_pTab = NULL;
  }

  void SetTable(NonRedKmerTable * p) {
    m_pTab = p;
  }

  void WeldableKmer(DNAVector & out, const DNAVector & a, int one, const DNAVector & b, int two) 
  {
    out.resize(m_kk);
    int flank = (m_kk - m_k)/2;
    
    int startA = one-flank;
    int stopA = one+m_k;
    int startB = two+m_k;
    int stopB = startB + flank;

    if (startA < 0 || stopB >= b.isize()) {
      out.resize(0);
      return;
    }

    int i;
    int j = 0;
    for (i=startA; i<stopA; i++) {
      out[j] = a[i];
      j++;
    }
    for (i=startB; i<stopB; i++) {
      out[j] = b[i];
      j++;
    }

  }
  
  
  bool Weldable(const DNAVector & a, int one, const DNAVector & b, int two, int thresh = 0) 
  {
    int i;
    DNAVector d;
    WeldableKmer(d, a, one, b, two);
    if (d.isize() == 0)
      return false;

    int count = m_pTab->GetCount(d, 0);

    //cout << "Solder: ";
    //for (i=0; i<d.isize(); i++)
    //  cout << d[i];
    //cout << endl;

    if (count > thresh)
      return true;
    else
      return false;
  }

private:
  NonRedKmerTable * m_pTab;

  int m_k;
  int m_kk;
};

void Add(vecDNAVector & all, DNAVector & add, int & counter) 
{
  if (counter >= all.isize()) {
    all.resize(all.isize() + 1000000);
    cout << "Resizing to : " << all.isize() << endl;
  }
  all[counter] = add;
  counter++;


}


int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-i","input fasta");
  commandArg<string> readStringCmmd("-r","read fasta");
  commandArg<bool> strandCmmd("-strand","strand specific", false);
  commandArg<int> kCmmd("-k","k-mer size for pooling", 24);
  commandArg<int> kkCmmd("-kk","k-mer size for welding", 48);
  //commandArg<string> bStringCmmd("-o","output fasta");
  commandLineParser P(argc,argv);
  P.SetDescription("Makes a graph out of a fasta");
  P.registerArg(aStringCmmd);
  P.registerArg(readStringCmmd);
  P.registerArg(strandCmmd);
  P.registerArg(kCmmd);
  P.registerArg(kkCmmd);
  
  //P.registerArg(bStringCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  bool bStrand = P.GetBoolValueFor(strandCmmd);
  string readString = P.GetStringValueFor(readStringCmmd);
  int k = P.GetIntValueFor(kCmmd);
  int kk = P.GetIntValueFor(kkCmmd);
  
  vecDNAVector dna;
  
  cout << "Reading file..." << endl;
  dna.Read(aString);
  cout << "done!" << endl;
 
  KmerAlignCore core;
  core.SetNumTables(2);
  TranslateBasesToNumberExact trans;
  core.SetTranslator(&trans);
  core.AddData(dna);
  core.SortAll();
 
  int i, j;

  if (k != 24) {
    cout << "The only size of k supported is 24! Exiting!" << endl;
    return -1;
  }
  int thresh = 1;


  bool bNoWeld = false;
  
  NonRedKmerTable kmers(kk);
  
 
  //string dbname;
  //TempFile(dbname, aString);
  //kmers.Read(dbname);
  //kmers.Add(seq);


  Welder weld(k, kk);

  svec<int> mapped;
  mapped.resize(dna.isize(), -1);
  svec<Pool> pool;

  int max_count = 10;

  int count = 0;

  vecDNAVector crossover;
  //crossover.reserve(100000);

  cout << "Collecting k-mers..." << endl;

  int counter = 0;
 

  for (i=0; i<dna.isize(); i++) {
    DNAVector & d = dna[i];
    if (i % 10000 == 0) {
      cout << "Processed: " << 100.*(double)i/(double)dna.isize() << " % of sequences. " << endl;
    }

    for (j=0; j<=d.isize()-k; j++) {
      DNAVector sub;
      sub.SetToSubOf(d, j, k);

      if (IsSimple(sub))
	continue;

       svec<KmerAlignCoreRecord> matchesFW, matchesRC;   

      core.GetMatches(matchesFW, sub);
      if (!bStrand) {
	sub.ReverseComplement();
	core.GetMatches(matchesRC, sub);
      }

      int x;

      for (x=0; x<matchesFW.isize(); x++) {
	int c = matchesFW[x].GetContig();
	if (c == i)
	  continue;

	DNAVector & dd = dna[c];	
	int start = matchesFW[x].GetPosition();

	DNAVector add;

	weld.WeldableKmer(add, d, j, dd, start);
	if (add.isize() > 0)
	  Add(crossover, add, counter);      

	weld.WeldableKmer(add, dd, start, d, j);
	if (add.isize() > 0) {
	  Add(crossover, add, counter);
	}
      }
      for (x=0; x<matchesRC.isize(); x++) {
	int c = matchesRC[x].GetContig();
	if (c == i)
	  continue;
	DNAVector dd = dna[c];
	dd.ReverseComplement();

	int start = dd.isize() - matchesRC[x].GetPosition() - k;
	DNAVector add;

	weld.WeldableKmer(add, d, j, dd, start); 
	if (add.isize() > 0)
	  Add(crossover, add, counter);	 


	weld.WeldableKmer(add, dd, start, d, j);
	if (add.isize() > 0)
	  Add(crossover, add, counter);	 
      }
    }
  }

  crossover.resize(counter);

  cout << "Setting up/sorting structure." << endl;
  kmers.SetUp(crossover);


  if (!bNoWeld) {
    cout << "Reading reads..." << endl;
    vecDNAVector seq;
    seq.Read(readString);
    
    cout << "Filtering reads..." << endl;
    kmers.AddData(seq);
    cout << "Done!" << endl;
  
    //cout << "Reading & filtering reads." << endl;
    //vecDNAVectorStream seqStream;
    //seqStream.ReadStream(readString);
    //kmers.AddData(seqStream);
    //cout << "Done!" << endl;
    
    
    weld.SetTable(&kmers);
  }

  //=================================================================
  cout << "Sequences: " << dna.isize() << endl;
  for (i=0; i<dna.isize(); i++) {

    int cutoff = 0;
    DNAVector & d = dna[i];
    if (d.isize() < k)
      continue;

    if (mapped[i] == -2)
      continue;

    Pool * pCurr = NULL;
    int poolID = -1;
    if (mapped[i] == -1) {
      Pool tmp;
      pool.push_back(tmp);
      pCurr = &pool[pool.isize()-1];
      poolID = pool.isize()-1;
      pCurr->add(i);
    } else {
      pCurr = &pool[mapped[i]];
      poolID = mapped[i];
    }

    //cout << "Added " << i << " to pool " << pool.isize()-1 << endl;
    mapped[i] = pool.isize()-1;

    for (j=0; j<=d.isize()-k; j++) {
      //if (cutoff > 200)
      //break;

      DNAVector sub;
      sub.SetToSubOf(d, j, k);

      svec<KmerAlignCoreRecord> matchesFW, matchesRC;   

      core.GetMatches(matchesFW, sub);
      if (!bStrand) {
	sub.ReverseComplement();
	core.GetMatches(matchesRC, sub);
      }
 
      if (IsSimple(sub) && matchesFW.isize() + matchesRC.isize() > 1)
	continue;      
      
      //if (matchesFW.isize() + matchesRC.isize() > 10)
      //continue;

      int x;

      //if (matchesFW.isize() + matchesRC.isize() >= max_count)
      //continue;
      double coverage = Coverage(dna.Name(i));

     
      int minCov = (int)(coverage / 25.);

      for (x=0; x<matchesFW.isize(); x++) {
	int c = matchesFW[x].GetContig();
	if (c == i)
	  continue;
	if (mapped[c] != -1) 
	  continue;
	double coverage_other = Coverage(dna.Name(c));


	DNAVector & dd = dna[c];	
	int start = matchesFW[x].GetPosition();


	if (!bNoWeld && !(weld.Weldable(d, j, dd, start, minCov) || weld.Weldable(dd, start, d, j, minCov))) {
	  //cout << "Reject, no sauter (fw)..." << endl;
	  continue;
	}
	

	if (IsShadow(d, dd, j, start, k)) {
	  cout << "Toasting shadow: " << dna.Name(c) << endl;
	  mapped[c] = -2;
	  //continue;
	} else {
	  if (!IsGoodCoverage(coverage, coverage_other)) {
	    cout << "Rejecting fw merge between " << dna.Name(i);
	    cout << " and " << dna.Name(c) << endl;
	    continue;
	  }
	  cout << "Accept (fw)!!" << endl;
	  pCurr->add(c);
	  mapped[c] = poolID;
	  cutoff++;
	  // (cutoff > 10)
	  //eak;
	}
	//cout << "Mapped sequence " << c << " to pool " << mapped[c] << " + ";
	//for (int y=0; y<sub.isize(); y++)
	//  cout << sub[y];
	//cout << endl;
      }
      for (x=0; x<matchesRC.isize(); x++) {
	int c = matchesRC[x].GetContig();
	if (mapped[c] != -1) 
	  continue;
	if (c == i)
	  continue;

	double coverage_other = Coverage(dna.Name(c));


	DNAVector dd = dna[c];
	dd.ReverseComplement();

	int start = dd.isize() - matchesRC[x].GetPosition() - k;

	if (!bNoWeld && !(weld.Weldable(d, j, dd, start, minCov) || weld.Weldable(dd, start, d, j, minCov))) {
	  //cout << "Reject, no sauter (rc)..." << endl;
	  //dd.ReverseComplement();
	  continue;
	}
	if (!IsGoodCoverage(coverage, coverage_other)) {
	  cout << "Rejecting rc merge between " << dna.Name(i);
	  cout << " and " << dna.Name(c) << endl;
	  continue;
	}
	cout << "Accept (rc)!!" << endl;
	pCurr->add(c);
	cutoff++;
	dna[c].ReverseComplement();
 	mapped[c] = poolID;
	//cout << "Mapped sequence " << dna.NameClean(c) << " to pool " << mapped[c] << " -" << endl;
      }
      
    }
  }

  for (i=0; i<pool.isize(); i++) {
    const Pool & p = pool[i];
    cout << "COMPONENT " << i << "\t" << p.size() << endl;
    for (j=0; j<p.size(); j++) {
      int z = p.get(j);
      cout << ">Component_" << i << " " << p.size() << " " << z << endl;
      PrintSeq(dna[z]);
    }
    cout << "END" << endl;
  }

  return 0;

}
  
