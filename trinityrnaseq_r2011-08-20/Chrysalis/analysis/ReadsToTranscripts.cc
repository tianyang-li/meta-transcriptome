#include <string>

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "analysis/NonRedKmerTable.h"

#include "analysis/CompMgr.h"

class Assignment
{
public:
    Assignment() {
        m_comp = -1;
        m_read = -1;
    }
    
    void SetComp(int c) {
        m_comp = c;
    }
    void SetRead(int r) {
        m_read = r;
    }
    
    int Read() const {return m_read;}
    int Comp() const {return m_comp;}
    
    bool operator < (const Assignment & a) const {
        return (m_comp < a.m_comp);
    }
    
private:
    int m_comp;
    int m_read;
};



int main(int argc,char** argv)
{
    
    
    commandArg<string> aStringCmmd("-i","reads fasta");
    commandArg<string> fStringCmmd("-f","fasta input file (concatenated flat components)");
    commandArg<string> bStringCmmd("-o","output directory");
    commandArg<bool> strandCmmd("-strand","strand specific data", false);
    commandLineParser P(argc,argv);
    P.SetDescription("Assigns reads to graph components.");
    P.registerArg(aStringCmmd);
    P.registerArg(fStringCmmd);
    P.registerArg(bStringCmmd);
    P.registerArg(strandCmmd);
    
    P.parse();
    
    string aString = P.GetStringValueFor(aStringCmmd);
    string bString = P.GetStringValueFor(bStringCmmd);
    string fString = P.GetStringValueFor(fStringCmmd);
    bool bStrand = P.GetBoolValueFor(strandCmmd);
    
    vecDNAVector reads, dna;
    
    cout << "Reading file..." << endl;
    reads.Read(aString);
    dna.Read(fString);
    cout << "done!" << endl;
    //test.ReverseComplement();
    
    int i, j;
    
    int k = 25;
    
    ComponentFileMgr compMgr(bString);
    
    string readCountFile = bString + "/readcounts.out";
    
    
    NonRedKmerTable kt(k);
    kt.SetUp(dna, true);
    kt.SetAllCounts(-1);
    
    for (i=0; i<dna.isize(); i++) {
        const DNAVector & d = dna[i];
        for (j=0; j<=d.isize()-k; j++) {
            kt.SetCount(d, j, i);
        }
    }
    
    svec<Assignment> assign;
    assign.resize(reads.isize());
    
    for (i=0; i<assign.isize(); i++)
        assign[i].SetRead(i);
    
    
    for (i=0; i<reads.isize(); i++) {

        if (i % 1000 == 0) {
            cerr << "\r[" << i << "] reads processed, " << (float)i/(reads.isize())*100 << "% finished   ";
        }
        
        const DNAVector & d = reads[i];
        svec<int> comp;
        comp.reserve(4000);
        for (j=0; j<=d.isize()-k; j++) {
            int c = kt.GetCountReal(d, j);
            if (c >= 0)
                comp.push_back(c);
        }
        
        if (!bStrand) {
            DNAVector dd = d;
            dd.ReverseComplement();
            for (j=0; j<=dd.isize()-k; j++) {
                int c = kt.GetCountReal(dd, j);
                if (c >= 0)
                    comp.push_back(c);
            }
        }
        
        
        Sort(comp);
        int best = -1;
        int max = 0;
        int run = 0;
        for (j=1; j<comp.isize(); j++) {
            if (comp[j] != comp[j-1] || j+1 == comp.isize()) {
                if (run > max) {
                    max = run;
                    best = comp[j-1];
                }
                run = 0;
            } else {
                run++;
            }
        }
        //cout << "Read " << i << " maps to " << best << " with " << max << " k-mers" << endl;
        assign[i].SetComp(best);
        
    }
    
    Sort(assign);
    
    int last = -2;
    FILE * p = NULL;
    
    int readCount = 0;
    for (i=0; i<assign.isize(); i++) {
        // Print it out...
        if (assign[i].Comp() == -1)
            continue;
        
        if (assign[i].Comp() != last) {
            if (p != NULL) {
                fclose(p);
            }
            string name = compMgr.GetFileName(assign[i].Comp(), ".raw.fasta");
            cout << "Writing: " << name << endl;
            p = fopen(name.c_str(), "w");
            
            last = assign[i].Comp();
        }
        int r = assign[i].Read();
        fprintf(p, "%s\n", reads.Name(r).c_str());
        readCount++;
        const DNAVector & d = reads[r];
        for (j=0; j<d.isize(); j++)
            fprintf(p, "%c", d[j]);
        
        fprintf(p, "\n");
    }
    
    fclose(p);
    
    
    
    FILE * pReadCount = fopen(readCountFile.c_str(), "w");
    fprintf(pReadCount, "%d\n", readCount/2);
    fclose(pReadCount);
    
    return 0;
}

