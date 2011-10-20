
#include <string>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "analysis/TranscriptomeGraph.h"
#include "analysis/DNAVector.h"
#include "analysis/CompMgr.h"


void Execute(const char * command) {
    int ret = system(command);
    if (ret != 0) {
        cout << "COMMAND: " << command << endl;
        cout << "Died with exit code " << ret << endl;
        cout << "Exiting." << endl;
        exit(-1);
    }
    
}

bool Exists(const string & s) 
{
    FILE * p = fopen(s.c_str(), "r");
    if (p != NULL) {
        fclose(p);
        return true;
    }
    // cout << "FATAL ERROR: Could not open file for read: " << s << endl;
    // cout << "Please make sure to enter the correct file name(s). Exiting now." << endl;
    
    return false;
}


int main(int argc,char** argv)
{
    char execPath[512];
    strcpy(execPath, argv[0]);
    execPath[strlen(execPath)-9] = 0;
    string exec = execPath;
    cout << "Path to executables: " << exec << endl;
    
    if (exec == "")
        exec = "./";
    commandArg<string> iStringCmmd("-i","read fasta file");
    commandArg<string> iwormStringCmmd("-iworm","inchworm file", "");
    commandArg<string> oStringCmmd("-o","output directory");
    commandArg<string> butterflyCmmd("-butterfly","butterfly executable", "../Butterfly/Butterfly.jar");
    commandArg<bool> skipCmmd("-skip","skip initial 2 steps", false);
    commandArg<bool> strandCmmd("-strand","strand-specific data", false);
    commandArg<bool> nobreakCmmd("-nobreak","skip breaking", false);
    commandArg<int> minCmmd("-min","minimum sequence length", 300);
    commandArg<int> cpuCmmd("-cpu","number of CPUs to use", 10);
    commandArg<int> distCmmd("-dist","size of a read pair insert", 350);
    commandArg<int> minDumpLenCmmd("-min_all","skip components for which all seqs < min_all", 110);
    commandArg<bool> buttCmmd("-run_butterfly","runs butterfly locally", false);
    commandArg<int> iniKCmmd("-ini_k","value of k for pre-processing", 25);
    
    commandLineParser P(argc,argv);
    P.SetDescription("Assemble transcriptomes from reads.");
    P.registerArg(iStringCmmd);
    P.registerArg(iwormStringCmmd);
    P.registerArg(butterflyCmmd);
    P.registerArg(oStringCmmd);
    P.registerArg(skipCmmd);
    P.registerArg(strandCmmd);
    P.registerArg(nobreakCmmd);
    P.registerArg(minCmmd);
    P.registerArg(cpuCmmd);
    P.registerArg(buttCmmd);
    P.registerArg(distCmmd);
    P.registerArg(minDumpLenCmmd);
    P.registerArg(iniKCmmd);
    
    P.parse();
    
  


    string readString = P.GetStringValueFor(iStringCmmd);
    string iwormString = P.GetStringValueFor(iwormStringCmmd);
    string outDir = P.GetStringValueFor(oStringCmmd);
    string butterflyExec = P.GetStringValueFor(butterflyCmmd);
    
    bool bSkip = P.GetBoolValueFor(skipCmmd);
    bool bStrand = P.GetBoolValueFor(strandCmmd);
    int minDumpLen = P.GetIntValueFor(minDumpLenCmmd);
    int minLen = P.GetIntValueFor(minCmmd);
    int nCPU = P.GetIntValueFor(cpuCmmd);
    int pairDist = P.GetIntValueFor(distCmmd);
    int ini_k = P.GetIntValueFor(iniKCmmd);
    bool bBreak = true;
    if (P.GetBoolValueFor(nobreakCmmd))
        bBreak = false;
    
    bool bButt = P.GetBoolValueFor(buttCmmd);
    

    if (minDumpLen > minLen) {
        minDumpLen = minLen;
    }
    
    string command;
    
    command = "mkdir ";
    command += outDir;
    system(command.c_str());
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // alternative to running Inchworm 
    // (TranscriptomeFromVaryK is the predecessor to Inchworm, originally referred to as Ananas)
    string iniSeqs = outDir + "/sequences.out";
    if (iwormString == "") {
        command = exec + "TranscriptomeFromVaryK -k 25 -o "; 
        command += outDir + "/sequences.out -i ";
        command += readString;
        if (bStrand) {
            command += " -strand ";
        }
        
        char tmptmp[256];
        sprintf(tmptmp, " -k %d ", ini_k);
        command += tmptmp;
        
        if (!bSkip) {
            cout << "Running: " << command << endl;
            Execute(command.c_str());
        }
    } else {
        cout << "Processing Inchworm outputs!" << endl;
        iniSeqs = iwormString;
        
        if (!Exists(iwormString)) {
            return -1;
        }
        
    }
    
    if (!bStrand)
        bBreak = false;
    
    
    if (!Exists(readString)) {
        return -1;
    }
    
    
  ////////////////////////////////////////////////////////////////////////////////////
  // Original break/join code that was later discarded as part of the Trinity process
  /*  
  // Break...
  command = exec + "BreakTransByPairs -i ";
  command += readString;
  command += " -f ";
  command += iniSeqs;
  command += " -o ";
  command += outDir + "/broken.out";
  if (bBreak) {
    cout << "Running: " << command << endl;
    Execute(command.c_str());
  }


  //...and Join
  // Break...
  command = exec + "JoinTransByPairs -i ";
  command += readString;
  command += " -f ";
  if (bBreak) 
    command += outDir + "/broken.out";
  else
    command += iniSeqs;
    
  command += " -o ";
  command += outDir + "/scaffolds.out";
  
  cout << "Running: " << command << endl;
  Execute(command.c_str());
  */
  /////////////////////////////////////////////////////////////////////////////////////




    /////////////////////////////////////////////////////////////////////////////////////
    // GraphFromFasta:
    // Pool sequences (clustering of Inchworm contigs, join using read-support)
    // Create output file: components.out
    
    string graphFromFastaCompleteFile = outDir + "/GraphFromFasta.finished";

    command = exec + "GraphFromFasta -i ";
    command += iniSeqs;
    //command += outDir + "/scaffolds.out > ";
    command += " -r ";
    command += readString;
    
    if (bStrand) {
        command += " -strand ";
    }
    command +=  " > ";
    string components_file = outDir + "/components.out";
    command += components_file;

    if (Exists(components_file) && Exists(graphFromFastaCompleteFile)) {
        cerr << "File: " << components_file << " already exists. Using existing file assuming resume mode" << endl << endl;
    }
    else {
        cout << "Running: " << command << endl;
        Execute(command.c_str());

        // make completion-marking file
        string complete_file_cmd = "touch " + graphFromFastaCompleteFile;
        Execute(complete_file_cmd.c_str());
        
    }
    //
    ///////////////////////////////////////////////////////////////////////////////////////
    
    

    ///////////////////////////////////////////////////////////////////////////////////////
    // Read components.out, create chrysalis/bundles.fasta which represents clustered inchworm sequences
    // Create de Bruijn graphs based on clustered Inchworm seqs (partitioned chrysalis/RawComps.\d+/comp\d+.raw.graph files) 

    string iString = outDir + "/components.out";
    //string oString = outDir + "/graph.out";
    
    string tmpName = outDir + "/tmp.fasta";
    //FILE * p = fopen(oString.c_str(), "w");
    //fclose(p);
    
    FlatFileParser parser;  
    parser.Open(iString);  // read components.out file
    int k = 0;
    
    ComponentFileMgr mgr(outDir);
    
    
    vecDNAVector tmpSeq;
    
    vecDNAVector bundled;
    //bundled.reserve(1000000);
    string bundledName = outDir + "/bundled.fasta";
    DNAVector separator;
    separator.SetFromBases("X");
    
    FILE * pOut = NULL;
    
    while (parser.ParseLine()) {
        if (parser.GetItemCount() == 0)
            continue;
        
        if (parser.AsString(0) == "COMPONENT") {
            int n = parser.AsInt(2);
            int c = parser.AsInt(1);
            //p = fopen(oString.c_str(), "a");
            //fprintf(pOut, "Component %d\n", c);
            //fclose(p);
            
            
            //p = fopen(tmpName.c_str(), "w");
            
            tmpSeq.resize(0);
            
            while (parser.ParseLine()) {
                
                if (parser.GetItemCount() == 0)
                    continue;
                
                if (parser.AsString(0) == "END") {
                    break;	  
                }
                
                const char * ppp = parser.AsString(0).c_str();
                
                if (ppp[0] == '>') {
                    DNAVector tmp;
                    tmpSeq.push_back(tmp);
                    continue;
                }
                
                DNAVector app;
                app.SetFromBases(parser.AsString(0));
                // cout << "adding: " << parser.AsString(0) << endl;
                tmpSeq[tmpSeq.isize()-1] += app;
                
                //fprintf(p, "%s\n", parser.Line().c_str());
            }
            //fclose(p);
            
            if (tmpSeq.isize() == 1 && tmpSeq[0].isize() < minLen) {
                // cerr << "-discarding entry, too short." << endl;
                continue;
            }
            
            bool bGood = false;
            for (int x=0; x<tmpSeq.isize(); x++) {
                if (tmpSeq[x].isize() > minDumpLen) {
                    bGood = true;
                    break;
                }
            }
            
            if (!bGood) // no inchworm contig > minDumpLen
                continue;
            
            DNAVector toAdd;
            for (int x = 0; x<tmpSeq.isize(); x++) {
                toAdd += tmpSeq[x];
                toAdd += separator;
            }
            
            bundled.push_back(toAdd);
            
            //cout << "Printing component: " << k << endl;
            
            // Create and write de Bruijn graph for this component
            string name = mgr.GetFileName(k, ".raw.graph");
            pOut = fopen(name.c_str(), "w");
            
            fprintf(pOut, "Component %d\n", k);
            
            TranscriptomeGraph(tmpSeq,
                               pOut,
                               24);
            
            tmpSeq.resize(0);
            
            fclose(pOut);
            
            
            //string cmd = "./TranscriptomeFromKmers -k 24 -append -o " + oString + " -i " + tmpName;
            //    system(cmd.c_str());
            k++;
            //if (k == 2)
            //break;
        }
        
    }
    
    if (tmpSeq.isize() > 0) {
        
        if (tmpSeq.isize() > 1 || tmpSeq[0].isize() >= minLen) {
            DNAVector toAdd;
        
            for (int x = 0; x<tmpSeq.isize(); x++) {
                toAdd += tmpSeq[x];
                toAdd += separator;
            }
      
            bundled.push_back(toAdd);
            
            string name = mgr.GetFileName(k, ".raw.graph");
            pOut = fopen(name.c_str(), "w");
            
            fprintf(pOut, "Component %d\n", k);
            
            TranscriptomeGraph(tmpSeq,
                               pOut,
                               24);
            fclose(pOut);
        }
    }
    //  fclose(pOut);
    
    bundled.Write(bundledName);
    //
    ////////////////////////////////////////////////////////////////////
    

    
    /////////////////////////////////////////////////////////////////////
    // Make fasta files for each component...
    string readcounts_file = outDir + "/readcounts.out";
    string readsToTranscriptsCompleteFile = outDir + "/readsToTranscripts.finished";
    
    if (! (Exists(readcounts_file) && Exists(readsToTranscriptsCompleteFile) ) ) {

            command = exec + "ReadsToTranscripts -i ";
            command += readString;
            command += " -f  ";
            command += bundledName;
            command += " -o ";
            command += outDir;
            if (bStrand) {
                command += " -strand";
            }
            //if (bStrand)
            //command += " -strand";
            cout << "Running: " << command << endl;
            Execute(command.c_str());
            
            string complete_cmd = "touch " + readsToTranscriptsCompleteFile;
            Execute(complete_cmd.c_str());
    }
    else {
        cerr << "File: " << readsToTranscriptsCompleteFile << " already exists. Using existing read/transcript mappings assuming resume mode." << endl << endl;
    }
    
    
    // readcounts_file created by the above.
    FlatFileParser readCount;  
    readCount.Open(readcounts_file);
    readCount.ParseLine();
    string numReads = readCount.Line();
        
    
    
    //////////////////////////////////////////////////////////////////////
    // QuantifyGraph
    //////////////////////////////////////////////////////////////////////

    //svec<string> targets;
    string butterflyCommands = outDir + "/butterfly_commands";
    FILE * pButterfly = fopen(butterflyCommands.c_str(), "w");

    string quantifyGraphCommands = outDir + "/quantifyGraph_commands";
    FILE * pQGraphCmds = fopen(quantifyGraphCommands.c_str(), "w");
    
    for (int i=0; i<k; i++) {
        string graph = mgr.GetFileName(i, ".raw.graph");
        string fasta = mgr.GetFileName(i, ".raw.fasta");
        string finalgraph = mgr.GetFileName(i, ".out");
        //string fnialreads = mgr.GetFileName(k, ".finalgraph.reads");
        
        
        FILE * pFasta = fopen(fasta.c_str(), "r");
        if (pFasta == NULL)
            continue;
       
        fclose(pFasta);
        
        //string target = "java -jar /seq/vsag_scratch01/manfred/data/Trinity/Butterfly/Butterfly.jar -A -N 22000000 -L 350 -F 300 -C ";
        string target = "java -jar " + butterflyExec + " -N " + numReads;
        char tmp[256];
        sprintf(tmp, " -L %d -F %d -C ", minLen, pairDist);
        
        target += tmp;
        
        
        target += mgr.GetFileName(i, "");
        fprintf(pButterfly, "%s\n", target.c_str());
        
        
        command = exec + "QuantifyGraph -g ";
        command += graph + " -o ";
        command += finalgraph + " -i ";
        command += fasta;
        if (bStrand)
            command += " -strand";
        
        //cout << "Running: " << command << endl;
        // Execute(command.c_str());
    
        fprintf(pQGraphCmds, "%s\n", command.c_str());
    }
    fclose(pButterfly);
    fclose(pQGraphCmds);

    
    //////////////////////////////////////////////
    // can run butterfly directly (typically not done here, but done as part of the Trinity.pl wrapper)
    
    // not any more...   both quantify graph and butterfly commands are now written to files, to be processed in parallel by cmd_process_forker.pl

    /*
    if (bButt) {
        command = exec + "RunButterfly -i ";
        command += butterflyCommands;
        char tmptmp[256];
        sprintf(tmptmp, "%d", nCPU);
        command += " -n ";
        command += tmptmp;
        cout << "Running: " << command << endl;
        system(command.c_str());
        
    }
    //////////////////////////////////////////////
    */
    
    return(0);

}
