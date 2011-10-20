#include <stdlib.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <math.h>

#include "IRKE.hpp"
#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"
#include "argProcessor.hpp"


using namespace std;


int run_IRKE(int argc, char* argv[]);

void examine_CDS_paths(IRKE& irke, string cds_fasta_filename, unsigned int min_cov, float min_connectivity, float min_entropy,
					   bool write_coverage_flag, string coverage_filename);


void thread_sequences_through_graph(IRKE& irke, string cds_fasta_filename);

string usage(ArgProcessor args);

int main (int argc, char* argv[]) {
    
    try {
        return(run_IRKE(argc, argv));
    }
    
    catch (string err) {
        cerr << err << endl;
    }
    
    return(1);
    
}


int run_IRKE(int argc, char* argv[]) {
    
    string progname (argv[0]);
    
    // defaults:
    unsigned int MAX_RECURSION = 1;  // inchworm step size.
    
    
    unsigned int MONITOR = 0; // verbose... beware, becomes I/O bound.
    
    bool describe_kmers = false; // if true, just dumps a report of all kmers in the graph and then exits.
    
    string reads_fasta_filename;
    string kmers_fasta_filename;
    unsigned int kmer_length = 25;
    unsigned int min_kmer_count = 1;
    float MIN_CONNECTIVITY_RATIO = 0.0;
    float min_any_entropy = 0.0;
    
    bool double_stranded_flag = false;  // strand-specific by default
    
    /*************************************/
    // Report inchworm assemblies option:
    bool run_inchworm = false;
    unsigned int MIN_ASSEMBLY_LENGTH = 48; // minimum length of an inchworm assembly for reporting.
    unsigned int MIN_ASSEMBLY_COVERAGE = 2; // minimum average kmer coverage for assembly to be reported.
    float MIN_SEED_ENTROPY = 1.5; // minimum entropy for a Kmer to seed in inchworm assembly construction.
    unsigned int MIN_SEED_COVERAGE = 2; // minimum kmer coverage for a seed.
    
    bool PACMAN = false;  // disables revisiting nodes previously visited within an inchworm session
    bool CRAWL = false;
    unsigned int crawl_length = 1;
    
    bool reassembleIworm = false;
    unsigned long pruneSingletonKmersReadInterval = 0;
    
    /***********************************
    // Report the debruijn graph option:
    bool report_components = false;
    unsigned int MIN_COMPONENT_SIZE = 100; // min number of kmers in a debruijn isolated graph for reporting it.
    bool summarize_component_size_only_flag = false; // don't report kmers, just the size of the components.
    */
    
    /***************************************/
    // search for CDS in the graph options:
    bool search_CDS_flag = false;
    bool thread_CDS_flag = false;
    string cds_fasta_filename;
    
    
    
    // for iterating through minimum coverage statistics
    unsigned int max_min_coverage = 0; // will be adjusted below.
    unsigned int cov_step = 1; 
    
    float max_con_ratio = 0.0;
    float con_step = 0.05;
    
    float max_min_entropy = 0.0;
    float entropy_step = 0.1;
    
    bool WRITE_COVERAGE = false;
    string COVERAGE_OUTPUT_FILENAME;
    
    
    // error-containing kmer removal:
    bool prune_error_kmers = true;
    float min_ratio_non_error = 0.05f;
    
    long start, end;

    
    /***************************************/
    /***   Option Processing   ************/
    /**************************************/
    
    
    try {
        
        ArgProcessor args(argc, argv);
        
        /* Check for essential options */
        if (args.isArgSet("--help") 
            ||  (! (args.isArgSet("--reads") || args.isArgSet("--kmers")))
            
            || (! (args.isArgSet("--run_inchworm") 
//                   || args.isArgSet("--component_mode")
                   || args.isArgSet("--checkFastaPath")
                   || args.isArgSet("--threadFasta")
                   || args.isArgSet("--describe_kmers")
                   
                   ) 
                )
            ) {
            cerr  << usage(args) << endl << endl;
            
            return 1;
        }
        
        /* Process option values */
        
        if (args.isArgSet("--reads"))
        		reads_fasta_filename = args.getStringVal("--reads");
        
        if (args.isArgSet("--kmers"))
	        kmers_fasta_filename = args.getStringVal("--kmers");

        if (args.isArgSet("-K")) {
            kmer_length = args.getIntVal("-K"); 
            cerr << "Kmer length set to: " << kmer_length << endl;
        }
        
        if (args.isArgSet("--minKmerCount")) {
            
            min_kmer_count = args.getIntVal("--minKmerCount");
            cerr << "Min coverage set to: " << min_kmer_count << endl;
            
        }
        
        if (args.isArgSet("-L")) {
            MIN_ASSEMBLY_LENGTH = args.getIntVal("-L");
            cerr << "Min assembly length set to: " << MIN_ASSEMBLY_LENGTH << endl; 
        }
        if (args.isArgSet("--min_assembly_coverage")) {
            MIN_ASSEMBLY_COVERAGE = args.getIntVal("--min_assembly_coverage");
            cerr << "Min assembly coverage set to: " << MIN_ASSEMBLY_COVERAGE << endl;
        }
        
        if (args.isArgSet("--crawl")) {
            CRAWL = true;
            cerr << "Crawling enabled" << endl;
            if (args.isArgSet("--crawl_length")) {
                crawl_length = args.getIntVal("--crawl_length");
                cerr <<  "Crawl length set to: " << crawl_length << endl;
            }
        }
        
        if (args.isArgSet("--pacman")) {
            PACMAN = true;
            cerr << "PACMAN enabled." << endl;
        }
        
        if (args.isArgSet("--monitor")) {
            MONITOR = args.getIntVal("--monitor");
            cerr << "Monitor turned on, set to: " << MONITOR << endl;
        }
        
        if (args.isArgSet("--min_con_ratio")) {
            MIN_CONNECTIVITY_RATIO = args.getFloatVal("--min_con_ratio");
        }
        
        if (args.isArgSet("--DS")) {
            double_stranded_flag = true;
            cerr << "double stranded mode set" << endl;
        }
        
        if (args.isArgSet("--min_seed_entropy")) {
            MIN_SEED_ENTROPY = args.getFloatVal("--min_seed_entropy");
            cerr << "Min seed entropy set to: " << MIN_SEED_ENTROPY << endl;
        }
        if (args.isArgSet("--min_seed_coverage")) {
            MIN_SEED_COVERAGE = args.getIntVal("--min_seed_coverage");
            cerr << "min seed coverage set to: " << MIN_SEED_COVERAGE << endl;
        }
        
        if (args.isArgSet("--min_any_entropy")) {
            min_any_entropy = args.getFloatVal("--min_any_entropy");
            cerr << "min entropy set to: " << min_any_entropy << endl;
        }
        
        if (args.isArgSet("--run_inchworm")) {
            run_inchworm = true;
            
            if (args.isArgSet("-R")) {
                MAX_RECURSION = args.getIntVal("-R");
                cerr << "Setting max recursion to: " << MAX_RECURSION << endl;
            }
            
        }
/*
        if (args.isArgSet("--component_mode")) {
            report_components = true;
            
            if (args.isArgSet("--min_component_size")) {
                MIN_COMPONENT_SIZE = args.getIntVal("--min_component_size");
            }
            
            if (args.isArgSet("--component_size_only")) {
                summarize_component_size_only_flag = true;
            }
        }
*/
        
        if (args.isArgSet("--checkFastaPath")) {
            cds_fasta_filename = args.getStringVal("--checkFastaPath");
            search_CDS_flag = true;
            
            if (args.isArgSet("--max_min_cov")) {
                max_min_coverage = args.getIntVal("--max_min_cov");
                
                if (args.isArgSet("--cov_step")) {
                    cov_step = args.getIntVal("--cov_step");
                }
            }
            
            if (args.isArgSet("--max_con_ratio")) {
                max_con_ratio = args.getFloatVal("--max_con_ratio");
                
                if (args.isArgSet("--con_step")) {
                    con_step = args.getFloatVal("--con_step");
                }
            }
            
            if (args.isArgSet("--max_min_entropy")) {
                max_min_entropy = args.getFloatVal("--max_min_entropy");
                
                if (args.isArgSet("--entropy_step")) {
                    entropy_step = args.getFloatVal("--entropy_step");
                }
                
            }
        }
        
        
        if (args.isArgSet("--threadFasta")) {
            cds_fasta_filename = args.getStringVal("--threadFasta");
            thread_CDS_flag = true;
        }
        
        if (args.isArgSet("--describe_kmers")) {
            describe_kmers = true;
        }
        
        if (args.isArgSet("--coverage_outfile")) {
            WRITE_COVERAGE = true;
            COVERAGE_OUTPUT_FILENAME = args.getStringVal("--coverage_outfile");
        }
        
        if (args.isArgSet("--reassembleIworm")) {
            reassembleIworm = true;
            cerr << "Running in re-assembly mode." << endl;
        }
        
        if (args.isArgSet("--pruneSingletonKmersReadInterval")) {
            pruneSingletonKmersReadInterval = args.getLongVal("--pruneSingletonKmersReadInterval");
        }
        
        
        // kmer error removal options
        if (args.isArgSet("--no_prune_error_kmers")) {
            prune_error_kmers = false;
        }
        
        if (prune_error_kmers && args.isArgSet("--min_ratio_non_error")) {
            min_ratio_non_error = args.getFloatVal("--min_ratio_non_error");
            cerr << "Set to prune kmers below min ratio non-erro: " << min_ratio_non_error << endl;
        }
        
		
        /* adjust values of max if not set directly or set improperly */
        if (max_min_coverage < min_kmer_count) {
            max_min_coverage = min_kmer_count;
        }
        
        if (max_con_ratio < MIN_CONNECTIVITY_RATIO) {
            max_con_ratio = MIN_CONNECTIVITY_RATIO;
        }
        
        if (max_min_entropy < min_any_entropy) {
            max_min_entropy = min_any_entropy;
        }
        
        
    }
    
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    
    /**********************************************/
    /********  Run IRKE  **************************/
    /**********************************************/
    
    
    IRKE irke(kmer_length, MAX_RECURSION, MIN_SEED_ENTROPY, MIN_SEED_COVERAGE, min_any_entropy, PACMAN, CRAWL, crawl_length, double_stranded_flag);
    
    if (MONITOR > 0) {
        irke.monitor_progress(MONITOR); // helpful to see what's happening, but slows it down considerably.
    }
    
    if (pruneSingletonKmersReadInterval > 0) {
        irke.set_prune_singleton_read_interval(pruneSingletonKmersReadInterval);
    }
    
    if (reads_fasta_filename.length())
    		irke.build_graph(reads_fasta_filename, reassembleIworm, false);
    else if (kmers_fasta_filename.length())
    		irke.build_graph(kmers_fasta_filename, reassembleIworm, true);
    
/*
    if (min_kmer_count > 1) {
        cerr << endl << "-pruning kmers below coverage: " << min_kmer_count << " ...."; 
        start = time(NULL);
        irke.prune_kmers_min_count(min_kmer_count);
        end = time(NULL);
        cerr << "done, taking " << end-start << " seconds. " << irke.get_graph_size() << " kmers remain in graph." << endl << endl;
    }
    
    if (min_any_entropy > 0) {
        cerr << "-pruning kmers with entropy below: " << min_any_entropy << "....";
        start = time(NULL);
        irke.prune_kmers_min_entropy(min_any_entropy);
        end = time(NULL);
        cerr << "done, taking " << end-start << " seconds. " << irke.get_graph_size() << " kmers remain in graph." << endl << endl;
    }

    if (prune_error_kmers) {
        cerr << "Pruning kmers below min ratio of dominant kmer" << endl;
        start = time(NULL);
        irke.prune_kmer_extensions(min_ratio_non_error);
        end = time(NULL);
        cerr << "done, taking " << end-start << " seconds. " << irke.get_graph_size() << " kmers remain in graph." << endl << endl;
    }
*/
    if (min_kmer_count > 1 || min_any_entropy > 0 || prune_error_kmers) {
        cerr << "Pruning kmers (min_kmer_count=" << min_kmer_count
             << " min_any_entropy=" << min_any_entropy
             << " min_ratio_non_error=" << min_ratio_non_error << ")"<< endl;
        irke.prune_some_kmers(min_kmer_count, min_any_entropy, prune_error_kmers, min_ratio_non_error);
    }

    if (describe_kmers) {
        irke.describe_kmers();
    }
    
    else if (search_CDS_flag) {
        
        // iterate over min coverage values
        for (unsigned int minC = min_kmer_count; minC <= max_min_coverage; minC += cov_step) {
            
            // iterate over min connectivity values
            for (float min_conn = MIN_CONNECTIVITY_RATIO; min_conn <= max_con_ratio; min_conn += con_step) {
                
                for (float min_entropy = min_any_entropy; min_entropy <= max_min_entropy; min_entropy += entropy_step) {
                    
                    examine_CDS_paths(irke, cds_fasta_filename, minC, min_conn, min_entropy, WRITE_COVERAGE, COVERAGE_OUTPUT_FILENAME);
                    
                }
                
            }
        }
    }
    
    else if (thread_CDS_flag) {
        
        thread_sequences_through_graph(irke, cds_fasta_filename);
        
    }
/*
    else if (report_components) {
		
        cerr << "-reporting connected components ...";
        irke.report_connected_components(MIN_CONNECTIVITY_RATIO, MIN_COMPONENT_SIZE, summarize_component_size_only_flag); 
        cerr << " done." << endl;
    }
*/
    else if (run_inchworm) { // note this would have to be true since either or both are required.
        irke.compute_sequence_assemblies(MIN_CONNECTIVITY_RATIO, MIN_ASSEMBLY_LENGTH, MIN_ASSEMBLY_COVERAGE, WRITE_COVERAGE, COVERAGE_OUTPUT_FILENAME);
    }
    
    
    return (0);
    
}



void examine_CDS_paths(IRKE& irke, string cds_fasta_filename, unsigned int min_cov, float min_connectivity, float min_entropy,
					   bool WRITE_COVERAGE, string COVERAGE_OUTPUT_FILENAME) {
    
    
    Fasta_reader fasta_reader(cds_fasta_filename);
    
    ofstream coverage_writer;
    if (WRITE_COVERAGE) {
        coverage_writer.open(COVERAGE_OUTPUT_FILENAME.c_str());
        if (! coverage_writer.is_open()) {
            throw(stacktrace() + "Error, cannot write to file: " + COVERAGE_OUTPUT_FILENAME);
        }
    }
    
    
    while (fasta_reader.hasNext()) {
        
        Fasta_entry fe = fasta_reader.getNext();
        
        string accession = fe.get_accession();
        
        string seq = fe.get_sequence();
        
        vector<unsigned int> coverage_counter;
        
        cout << accession << "\tCov: " << min_cov << "\tCon: " << min_connectivity << "\tE: " << min_entropy;
        
        if (irke.sequence_path_exists(seq, min_cov, min_entropy, min_connectivity, coverage_counter)) {
            
            cout << "\tT" << endl;
            
        }
        else {
            
            cout << "\tF" << endl;
        }
        
        if (WRITE_COVERAGE) {
            
            coverage_writer << ">" << accession << endl;
            
            for (unsigned int i = 0; i < coverage_counter.size(); i++) {
                coverage_writer << coverage_counter[i];
                if ( (i +1) % 30 == 0) {
                    coverage_writer << endl;
                }
                else {
                    coverage_writer << " ";
                }
            }
            coverage_writer << endl;
        }
        
        
    }
    
    if (WRITE_COVERAGE) {
        coverage_writer.close();
    }
    
    
    return;
}



void thread_sequences_through_graph(IRKE& irke, string cds_fasta_filename) {
    
    Fasta_reader fasta_reader(cds_fasta_filename);
    
    while (fasta_reader.hasNext()) {
        
        Fasta_entry fe = fasta_reader.getNext();
        
        string accession = fe.get_accession();
        
        string seq = fe.get_sequence();
        
        cout << "// " << accession << endl;
        
        cout << irke.thread_sequence_through_graph(seq) << endl << endl;
        
    }
    
    return;
}

string usage (ArgProcessor args) {
    
    bool show_advanced = args.isArgSet("--show_advanced");
    
    stringstream usage_info;
    
    usage_info 
        << endl << endl
        << "Usage:" << endl
	<< "  inchworm --reads <filename.fasta> --run_inchworm [opts] " << endl
        << "  inchworm --kmers <filename.fasta> --run_inchworm [opts] " << endl << endl
        
        << "**Required" << endl
        << "--reads  <str>            " << ":fasta file containing reads" << endl
        << "--kmers  <str>            " << ":fasta file containing kmers" << endl
        << endl;
    
    
    
    usage_info
        << endl
        << "** Run modes:" << endl
        << "--run_inchworm           " << ":run inchworm, report sequences" << endl;
    
    if (show_advanced) {
        usage_info 
            << "    * advanced: " << endl
            << "    --describe_kmers         " << ":just describe kmers in the graph" << endl
            << "    --checkFastaPath <str>   " << ":check fastea seqs for path in graph" << endl
            << "    --threadFasta <str>      " << ":thread sequences through the graph, report kmer" << endl;
        
    }
    
    
    usage_info
        << endl
        << "** General opts" << endl
        << "-K <int>                         " << ":kmer length (default: 25, meaning 24mer overlaps)  max = 32 (stored as 64-bit integers, 2-base char encoding)" << endl
        << "--minKmerCount <int>             " << ":min kmer count, lower abundant kmers are pruned from graph (default: 1)" << endl
        << "-L <int>                         " << ":min contig length to be reported (default: 48)" << endl
        << "--min_assembly_coverage  <int>   " << ":min kmer coverage of an assembled contig to be reported (default: 2)" << endl
        << "--coverage_outfile <str>         " << ":file to store per-base kmer coverage for contigs" << endl
        << "--DS                             " << ":double-stranded RNA-Seq mode (not strand-specific)" << endl
        << "--no_prune_error_kmers           " << ":disable pruning of kmers that occur at below --min_ratio_non_error " << endl
        << "--min_ratio_non_error <float>    " << ":min ratio for low/high alternative extension that is not an error (default: 0.05" << endl
        
        ;
    
    if (show_advanced) {
        usage_info 
            << "    * advanced opts:" << endl
            << "    -R <int>                  " << ":maximum recursion length (default: 1)" << endl
            << "    --min_any_entropy <float> " << ":kmer must have at least this entropy to be added to the graph (default: 0)" << endl
            << "    --min_seed_entropy <float> " << ":min seed kmer entropy (default: 1.5)" << endl
            << "    --min_seed_coverage <int> " << ":min seed kmer coverage (default: 2)" << endl
            << "    --crawl                   " << ":drag tail of inchworm by --crawl_length" << endl
            << "    --crawl_length  <int>     " << ":length of tail end crawl (default: 1)" << endl
            << "    --pacman                  " << ":exclude kmers as they're visited in descending order of counts" << endl
            << "    --pruneSingletonKmersReadInterval <long> " << ":at num reads parsed interval, prune singleton kmers (default: 0)" << endl
            << "    --min_con_ratio <float> " << ":minimum connectivity ratio (min/max) between neighboring kmers (default: 0)" << endl
            << "    --reassembleIworm       " << ":redo assembly using previous inchworm outputs. Post-process to vary-K" << endl;
        
        
        
    }
    
    
    usage_info 
        << endl
        << "** Misc opts." << endl
        << "--help                   " << ":produce this menu." << endl
        << "--monitor <int>          " << ":verbosity. ( '1': recommended, '2': for debugging ) " << endl
        << "--show_advanced          " << ":more advanced options (mostly experimental)" << endl 
        << endl << endl << endl;
    
    
    return(usage_info.str());
	
}
