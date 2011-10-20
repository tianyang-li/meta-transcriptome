#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use File::Basename;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<_EOUSAGE_;

##################################################################
#
# Required:
#
#  --seqType <string>  :type of reads: (fq or fa)
#
#  If paired reads:
#
#      --left  <string>    :left reads
#      --right <string>    :right reads
# 
#  Or, if unpaired reads:
#
#      --single <string>   :single reads
#

#  --output <string>     :name of directory for output (will be created if it doesn't already exist) 
#                                 default( "trinity_out_dir" )
#
#  if strand-specific data, set:
#  
#      --SS_lib_type <string>  :if paired: RF or FR,  if single: F or R
#
#  
#
#  Butterfly-related options:
# 
#      --run_butterfly                 :executes butterfly commands.  Do not set this if you want to spawn them on a computing grid.   
#
#      --bfly_opts <string>            :parameters to pass through to butterfly (see butterfly documentation).
#
#      --bflyHeapSpace <string>        :java heap space setting for butterfly (default: 1000M) => yields command java -Xmx1000M -jar Butterfly.jar ... \$bfly_opts
#   
#  Inchworm-related options:
#
#      --no_meryl                      :do not use meryl for computing the k-mer catalog (default: uses meryl, providing improved runtime performance)
#      --min_kmer_cov <int>            :min count for K-mers to be assembled by Inchworm (default: 1)
#
#
# Misc:
#
#  --CPU <int>               :number of CPUs to use, default: 2 
#
#  --min_contig_length <int> :minimum assembled contig length to report (def=200)
#
#  --paired_fragment_length <int>  :maximum length expected between fragment pairs (aim for 90% percentile)  (def=300)
#
#  --jaccard_clip     :option, set if you have paired reads and you expect high gene density with UTR overlap (use FASTQ input file format for reads).
#
#  --run_ALLPATHSLG_error_correction :runs the read error correction process built into ALLPATHSLG.
#                                     (requires ALLPATHSLG to be installed, and installation directory indicated
#                                      by the env variable 'ALLPATHSLG_BASEDIR')
#
#
#####################################################################################################################################



_EOUSAGE_

	;



my $advanced_options = <<__ADVANCED__;

#
#  --inchworm <string>   :start with existing inchworm assemblies (skips the inchworm step given existing output)  
#  --USE_READS <string>  :use the fasta reads already extracted


__ADVANCED__
	
	;


my $ROOTDIR = "$FindBin::Bin";
my $UTILDIR = "$ROOTDIR/util";
my $INCHWORM_DIR = "$ROOTDIR/Inchworm";
my $CHRYSALIS_DIR = "$ROOTDIR/Chrysalis";
my $BUTTERFLY_DIR = "$ROOTDIR/Butterfly";
my $MERYL_DIR = "$ROOTDIR/trinity-plugins/kmer/meryl";

# option list:
my ($seqType, $left_file, $right_file, $single_file, $SS_lib_type, $min_contig_length,
	$paired_fragment_length, $jaccard_clip, $inchworm_file, $show_advanced_options,
	$output_directory, $run_ALLPATHSLG_error_correction_flag,
	);

# defaults:
$min_contig_length = 200;
$paired_fragment_length = 300;


my $RUN_BUTTERFLY_FLAG = 0;
my $bfly_opts = "";
my $bflyHeapSpace = "1G";

my $USE_MERYL = 1;
my $min_kmer_cov = 1;

my $CPU = 2;

## ADVANCED OPTIONS:
     # todo  add some.

my $use_reads;

my $no_meryl_flag = 0;



&GetOptions( 
             
             ## general opts
             "seqType=s" => \$seqType,
             "left=s" => \$left_file,
             "right=s" => \$right_file,
             "single=s" => \$single_file,
             
             "SS_lib_type=s" => \$SS_lib_type,
             
             "output=s" => \$output_directory,
             
             "min_contig_length=i" => \$min_contig_length,
             "paired_fragment_length=i" => \$paired_fragment_length,
             "jaccard_clip" => \$jaccard_clip,
             
             
             # Butterfly opts
             'run_butterfly' => \$RUN_BUTTERFLY_FLAG,
             'bfly_opts=s' => \$bfly_opts,
             'bflyHeapSpace=s' => \$bflyHeapSpace,
             
             'CPU=i' => \$CPU,
             'run_ALLPATHSLG_error_correction' => \$run_ALLPATHSLG_error_correction_flag,
             
             
             # Inchworm opts
             'no_meryl' => \$no_meryl_flag,
             'min_kmer_cov=i' => \$min_kmer_cov,
             
             "show_advanced_options" => \$show_advanced_options,
             
             
             ## Advanced options
             
             "inchworm=s" => \$inchworm_file,
             "USE_READS=s" => \$use_reads,
             
             );


if (@ARGV) {
    die "Error, do not understand options: @ARGV";
}


if ($no_meryl_flag) {
    $USE_MERYL = 0;
}

## Check options set:

unless ($seqType && 
		( 
		  ($left_file && $right_file) || $single_file || $use_reads) 
		) {
	die $usage;
}


## Check Java version:
my $java_version = `java -version 2>&1 `;
unless ($java_version =~ /java version \"1\.6\./) {
    die "Error, Trinity requires access to Java version 1.6.  Currently installed version is: $java_version";
}


if ($USE_MERYL && ! -e "$MERYL_DIR/meryl") {
    die "Error, cannot find meryl at $MERYL_DIR/meryl";
}


unless ($output_directory) {
	$output_directory = "trinity_out_dir";
}

if ($bflyHeapSpace !~ /^\d+[MG]$/) {
    die "Error, bflyHeapSpace must be set to a value of format: \\d+G or \\d+M  (eg. 1G or 1000M)";
}


if ($run_ALLPATHSLG_error_correction_flag) {
	unless ($ENV{ALLPATHSLG_BASEDIR}) {
		die "fatal: for error-correction by ALLPATHSLG, you must have ALLPATHSLG installed and the environmental"
			. " variable 'ALLPATHSLG_BASEDIR' set. ";
	}
	
	unless ($seqType eq "fq") {
		die "Error, the error-correction by AllpathsLG requires fastq files as a starting point. ";
	}

}


my $CPU_MAX = 10;
if ($CPU > $CPU_MAX) {
    ## comment this out if you really do want to use more processors.  I suspect it may only cause trouble, though.
    print STDERR "Warning, --CPU $CPU might be excessive.  Limiting it to $CPU_MAX for now.\n";
    $CPU = $CPU_MAX;
}


$ENV{OMP_NUM_THREADS} = $CPU; ## for Inchworm


main: {

    if ($RUN_BUTTERFLY_FLAG) {
        print STDERR "-since butterfly will eventually be run, lets test for proper execution of java\n";
        &test_java_failure_capture();
        print STDERR "-java tests succeeded.\n\n\n";
    }
    
    
	my $start_dir = cwd();

	## create complete paths for input files:
	$left_file = &create_full_path($left_file) if $left_file;
	$right_file = &create_full_path($right_file) if $right_file;
	$single_file = &create_full_path($single_file) if $single_file;
	$inchworm_file = &create_full_path($inchworm_file) if $inchworm_file;
	$output_directory = &create_full_path($output_directory);
	
	$use_reads = &create_full_path($use_reads) if $use_reads;

	unless (-d $output_directory) {
		
		mkdir $output_directory or die "Error, cannot mkdir $output_directory";
	}
	
	chdir ($output_directory) or die "Error, cannot cd to $output_directory";
	
	my $trinity_target_fa = "";
	
	if ($use_reads) {
		$trinity_target_fa = "reads.fa";
		&process_cmd("ln -s $use_reads $trinity_target_fa");
	}
	elsif ($left_file && $right_file) {
		
		if ($run_ALLPATHSLG_error_correction_flag) {
			&process_cmd("$ROOTDIR/util/run_AllpathsLG_error_correction.pl $left_file $right_file");
			$left_file = "$left_file.ErrCor.fq";
			$right_file = "$right_file.ErrCor.fq";
		}
		
		my ($left_SS_type, $right_SS_type);
		if ($SS_lib_type) {
			($left_SS_type, $right_SS_type) = split(//, $SS_lib_type);
		}
		
		&prep_seqs($left_file, $seqType, "left", $left_SS_type) unless (-s "left.fa");
		&prep_seqs($right_file, $seqType, "right", $right_SS_type) unless (-s "right.fa");
		
		$trinity_target_fa = "both.fa";
		&process_cmd("cat left.fa right.fa > $trinity_target_fa") unless (-s $trinity_target_fa);
		
	}
	elsif ($single_file) {
		if ($run_ALLPATHSLG_error_correction_flag) {
			&process_cmd("$ROOTDIR/util/run_AllpathsLG_error_correction.pl $single_file");
			$single_file = "$single_file.ErrCor.fq";
		}
		
		&prep_seqs($single_file, $seqType, "single", $SS_lib_type) unless (-s "single.fa");
		$trinity_target_fa = "single.fa";
	}
	
	else {
		die "not sure what to do. "; # should never get here.
	}
	
	#################
	## Inchworm step:
	
	unless ($inchworm_file) {
		$inchworm_file = &run_inchworm($trinity_target_fa, $SS_lib_type);
    }
	
	if ($jaccard_clip && $left_file && $right_file) {
		$inchworm_file = &run_jaccard_clip($inchworm_file, $left_file, $right_file, $seqType, $SS_lib_type);
	}
	
    unless (-s $inchworm_file) {
        die "Error, no Inchworm output is detected at: $inchworm_file";
    }
	
	
	##################
	## Chrysalis step:
	
	my $butterfly_cmds = &run_chrysalis($inchworm_file, $trinity_target_fa,
										$min_contig_length, $paired_fragment_length, $SS_lib_type);



	print "Inchworm and Chrysalis complete.  Butterfly commands to execute are provided here:\n"
		. $butterfly_cmds . "\n\n";
	


	if ($RUN_BUTTERFLY_FLAG) {

		my $cmd = "$ROOTDIR/util/cmd_process_forker.pl -c $butterfly_cmds --CPU $CPU --shuffle";  # shuffle them since the first ones are usually the longest-running ones.
		&process_cmd($cmd);

		
		## capture results:
		$cmd = 'find ./chrysalis -name "*allProbPaths.fasta" -exec cat {} \; > Trinity.fasta';

		&process_cmd($cmd);

        
		print "\n\n";
        print "###################################################################\n";
        print "Butterfly assemblies are written to $output_directory/Trinity.fasta\n";
        print "###################################################################\n\n\n";
		
        
	}
	
	exit(0);
	
}


####
sub run_chrysalis {
	my ($inchworm_file, $reads_file,
		$min_contig_length, $paired_fragment_length, $SS_lib_type) = @_;


    my $butterfly_cmds = &create_full_path("chrysalis/butterfly_commands");
    
    my $quantify_graph_cmds = &create_full_path("chrysalis/quantifyGraph_commands");
    
    my $adjusted_butterfly_cmds = "$butterfly_cmds.adj";
    
    if (-s $butterfly_cmds && -s $adjusted_butterfly_cmds) {
                
        print "#### WARNING: Chrysalis results already exist. Not rerunning Chrysalis. ####\n\n";
            
    }
    else {
        ## run Chrysalis
        
        my $cmd = "$CHRYSALIS_DIR/Chrysalis -i $reads_file -iworm $inchworm_file -o chrysalis "
            . " -min $min_contig_length -dist $paired_fragment_length ";
        
        if ($SS_lib_type) {
            $cmd .= " -strand 1 ";
        }
        
        $cmd .= " -butterfly $BUTTERFLY_DIR/Butterfly.jar ";
        
        eval {

            &process_cmd($cmd);
        };

        if ($@) {
            print "Error, the Chrysalis process failed:\n$@\n";
            print "In nearly all cases, this is related to not having the stacksize set to unlimited, a prerequisite to running Trinity.\n";
            print "Please visit:\n";
            print "\n       http://trinityrnaseq.sourceforge.net/trinity_faq.html#ques_E\n\n";
            print "for details.\n\n";
        }
    }
    
    
	unless (-s $butterfly_cmds) {
		croak "Error, chrysalis did not report butterfly commands file: $butterfly_cmds";
	}
	

    

    # see if we need to run the quantifyGraph commands:
    my $quantify_graph_cmds_finished = &create_full_path("chrysalis/quantifyGraph_commands.run.finished");
    if (! -e $quantify_graph_cmds_finished) {
        ## run it
        my $cmd = "$ROOTDIR/util/cmd_process_forker.pl -c $quantify_graph_cmds --CPU $CPU --shuffle";
        &process_cmd($cmd);
        
        &process_cmd("touch $quantify_graph_cmds_finished");
    }
    

    ## Rewrite the Butterfly commands

    ## add additional butterfly opts.
    open (my $fh, $butterfly_cmds) or die "Error, cannot read file $butterfly_cmds";
    open (my $ofh, ">$adjusted_butterfly_cmds") or die "Error, cannot write to $adjusted_butterfly_cmds";
    while (<$fh>) {
        my $line = $_;
        chomp $line;
        $line =~ s/^java /java -Xmx$bflyHeapSpace / or die "Error, cannot modify command";
        if ($bfly_opts) {
            $line .= " $bfly_opts";
        }
        print $ofh $line . "\n";
        #print STDERR $line . "\n";
    }
    close $ofh;
    close $fh;
    
        
	return($adjusted_butterfly_cmds);
	
}


####
sub run_inchworm {
	my ($reads, $strand_specific_flag) = @_;

    my $inchworm_outfile = "inchworm.K25.L48";

    my $inchworm_cmd;
    
    if ($USE_MERYL) {
        
        ## Strip header info out of fasta; meryl has limits as to how much header info it can store
        my $cmd = "$ROOTDIR/util/strip_fasta_header.pl $reads > $reads.headless";
        &process_cmd($cmd) unless (-s "$reads.headless");
        
        
        unless (-s "meryl_kmer_db.mcdat") {
            ## build k-mer db using meryl
            $cmd = "$MERYL_DIR/meryl -v -B -m 25 -s $reads.headless -o meryl_kmer_db";
            if ($strand_specific_flag) {
                $cmd .= " -f"; # forward strand k-mers only
            }
            else {
                $cmd .= " -C"; # canonical (one or other of the potential DS k-mer)
            }
            &process_cmd($cmd);
        }
        
        my $meryl_kmer_file = "meryl.kmers.min${min_kmer_cov}.fa";

        ## output k-mers
        $cmd = "$MERYL_DIR/meryl -Dt -n $min_kmer_cov -s meryl_kmer_db > $meryl_kmer_file";
        &process_cmd($cmd) unless (-s $meryl_kmer_file);
        
        
        $inchworm_cmd = "$INCHWORM_DIR/bin/inchworm --kmers $meryl_kmer_file --run_inchworm -K 25 -L 48 --monitor 1 ";
    }
    else {
        $inchworm_cmd = "$INCHWORM_DIR/bin/inchworm --reads $reads --run_inchworm -K 25 -L 48 --monitor 1 ";
        if ($min_kmer_cov > 1) {
            $inchworm_cmd .= " --min_kmer_coverage $min_kmer_cov ";
        }
    }
    
    
    unless ($strand_specific_flag) {
		$inchworm_cmd .= " --DS ";
		$inchworm_outfile .= ".DS";
	}
	
	$inchworm_outfile .= ".fa";
	
	$inchworm_cmd .= " 2>monitor.out > $inchworm_outfile";
	

    if (-s $inchworm_outfile) {
        print STDERR "####### WARNING:  Inchworm results: $inchworm_outfile already exist.  Reusing them rather than rerunning inchworm.#######\n";
    }
    else {
        
        eval {
            &process_cmd($inchworm_cmd);;
        };

        if ($@) {

            print STDERR "$@\n";
            print "** The inchworm process failed.  Below is the tail end of the log file:\n\n\n    $output_directory/monitor.out \n\n";
            system("tail -n10 monitor.out");
            print STDERR "\n\nIf it indicates bad_alloc(), then Inchworm ran out of memory.  You'll need to either reduce the size of your data set or run Trinity on a server with more memory available.\n\n";
            exit(1);
        }
    }
    

	return($inchworm_outfile);
	
}

####
sub prep_seqs {
	my ($initial_file, $seqType, $file_prefix, $SS_lib_type) = @_;

	if ($seqType eq "fq") {
		# make fasta
		
		my $cmd = "$UTILDIR/fastQ_to_fastA.pl -I $initial_file ";
		if ($SS_lib_type && $SS_lib_type eq "R") {
			$cmd .= " --rev ";
		}
		$cmd .= "> $file_prefix.fa";
		
		&process_cmd($cmd) unless (-e "$file_prefix.fa");
	}
	elsif ($seqType eq "fa") {
		if ($SS_lib_type && $SS_lib_type eq "R") {
			my $cmd = "$UTILDIR/revcomp_fasta.pl $initial_file > $file_prefix.fa";
			&process_cmd($cmd) unless (-s "$file_prefix.fa");
		}
		else {
			## just symlink it here:
			my $cmd = "ln -s $initial_file $file_prefix.fa";
			&process_cmd($cmd) unless (-s "$file_prefix.fa");
		}
	}
	
	return;
}



###
sub create_full_path {
	my ($file) = @_;

	my $cwd = cwd();
	if ($file !~ m|^/|) { # must be a relative path
		$file = $cwd . "/$file";
	}
	
	return($file);
}



####
sub process_cmd {
	my ($cmd) = @_;

	print "CMD: $cmd\n";

	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	
	return;
}


####
sub run_jaccard_clip {
	my ($inchworm_file, $left_file, $right_file, $seqType, $SS_lib_type) = @_;

	my $output_file = "$inchworm_file.clipped.fa";

    if (-s $output_file) {
        print STDERR "###### WARNING: $output_file already exists, skipping the jaccard-clip step, using already existing output: $output_file\n";
        return($output_file);
    }
    
	my $cmd = "$UTILDIR/inchworm_transcript_splitter.pl --iworm $inchworm_file "
		. " --left $left_file --right $right_file --seqType $seqType";

	if ($SS_lib_type) {
		$cmd .= " --SS_lib_type $SS_lib_type ";
	}
	
	&process_cmd($cmd);



	unless (-s $output_file) {
		croak "Error, jaccard clipping didn't produce the expected output file: $output_file";
	}

	return($output_file);
}


####
sub test_java_failure_capture {
    
    print "#######################################\n";
    print "Running Java Tests\n";
    
    my $cmd = "java -jar $UTILDIR/ExitTester.jar 0";
    eval {
        &process_cmd($cmd);
    };
    if ($@) {
        print STDERR "Error encountered in testing for running of a simple java application. ";
        print "$@\n\n";
        print STDERR "Please check your java configuration.\n";
        exit(1);
        
    }
    
    $cmd = "java -jar $UTILDIR/ExitTester.jar 1";
    eval {
        &process_cmd($cmd);
    };

    if ($@) {
        print "-we properly captured the java failure status, as needed.  Looking good.\n";
    }
    else {
        print STDERR "-we are unable to properly capture java failure status.  Please be sure that java (or any wrapper around java that's being used) can properly capture and propagate failure status before proceeding.\n";
        exit(1);
    }

    print "Java tests succeeded.\n";
    print "###################################\n\n";
    
    return;
}

