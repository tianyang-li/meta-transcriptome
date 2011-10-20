#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);

use FindBin;

use Cwd;

$ENV{LC_ALL} = 'C';

my $util_dir = "$FindBin::Bin/../util";

{
	my @progs_req = qw (bowtie-build bowtie);
	
	foreach my $prog (@progs_req) {
		unless (`which $prog` =~ /^\//) {
			die "Error, cannot find path to required program: $prog";
		}
	}
}

my $usage = <<_EOUSAGE_;

Note:  Uses bowtie to map reads.

#########################################################################################
#
# Required:
#
#  --iworm             inchworm assembled contigs
#
#  --left_fq           left fragment fq file
#  --right_fq          right fragment fq file
#
# Optional (if strand-specific RNA-Seq):
#  
#  --SS_lib_type       RF or FR
#
#  --work_dir         directory to perform data processing (default: workdir.\$pid
#
#
#  --JUST_ALIGN       only bowtie read alignments are done, pipeline ends after that step.
#
###########################################################################################

_EOUSAGE_

	;


my $inchworm_contigs;
my $left_fq;
my $right_fq;
my $SS_lib_type;
my $work_dir;
my $JUST_ALIGN_FLAG;

&GetOptions( 'iworm=s' => \$inchworm_contigs,
			 'left_fq=s' => \$left_fq,
			 'right_fq=s' => \$right_fq,
			 'SS_lib_type=s' => \$SS_lib_type,
			 'work_dir=s' => \$work_dir,
			 'JUST_ALIGN' => \$JUST_ALIGN_FLAG,
			 );


unless ($inchworm_contigs && $left_fq && $right_fq) {
	die $usage;
}


unless ($work_dir) {
	$work_dir = "jaccard_clip_workdir.$$";
}

main: {
	
	my $curr_dir = cwd();
	
	# create full paths to inputs, if not set.
	foreach my $file ($inchworm_contigs, $left_fq, $right_fq) {
		
		unless ($file =~ /^\//) {
			$file = "$curr_dir/$file";
		}
	}

	my $outdir = $work_dir;
	unless (-d $outdir) {
		mkdir ($outdir) or die "Error, cannot mkdir $outdir";
	}
	
	&process_cmd("ln -s $inchworm_contigs $outdir/iworm.fa") unless (-e "$outdir/iworm.fa");
	&process_cmd("ln -s $left_fq $outdir/left.fq") unless (-e "$outdir/left.fq");
	&process_cmd("ln -s $right_fq $outdir/right.fq") unless (-e "$outdir/right.fq");

	chdir $outdir or die "Error, cannot cd to $outdir";
	
	&process_cmd("bowtie-build iworm.fa iworm") unless (-e "iworm.1.ebwt");
	
	## align left reads:
	&process_cmd("bowtie -q -S --sam-nohead iworm left.fq > left.fq.sam") unless (-s "left.fq.sam");
	&process_cmd("sort -T . -S 2G left.fq.sam > left.fq.nameSorted.sam") unless (-s "left.fq.nameSorted.sam");

	## align right reads:
	&process_cmd("bowtie -q -S -sam-nohead iworm right.fq > right.fq.sam") unless (-s "right.fq.sam");
	&process_cmd("sort -T . -S 2G right.fq.sam > right.fq.nameSorted.sam") unless (-s "right.fq.nameSorted.sam");
	
	## combine left and right into single file.
	&process_cmd("$util_dir/merge_left_right_nameSorted_SAMs.pl --left_sam left.fq.nameSorted.sam --right_sam right.fq.nameSorted.sam > combined.nameSorted.sam") unless (-s "combined.nameSorted.sam");

	## sort by coordinate.
	&process_cmd("sort -T . -S 2G -k 3,3 -k 4,4n combined.nameSorted.sam > combined.coordSorted.sam") unless (-s "combined.coordSorted.sam");
	
	
	my $final_sam_file = "combined.coordSorted.sam";
	
	if ($SS_lib_type) {
		## separate by strand:
		&process_cmd("$util_dir/SAM_strand_separator.pl combined.coordSorted.sam $SS_lib_type");
	
		$final_sam_file = "combined.coordSorted.sam.+.sam";
	}

	if ($JUST_ALIGN_FLAG) {
		print "Just aligning reads.  No splitting.  Alignments are now complete.\n";
		exit(0);
	}
	

	## run Jaccard computation:
	&process_cmd("$util_dir/SAM_ordered_pair_jaccard.pl $final_sam_file > $final_sam_file.jaccard.wig");

	## define the transcript clip points:
	&process_cmd("$util_dir/jaccard_wig_clipper.pl $final_sam_file.jaccard.wig > $final_sam_file.jaccard.wig.clips");

	## clip the inchworm transcripts:
	&process_cmd("$util_dir/jaccard_fasta_clipper.pl $inchworm_contigs  $final_sam_file.jaccard.wig.clips > $inchworm_contigs.clipped.fa");
	
	
	exit(0);
	
}


####
sub process_cmd {
	my ($cmd) = @_;
	
	print "CMD: $cmd\n";

	my $ret = system($cmd);
	
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}

	return($ret);
}
