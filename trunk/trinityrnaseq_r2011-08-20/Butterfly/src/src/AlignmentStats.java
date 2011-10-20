import jaligner.Alignment;
import java.util.*;

public class AlignmentStats {

	int matches = 0;
	int mismatches = 0;
	int gaps = 0; // internal gaps
	int alignment_length = 0; 
	int max_internal_gap_length = 0;
	int left_gap_length = 0;
	int right_gap_length = 0;
	
	int total_not_matched = 0; // includes all gaps (internal and external) and mismatches
	
	// positions of the bounds of the alignments inside the terminal gaps
	int left_bound = -1;
	int right_bound = -1;

	HashMap<String,Integer> alignedBaseCounter;
	
	
	public AlignmentStats() {};
	
	public AlignmentStats increment_alignment_stats (AlignmentStats a) {
		
		AlignmentStats b = new AlignmentStats();
		
		b.matches = this.matches + a.matches;
		b.mismatches = this.mismatches + a.mismatches;
		b.gaps = this.gaps + a.gaps;
		b.alignment_length = this.alignment_length + a.alignment_length;
		b.max_internal_gap_length = Math.max(this.max_internal_gap_length, a.max_internal_gap_length);
		
		b.total_not_matched = this.total_not_matched + a.total_not_matched;
		
		return(b);
		
	}
	
	
	
	public AlignmentStats (Alignment alignment) {
	
		char[] sequence1 = alignment.getSequence1();
        char[] sequence2 = alignment.getSequence2();
        // char[] markup = alignment.getMarkupLine();
	
        String nameA = alignment.getName1();
        String nameB = alignment.getName2();
        
        alignedBaseCounter = new HashMap<String,Integer>();
        alignedBaseCounter.put(nameA, new Integer(0));
        alignedBaseCounter.put(nameB, new Integer(0));
        
        left_bound = compute_left_bound(sequence1, sequence2);
        right_bound = compute_right_bound(sequence1, sequence2);  
        
        alignment_length = right_bound - left_bound + 1;

//        boolean align_started = false; //Moran - was never used..

        for (int i = left_bound; i <= right_bound; i++) {

        	char c1 = sequence1[i];
        	char c2 = sequence2[i];


        	// check for gap
        	if (c1 == Alignment.GAP || c2 == Alignment.GAP) {
        		gaps++;
        	}
        	
        	
        	// examine match types
        	if (c1 != Alignment.GAP || c2 != Alignment.GAP) {        		
        		
        		if (c1 != Alignment.GAP) {
        			alignedBaseCounter.put(nameA, new Integer(alignedBaseCounter.get(nameA) + 1));
        		}
        		if (c2 != Alignment.GAP) {
        			alignedBaseCounter.put(nameB, new Integer (alignedBaseCounter.get(nameB) + 1));
        		}

        		if (c1 != Alignment.GAP && c2 != Alignment.GAP) {

        			if (c1 == c2) {	
        				matches++;
        			}
        			else {
        				mismatches++;
        			}
        		}

        	}

        }

        // check for internal gap lengths
        
        max_internal_gap_length = 0;
        
        right_gap_length = sequence1.length - right_bound - 1; // right end of sequence joins up again in shared vertex
        left_gap_length = left_bound;
        
        
        int gap_length_1 = 0;
        int gap_length_2 = 0;
        
        for (int i = left_bound; i <= right_bound; i++) {
        
        	// check sequence 1
        	char c1 = sequence1[i];
        	
        	if (c1 == Alignment.GAP) {
        		gap_length_1++;
        	}
        	else {
        		if (gap_length_1 > max_internal_gap_length)
        			max_internal_gap_length = gap_length_1;

        		gap_length_1 = 0;
        	}
        	
        	
        	// check sequence 2
        	char c2 = sequence2[i];
        	
        	if (c2 == Alignment.GAP) {
        		gap_length_2++;
        	}
        	else {
        		if (gap_length_2 > max_internal_gap_length)
        			max_internal_gap_length = gap_length_2;
        		gap_length_2 = 0;
        	}
        	
        }
    	// update total number of nt not matched
        total_not_matched = mismatches + this.gaps + this.left_gap_length + this.right_gap_length;
        
	}
     
	
	public int get_count_of_bases_in_aligned_region (String seqName) {
		return(alignedBaseCounter.get(seqName).intValue());
	}
	
	private int compute_left_bound (char[] sequence1, char[] sequence2) {
		
		for (int i = 0; i < sequence1.length; i++) {
			char a = sequence1[i];
			char b = sequence2[i];
			if (a != Alignment.GAP && b != Alignment.GAP) {
				return(i);
			}
		}
		
		throw new RuntimeException("no aligned bases found");
	}
	
	private int compute_right_bound (char[] sequence1, char[] sequence2) {
		
		for (int i = sequence1.length-1; i >= 0; i--) {
			char a = sequence1[i];
			char b = sequence2[i];
			if (a != Alignment.GAP && b != Alignment.GAP) {
				return(i);
			}
		}
		
		throw new RuntimeException("no aligned bases found");
	}
	
	
	public String toString() {
		
		String text = "matches: " + this.matches 
					+ ", mismatches: " + this.mismatches
					+ ", internal gaps: " + this.gaps
					+ ", left_gap: " + this.left_gap_length
					+ ", right_gap: " + this.right_gap_length
					+ ", total_diffs: " + this.total_not_matched;
		
		return(text);
		
	}
	
}
