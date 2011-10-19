Gumbel parameters for sequence alignment statistics

The NCBI Gumbel parameters library computes Gumbel parameters for an
arbitrary scoring matrix and gap penalties. It also includes P-values
calculation using a new more accurate finite size correction method.
The parameters can be used for Basic Local Alignment Search Tool (BLAST).
Parameters for a gapless scoring scheme can be calculated as an option.

References
1. Park, Y., Sheetlin, S. and Spouge, J.L. (2005) Accelerated convergence
and robust asymptotic regression of the Gumbel scale parameter for gapped
sequence alignment. Journal of Physics A: mathematical and general, 38,
97-108.
2. Sheetlin, S., Park, Y. and Spouge, J.L. (2005) The Gumbel pre-factor k
for gapped local alignment can be estimated from simulations of global
alignment. Nucleic Acids Research, 33, 4987-4994.
3. Park, Y., Sheetlin, S. and Spouge, J.L. (2009) Estimating the Gumbel
Scale Parameter for Local Alignment of Random Sequences by Importance
Sampling with Stopping Times. Annals of Statistics.
4. Park, Y., Sheetlin, S. and Spouge, J.L. (2009) On-line computation of
statistical parameters for local alignment score distribution. In
preparation.

Additional information can be found here:
http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html.ncbi/blast/index.html
Description of the program.

To compile the C++ files in UNIX please replace the line

#define _MSDOS_

by the line

//#define _MSDOS_

in the file "sls_alp_data.hpp".

The program can be run in 2 different modes (separately or in a single run). 
1)Calculation of Gumbel parameters
2)Calculation of P-values using precalculated Gumbel parameters.

Mode 1. Calculation of Gumbel parameters.
USAGE:
alp13 <PARAMETERS>

8 mandatory PARAMETERS in the mode 1:
-scoremat <scoring matrix file name>
-freqs1 <background probabilities for the first sequence file name>
-freqs2 <background probabilities for the second sequence file name>
-gapopen <gap opening penalty>
-gapextend <gap extension penalty>
-K <relative accuracy for K calculation>
-lambda <relative accuracy for lambda calculation>
-time <time>

1 optional parameter in the mode 1:
-gumbelparout <output file name for calculated Gumbel parameters>

Background frequencies can be different for the sequence 1 and the sequence 2.

Example:

alp13 -scoremat test/blosum62 -freqs1 test/freqs.out -freqs2 test/freqs.out -gapopen 11 -gapextend 1 -K 0.005 -lambda 0.001 -max_time 10 -max_mem 2000 -gapped true -gumbelparout g.out

It will calculate Gumbel parameters for the matrix from the file test/blosum62; 
background probabilities from the file test/freqs.out; penalties 11/1; 
relative accuracy for K 0.5%; relative accuracy for lambda 0.1%; 
calculation time is 10 seconds; maximum allowed memory use is 2000 Mbs; 
the parameters will be outputted into the file named "g.out". 

Example of calculation of parameters for gapless scoring system:

alp13 -scoremat test/blosum62 -freqs1 test/freqs.out -freqs2 test/freqs.out -gapopen 11 -gapextend 1 -K 0.005 -lambda 0.001 -max_time 10 -max_mem 2000 -gapped false -gumbelparout g.out

Gap penalties will be ignored in this case.

Mode 2. Calculation of P-values using precalcualated Gumbel parameters.
USAGE:
alp13 <PARAMETERS>

4 mandatory PARAMETERS in the mode 2:
-score1 <the first score for P-values calculation within the range [score1,score2]>
-score2 <the second score for P-values calculation within the range [score1,score2]>
-seqlen1 <Length for the sequence 1>
-seqlen2 <Length for the sequence 2>
If the program does not have parameters sufficient for the mode 1 then 
the following fifth parameter is required:
-gumbelparin <Name of the file with precalculated Gumbel parameters>
One more optional parameter in the mode 2:
-pvalout <Output file name for calculated P-values>

The Gumbel parameters can be inputted from a file in the same format as 
generated in the mode 1 or can be calculated in the same run but in this case 
all parameters required for the mode 1 has to be set. 
For example the command:

alp13 -scoremat test/blosum62 -freqs1 test/freqs.out -freqs2 test/freqs.out -gapopen 11 -gapextend 1 -K 0.005 -lambda 0.001 -max_time 10 -max_mem 2000 -gapped true -seqlen1 200 -seqlen2 300 -score1 245 -score2 255

calculates Gumbel parameters at first and then uses them for calculation of 
P-values for the scores range [245,255] and sequences lengths 200 and 300.

Another example of usage is 

alp13 -gumbelparin g.out -seqlen1 200 -seqlen2 300 -score1 245 -score2 255 -pvalout p.out

where the program inputs all necessary parameters from the file "g.out" 
which has to be precalculated beforehand using the mode 1. 
The P-values will be written into the file "p.out".

The program uses parameters whenever it is possible and ignores otherwise. 
For example

alp13 -gumbelparin g.out -seqlen1 200 -seqlen2 300 -score1 245 -score2 255 -pvalout p.out -scoremat test/blosum62 -freqs1 test/freqs.out

does the same as 

alp13 -gumbelparin g.out -seqlen1 200 -seqlen2 300 -score1 245 -score2 255 -pvalout p.out

since not all required parameters are set for the mode 1.


