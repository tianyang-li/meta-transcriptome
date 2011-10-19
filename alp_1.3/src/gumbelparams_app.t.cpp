//File Description: Demo application for real time repeats Gumbel parameters computation


#include <vector>
#include <cfloat>

#include "sls_alp_data.hpp"
#include "sls_alp_sim.hpp"


#include "gumbel_params.hpp"

using namespace Sls;
using namespace std;


/////////////////////////////////////////////////////////////////////////////
void error_message()
{
	cout<<"The program can be run in two modes (separately or simultaneously):\nmode 1: calculation of Gumbel parameters for pairwise alignment\nmode 2: calculation of P-values\n\n";
	cout<<"USAGE:\nalp13 <PARAMETERS>\n\n8 mandatory PARAMETERS in the mode 1:\n-scoremat <scoring matrix file name>\n-freqs1 <background probabilities for the first sequence file name>\n-freqs2 <background probabilities for the second sequence file name>\n-gapopen <gap opening penalty>\n-gapextend <gap extension penalty>\n-K <relative accuracy for K calculation>\n-lambda <relative accuracy for lambda calculation>\n-time <time>\n1 optional parameter in the mode 1:\n-gumbelparout <output file name for calculated Gumbel parameters>\n\n";
	cout<<"USAGE:\nalp13 <PARAMETERS>\n\n4 mandatory PARAMETERS in the mode 2:\n-score1 <the first score for P-values calculation\nwithin the range [score1,score2]>\n-score2 <the second score for P-values calculation\nwithin the range [score1,score2]>\n-seqlen1 <Length for the sequence 1>\n-seqlen2 <Length for the sequence 2>\nif the program does not have parameters sufficient\nfor the mode 1 then the following fifth parameter is required:\n-gumbelparin <Name of the file with precalculated Gumbel parameters>\nOne more optional parameter in the mode 2:\n-pvalout <Output file name for calculated P-values>\n\n";

	throw error("Parameters of the program are wrong\n",4);
};

void command_line_interpreter(//extracts parameters from the command line
	int argc, char* argv[],//arguments of the program
	long int &rand_,//randomization number
	long int &gapopen_,//gap opening penalty
	long int &gapextend_,//gap extension penalty
	string &freqs1_file_name_,//probabilities file name
	string &freqs2_file_name_,//probabilities file name
	string &scoremat_file_name_,//scoring matrix file name
	double &eps_lambda_,//relative error for lambda calculation
	double &eps_K_,//relative error for K calculation
	double &max_time_,//maximum allowed calculation time in seconds
	double &max_mem_,//maximum allowed memory usage in MB
	string &gumbelparout_file_name_,//Gumbel parameters output file name
	bool &gapped_,//if true then performe estimation using gap penalties; 
			//if falce then gapless parameters will be estimated and
			//gap penalties will be ignored

	string &gumbelparin_file_name_,//Gumbel parameters input file name

	long int &seqlen1_,//length of the sequence 1
	long int &seqlen2_,//length of the sequence 2

	long int &score1_,
	long int &score2_,//P-values are calculated in the range [score1_,score2_]

	string &pvalout_file_name_,//P-values file name

	bool &Gumbel_mode_,//true if Gumbel parameters will be calculated
	bool &Pvalues_mode_)//true if P-values will be calculated
{

	if(argc%2!=1)
	{
		error_message();
	};

	Gumbel_mode_=false;
	Pvalues_mode_=false;

	rand_=-1;
	bool rand_flag=false;

	freqs1_file_name_="";
	bool freqs1_file_name_flag=false;

	freqs2_file_name_="";
	bool freqs2_file_name_flag=false;


	scoremat_file_name_="";
	bool scoremat_file_name_flag=false;

	gapopen_=-1;
	bool gapopen_flag=false;
	gapextend_=-1;
	bool gapextend_flag=false;

	eps_lambda_=DBL_MAX1;
	bool eps_lambda_flag=false;

	eps_K_=DBL_MAX1;
	bool eps_K_flag=false;

	max_time_=1;
	bool max_time_flag=false;

	max_mem_=500;
	bool max_mem_flag=false;

	gumbelparout_file_name_="";
	bool gumbelparout_file_name_flag=false;

	gapped_=true;
	bool gapped_flag=false;


//parameters for P-values calculation

	gumbelparin_file_name_="";
	bool gumbelparin_file_name_flag=false;

	seqlen1_=-1;
	bool seqlen1_flag=false;

	seqlen2_=-1;
	bool seqlen2_flag=false;

	score1_=-1;
	bool score1_flag=false;

	score2_=-1;
	bool score2_flag=false;

	pvalout_file_name_="pval.out";
	bool pvalout_file_name_flag=false;



	bool parameters_are_wrong_flag=false;
	long int i;
	for(i=0;i<floor(argc/2.0);i++)
	{

		string tmp=argv[2*i+1];

		if(tmp=="-gapped")
		{
			if(gapped_flag)
			{
				throw error("Error - parameter -gapped has been defined already\n",4);
			};
			string gapped_string=argv[2*i+2];
			if(gapped_string=="true")
			{
				gapped_=true;
			}
			else
			{
				if(gapped_string=="false")
				{
					gapped_=false;
				}
				else
				{
					throw error("Error - parameter -gapped can only take the values true or false\n",4);
				};
			};
			gapped_flag=true;
			continue;
		};

		if(tmp=="-freqs1")
		{
			if(freqs1_file_name_flag)
			{
				throw error("Error - parameter -freqs1 has been defined already\n",4);
			};
			freqs1_file_name_=argv[2*i+2];
			freqs1_file_name_flag=true;
			continue;
		};


		if(tmp=="-freqs2")
		{
			if(freqs2_file_name_flag)
			{
				throw error("Error - parameter -freqs2 has been defined already\n",4);
			};
			freqs2_file_name_=argv[2*i+2];
			freqs2_file_name_flag=true;
			continue;
		};


		if(tmp=="-scoremat")
		{
			if(scoremat_file_name_flag)
			{
				throw error("Error - parameter -scoremat has been defined already\n",4);
			};
			scoremat_file_name_=argv[2*i+2];
			scoremat_file_name_flag=true;
			continue;
		};

		if(tmp=="-rand")
		{
			if(rand_flag)
			{
				throw error("Error - parameter -rand has been defined already\n",1);
			};

			rand_=atol(argv[2*i+2]);
			rand_flag=true;
			continue;
		};


		if(tmp=="-gapopen")
		{
			if(gapopen_flag)
			{
				throw error("Error - parameter -gapopen has been defined already\n",1);
			};

			gapopen_=atol(argv[2*i+2]);
			gapopen_flag=true;
			continue;
		};

		if(tmp=="-gapextend")
		{
			if(gapextend_flag)
			{
				throw error("Error - parameter -gapextend has been defined already\n",1);
			};

			gapextend_=atol(argv[2*i+2]);
			gapextend_flag=true;
			continue;
		};

		//if(tmp=="-eps_lambda")
		if(tmp=="-lambda")
		{
			if(eps_lambda_flag)
			{
				throw error("Error - parameter -eps_lambda has been defined already\n",1);
			};

			eps_lambda_=atof(argv[2*i+2]);
			eps_lambda_flag=true;
			continue;
		};

		//if(tmp=="-eps_K")
		if(tmp=="-K")
		{
			if(eps_K_flag)
			{
				throw error("Error - parameter -eps_K has been defined already\n",1);
			};

			eps_K_=atof(argv[2*i+2]);
			eps_K_flag=true;
			continue;
		};

		if(tmp=="-max_time")
		{
			if(max_time_flag)
			{
				throw error("Error - parameter -max_time has been defined already\n",1);
			};

			// max_time_=atof(argv[2*i+2]);
			max_time_=DBL_MAX;
			max_time_flag=true;
			continue;
		};

		if(tmp=="-max_mem")
		{
			if(max_mem_flag)
			{
				throw error("Error - parameter -max_mem has been defined already\n",1);
			};

			max_mem_=atof(argv[2*i+2]);
			max_mem_flag=true;
			continue;
		};

		if(tmp=="-gumbelparin")
		{
			if(gumbelparin_file_name_flag)
			{
				throw error("Error - parameter -gumbelparin has been defined already\n",4);
			};
			gumbelparin_file_name_=argv[2*i+2];
			gumbelparin_file_name_flag=true;
			continue;
		};


		if(tmp=="-gumbelparout")
		{
			if(gumbelparout_file_name_flag)
			{
				throw error("Error - parameter -gumbelparout has been defined already\n",4);
			};
			gumbelparout_file_name_=argv[2*i+2];
			gumbelparout_file_name_flag=true;
			continue;
		};

		if(tmp=="-seqlen1")
		{
			if(seqlen1_flag)
			{
				throw error("Error - parameter -seqlen1 has been defined already\n",1);
			};

			seqlen1_=atol(argv[2*i+2]);
			seqlen1_flag=true;
			continue;
		};

		if(tmp=="-seqlen2")
		{
			if(seqlen2_flag)
			{
				throw error("Error - parameter -seqlen2 has been defined already\n",1);
			};

			seqlen2_=atol(argv[2*i+2]);
			seqlen2_flag=true;
			continue;
		};

		if(tmp=="-score1")
		{
			if(score1_flag)
			{
				throw error("Error - parameter -score1 has been defined already\n",1);
			};

			score1_=atol(argv[2*i+2]);
			score1_flag=true;
			continue;
		};

		if(tmp=="-score2")
		{
			if(score2_flag)
			{
				throw error("Error - parameter -score2 has been defined already\n",1);
			};

			score2_=atol(argv[2*i+2]);
			score2_flag=true;
			continue;
		};

		if(tmp=="-pvalout")
		{
			if(pvalout_file_name_flag)
			{
				throw error("Error - parameter -pvalout has been defined already\n",4);
			};
			pvalout_file_name_=argv[2*i+2];
			pvalout_file_name_flag=true;
			continue;
		};

	};

	Gumbel_mode_=(freqs1_file_name_flag&&freqs2_file_name_flag&&scoremat_file_name_flag&&gapopen_flag&&gapextend_flag&&eps_lambda_flag&&eps_K_flag);
	Pvalues_mode_=(seqlen1_flag&&seqlen2_flag&&score1_flag&&score2_flag);

	if(parameters_are_wrong_flag||!(Gumbel_mode_||Pvalues_mode_))
	{
		error_message();
	};

	if(Pvalues_mode_&&!Gumbel_mode_&&!gumbelparin_file_name_flag)
	{
		cout<<"Error - the parameter -gumbelparin is required for P-values calculation if the program does not have sufficient parameters for Gumbel calculation\n";
		error_message();
	};

	if(Pvalues_mode_&&Gumbel_mode_&&gumbelparin_file_name_flag)
	{
		cout<<"Warning - parameter -gumbelparin will be ignored\nsince the program computes Gumbel parameters\n\n";
	};

	if(Pvalues_mode_)
	{
		if(seqlen1_<=0||seqlen2_<=0)
		{
			cout<<"Error - please check the parameters -seqlen1 and -seqlen2\n";
		};

		if(score1_<0||score2_<0||score1_>score2_)
		{
			cout<<"Error - please check the parameters -score1 and -score2\n";
		};

	};

};

int main(int argc, char* argv[])
{
	long int max_number_of_unsuccessful_runs=1;

	long int init=max_number_of_unsuccessful_runs;//shows previous number of unsuccessful runs of the program
	long int rand;//randomization number
	long int gapopen;//gap opening 
	long int gapextend;//gap extension penalty
	string freqs1_file_name;//probabilities file name
	string freqs2_file_name;//probabilities file name
	string scoremat_file_name;//scoring matrix file name
	double eps_lambda;//relative error for lambda calculation
	double eps_K;//relative error for K calculation
	double max_time;//maximum allowed calculation time in seconds
	double max_mem;//maximum allowed memory usage in MB
	string gumbelparout_file_name;//probabilities file name
	bool gapped=true;

	string gumbelparin_file_name;//Gumbel parameters input file name

	long int seqlen1;//length of the sequence 1
	long int seqlen2;//length of the sequence 2

	long int score1;
	long int score2;//P-values are calculated in the range [score1_,score2_]

	string pvalout_file_name;//P-values file name

	bool Gumbel_mode=false;//true if Gumbel parameters will be calculated
	bool Pvalues_mode=false;//true if P-values will be calculated



	bool error_flag=false;

	try
	{
	try
	{


		command_line_interpreter(//extracts parameters from the command line
			argc, argv,//arguments of the program
			rand,//randomization number
			gapopen,//gap opening penalty
			gapextend,//gap extension penalty
			freqs1_file_name,//probabilities file name
			freqs2_file_name,//probabilities file name
			scoremat_file_name,//scoring matrix file name
			eps_lambda,//relative error for lambda calculation
			eps_K,//relative error for K calculation
			max_time,//maximum allowed calculation time in seconds
			max_mem,//maximum allowed memory usage in MB
			gumbelparout_file_name,//probabilities file name
			gapped,//if true then performe estimation using gap penalties; 
					//if falce then gapless parameters will be estimated and
					//gap penalties will be ignored
			gumbelparin_file_name,//Gumbel parameters input file name
			seqlen1,//length of the sequence 1
			seqlen2,//length of the sequence 2
			score1,
			score2,//P-values are calculated in the range [score1_,score2_]
			pvalout_file_name,//P-values file name
			Gumbel_mode,//true if Gumbel parameters will be calculated
			Pvalues_mode);//true if P-values will be calculated




		CGumbelParamsCalc::Run2(
			rand,//randomization number
			gapopen,//gap opening penalty
			gapextend,//gap extension penalty
			scoremat_file_name,//scoring matrix file name
			freqs1_file_name,//probabilities1 file name
			freqs2_file_name,//probabilities1 file name
			max_time,//maximum allowed calculation time in seconds
			max_mem,//maximum allowed memory usage in MB
			eps_lambda,//relative error for lambda calculation
			eps_K,//relative error for K calculation
			gumbelparout_file_name,
			gapped,

			gumbelparin_file_name,//Gumbel parameters input file name
			seqlen1,//length of the sequence 1
			seqlen2,//length of the sequence 2
			score1,
			score2,//P-values are calculated in the range [score1_,score2_]
			pvalout_file_name,//P-values file name
			Gumbel_mode,//true if Gumbel parameters will be calculated
			Pvalues_mode);//true if P-values will be calculated



		return 0;	

 	}
	catch (error er)
	{ 
		error_flag=true;


		if(er.error_code!=-1)
		{
			std::cout<<er.st;

		}
		else
		{
			if(er.error_code==2)
			{
				cout<<"Internal error in the program\n";
				//cout<<er.st<<endl;
			}
			else
			{
				cout<<"The previous attempt to estimate the parameters failed\n";
			};

			if(init+1<max_number_of_unsuccessful_runs)
			{
				cout<<"The program will be restarted\n";
			};

			throw error(er.st,er.error_code);
		};

		return 0;
	};
	}
	catch (...)
	{

		if(init+1>=max_number_of_unsuccessful_runs)
		{
			if(!error_flag)
			{
				cout<<"Internal error in the program\n";
			};
			std::cout<<"\nThe program cannot estimate the parameters\n";
			return 0;
		}
		else
		{

			string tmp=argv[0];
			long int i;
			for(i=1;i<argc;i++)
			{
				string tmp2=argv[i];
				if(tmp2=="-init")
				{
					break;
				};

				tmp+=" "+tmp2;
			};

			tmp+=" -init "+alp_data::long_to_string(init+1);

			cout<<tmp<<endl;

			system(tmp.data());

		};

		return 0;
	};



	return 0;
};

