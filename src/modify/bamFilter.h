#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include <htslib/sam.h>
#include <cstdlib>
#include <gzstream.h>
using namespace std;
typedef long long llong ;


class Para_FF {
	public:
		string input_Sam ;
		string OutPut ;
		llong start ;
		llong end ;
		string chr;
		int minLeng ;
		int minMapQ ;
		int maxHit ;
		bool TF ;

		Para_FF()
		{
			start=0;
			end=1000000000;
			minLeng=30 ;
			minMapQ =15 ;
			maxHit=1 ;
			TF=false;
			chr="";
		}
};




int  print_Ausage_FF()
{
	cout <<""
		"\n"
		"\tUsage: BamFilter -InPut <in.bam> -OutPut <out.bam>\n"
		"\n"
		"\t\t-InPut     <str>   InPut sam/bam File\n"
		"\t\t-OutPut    <str>   OutPut bam \n"
		"\n"
		"\t\t-MinMapQ   <int>   The min Bwa Mapping Q [15]\n" 
		"\t\t-MinLeng   <int>   the min length of read[30]\n" 
		"\t\t-Start     <int>   the start position out[0]\n" 
		"\t\t-End       <int>   the end position out[1e9]\n" 
		"\t\t-Chr       <str>   only the [chr] out[all]\n"
		"\n"
		"\t\t-help              show this help\n" 
		"\n";
	return 1;
}


int parse_Acmd_FF(int argc, char **argv, Para_FF * para_FF)
{
	if (argc <=2 ) {print_Ausage_FF();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FF->input_Sam=argv[i];
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_FF->OutPut=argv[i];
		}
		else if (flag  ==  "Chr")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_FF->chr=argv[i];
		}
		else if (flag  ==  "MinMapQ")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FF->minMapQ=atoi(argv[i]);
		}
		else if (flag  ==  "MinLeng")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_FF->minLeng=atoi(argv[i]);
		}
		///*////
		else if (flag  ==  "End")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FF->end=atol(argv[i]);
		}
		else if (flag  ==  "Start")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_FF->start=atol(argv[i]);
		}
		else if (flag  == "help")
		{
			print_Ausage_FF();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_FF->input_Sam).empty() ||   (para_FF->OutPut).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;    
	}

	return 1 ;
}



//programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int bam_Filter_main(int argc, char **argv)
{
	Para_FF * para_FF = new Para_FF;
	if (parse_Acmd_FF(argc, argv , para_FF )==0)
	{
		delete para_FF ; 
		return 0 ;
	}

	bam_hdr_t *header;
	bam1_t *aln = bam_init1();
	samFile *in = sam_open((para_FF->input_Sam).c_str(), "r");

	samFile *outR1 = sam_open(para_FF->OutPut.c_str(), "wb");

	header = sam_hdr_read(in);

	if (sam_hdr_write(outR1, header) < 0) 
	{
		fprintf(stderr, "Error writing output.\n");
		exit(-1);
	}
	int tmp=0;
	while (sam_read1(in, header, aln) >= 0)
	{
		if ( (aln->core).qual < (para_FF->minMapQ) ) 
		{
			continue ;
		}
		if ( (aln->core).pos < ((para_FF->start))  ||    (aln->core).pos > (para_FF->end) ) 
		{
			continue ;
		}
		if  ((aln->core).l_qseq   < ( para_FF->minLeng ))
		{
			continue ;
		}
		string chrID="*";
		if ((aln->core).tid >= 0)
		{ // chr
			chrID=header->target_name[(aln->core).tid];
		}

		if (!(para_FF->chr).empty()   &&   chrID!=(para_FF->chr))
		{
			continue ;
		}

		tmp=sam_write1(outR1, header, aln);
	}

	sam_close(in);
	sam_close(outR1);


	delete para_FF ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
