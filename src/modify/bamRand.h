#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include <htslib/sam.h>
#include <cstdlib>
using namespace std;
typedef long long llong ;


int  print_Ausage_Xam04()
{
	cout <<""
		"\n"
		"\tUsage: bamRand -InPut <in.bam> -OutPut <out.bam>\n"
		"\n"
		"\t\t-InPut        <str>   InPut sam/bam File\n"
		"\t\t-OutPut       <str>   OutPut bam \n"
		"\n"
		"\t\t-proportion <float>   read random out proportion [0.1]\n" 
		"\t\t-seed         <int>   rand seed,default time\n"
		"\n"
		"\t\t-help                 show this help\n" 
		"\n";
	return 1;
}


int parse_Acmd_Xam04(int argc, char **argv, In3str1v * para_Xam04 )
{
	if (argc <=2 ) {print_Ausage_Xam04();return 0;}

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
			para_Xam04->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_Xam04->InStr2=argv[i];
		}
		else if (flag  ==  "seed")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_Xam04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "proportion")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam04->InF=atof(argv[i]);
		}
		else if (flag  == "help")
		{
			print_Ausage_Xam04();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_Xam04->InStr1).empty() ||   (para_Xam04->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;    
	}
	return 1 ;
}



//programme entry
///////// swimming in the sky and flying in the sea ////////////
int bamRand_main(int argc, char **argv)
{
	In3str1v * para_Xam04 = new In3str1v;
	if (parse_Acmd_Xam04(argc, argv , para_Xam04 )==0)
	{
		delete para_Xam04 ; 
		return 0 ;
	}

	bam_hdr_t *header;
	bam1_t *aln = bam_init1();
	samFile *in = sam_open((para_Xam04->InStr1).c_str(), "r");

	samFile *outR1 = sam_open(para_Xam04->InStr2.c_str(), "wb");

	header = sam_hdr_read(in);

	if (sam_hdr_write(outR1, header) < 0) 
	{
		fprintf(stderr, "Error writing output.\n");
		exit(-1);
	}

	if (para_Xam04->InInt==0)
	{
		srand((unsigned)time(NULL));
	}
	else
	{
		srand(para_Xam04->InInt);
	}
	 
	int Y=int(para_Xam04->InF*100);
	int X=0;

	while (sam_read1(in, header, aln) >= 0)
	{
		X=(rand()%100);
		if (X>Y)  { continue; }		
		X=sam_write1(outR1, header, aln);
	}
	sam_close(in);
	sam_close(outR1);


	delete para_Xam04 ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
