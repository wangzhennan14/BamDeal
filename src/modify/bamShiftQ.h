#ifndef BamShiftQ_H_
#define BamShiftQ_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include "../ALL/dealfun.cpp"
#include <htslib/sam.h>
#include <cstdlib>
#include <sys/select.h>

using namespace std;
typedef long long llong ;



int  print_Ausage_shiftQ()
{
	cout <<""
		"\n"
		"\tUsage: BamShiftQ -InPut <in.bam> -OutPut <out.bam>\n"
		"\n"
		"\t\t-InPut     <str>   InPut sam/bam File\n"
		"\t\t-OutPut    <str>   OutPut bam\n"
		"\n"
		"\t\t-PhredQ    <int>   Final phred quality [1]\n"
		"\t\t                   [1:ASCII+33 2:ASCII+64]\n"
		"\t\t-MinMapQ   <int>   filter low min Bwa Mapping Q [10]\n" 
		"\t\t-MinLeng   <int>   filter low min length of read[30]\n" 
		"\n"
		"\t\t-help              show this help\n" 
		"\n";
	//"\t\t-ShiftQ    <int>   Final phred quality [+31/0/-31]\n" 
	return 1;
}

int parse_Acmd_Xam04(int argc, char **argv, Para_FF * para_Xam04)
{
	if (argc <=2 ) {print_Ausage_shiftQ();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","") ;

		if (flag  == "InPut" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam04->input_Sam=argv[i];
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_Xam04->OutPut=argv[i];
		}
		else if (flag  ==  "MinMapQ")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam04->minMapQ=atoi(argv[i]);
		}
		else if (flag  ==  "MinLeng")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_Xam04->minLeng=atoi(argv[i]);
		}
		else if (flag  ==  "PhredQ")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_Xam04->maxHit=atoi(argv[i]);
		}
		///*////
		else if (flag  == "help")
		{
			print_Ausage_shiftQ();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_Xam04->input_Sam).empty() ||   (para_Xam04->OutPut).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;    
	}

	string path=(para_Xam04->OutPut);

	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	if (ext != "bam")
	{
		(para_Xam04->OutPut)=path+".bam";
	}

	return 1 ;


}


//programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int bam_ShiftQ_main(int argc, char **argv)
{
	Para_FF * para_Xam04 = new Para_FF;
	if (parse_Acmd_Xam04(argc, argv , para_Xam04 )==0)
	{
		delete para_Xam04 ; 
		return 0 ;
	}

	samFile *outR1 = sam_open(para_Xam04->OutPut.c_str(), "wb");



	int ASCII_raw=Get_qual_Data ( (para_Xam04->input_Sam) );



	bam_hdr_t *header_SS;
	bam1_t *aln_SS = bam_init1();
	samFile *in_SS = sam_open((para_Xam04->input_Sam).c_str(), "r");


	header_SS = sam_hdr_read(in_SS);

	if (sam_hdr_write(outR1, header_SS) < 0) 
	{
		fprintf(stderr, "Error writing output.\n");
		exit(-1);
	}
	int tmp=1;
	if  (( para_Xam04->maxHit)==1   &&  (ASCII_raw!=0)   )
	{
		cout <<"PhredQ : ASCII+64   --->  ASCII+33 "<<endl;
		while  ( (sam_read1(in_SS, header_SS, aln_SS) >= 0) )
		{
			if ( (aln_SS->core).qual < (para_Xam04->minMapQ) ) 
			{
				continue ;
			}
			if  ((aln_SS->core).l_qseq   < ( para_Xam04->minLeng ))
			{
				continue ;
			}
			uint8_t  * seqQ=bam_get_qual(aln_SS);
			for(int i=0; i < (aln_SS->core).l_qseq ; i++)
			{
				seqQ[i]-=ASCII_raw ;
			}	
			tmp=bam_write1(outR1->fp.bgzf, aln_SS);
		}
	}
	else if ( ( para_Xam04->maxHit)!=1   &&  (ASCII_raw==0) )  
	{
		cout <<"PhredQ : ASCII+33   --->  ASCII+64 "<<endl;
		while  ( (sam_read1(in_SS, header_SS, aln_SS) >= 0) )
		{
			if ( (aln_SS->core).qual < (para_Xam04->minMapQ) ) 
			{
				continue ;
			}
			if  ((aln_SS->core).l_qseq   < ( para_Xam04->minLeng ))
			{
				continue ;
			}
			uint8_t  * seqQ=bam_get_qual(aln_SS);
			for(int i=0; i < (aln_SS->core).l_qseq ; i++)
			{
				seqQ[i]+=31 ;
			}	
			tmp=bam_write1(outR1->fp.bgzf, aln_SS);
		}
	}
	else
		//	if ( ( para_Xam04->maxHit)==1   &&  (ASCII_raw==0) )
	{
		if  (( para_Xam04->maxHit)==1 )
		{

			cout <<"PhredQ : ASCII+33   --->  ASCII+33\n Do you want to write out [y/n] "<<endl;
		}
		else
		{
			cout <<"PhredQ : ASCII+64   --->  ASCII+64\n Do you want to write out [y/n] "<<endl;
		}
		string YRun ="N" ;
		fd_set rfds;
		struct timeval tv;
		FD_ZERO(&rfds);
		FD_SET(0, &rfds);
		int maxfd =0 + 1;  
		int ret ;
		while(1) 
		{
			tv.tv_sec = 8;  
			tv.tv_usec = 0;
			ret = select(maxfd, &rfds, NULL, NULL, &tv);  
			if(ret<0)
			{
				printf("select error, process will eixt\n");  
				exit(0);  
			}  
			else if(FD_ISSET(STDIN_FILENO, &rfds))//测试是否有数据  
			{  
				char buf[100];
				fgets(buf, 100, stdin);
				YRun=buf[0];
				cout <<"your select is "<<YRun<<endl;
				break ;
			}
			else  
			{  
				cout <<"Time Out,default N "<<endl;
				break ;
			}  
		}

		if  (YRun == "N" ||  YRun == "n")
		{
			sam_close(in_SS);
			sam_close(outR1);
			bam_destroy1(aln_SS);
			bam_hdr_destroy(header_SS);
			string cc=" rm  -rf  "+ (para_Xam04->OutPut);
			std::system(cc.c_str()) ;
			cout <<"you can use linux command :\t ln -s "<<(para_Xam04->input_Sam)<<"\t"<<para_Xam04->OutPut<<endl;
			delete para_Xam04 ;
			return 0;
		}

		while  ( (sam_read1(in_SS, header_SS, aln_SS) >= 0) )
		{
			if ( (aln_SS->core).qual < (para_Xam04->minMapQ) ) 
			{
				continue ;
			}
			if  ((aln_SS->core).l_qseq   < ( para_Xam04->minLeng ))
			{
				continue ;
			}
			tmp=bam_write1(outR1->fp.bgzf, aln_SS);
		}
	}


	sam_close(in_SS);

	sam_close(outR1);
	bam_destroy1(aln_SS);
	bam_hdr_destroy(header_SS);

	delete para_Xam04 ;
	return 0;
}
#endif
///////// swimming in the sky and flying in the sea ////////////
