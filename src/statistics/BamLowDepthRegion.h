#ifndef bamLDR_H_
#define bamLDR_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <list>
#include <map>
#include <iomanip>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include <cstdlib>
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <stdio.h>

using namespace std;
typedef unsigned long long ubit64_t;
//KSEQ_INIT(gzFile, gzread)


int  bamLDR_help()
{
	cout <<""
		"\n"
		"\tUsage: LowDepth  -InList     <bam.list>  -OutPut  <out.bed>\n"
		"\tUsage: LowDepth  -InDepthFa  <Depth.fa.gz>  -OutPut  <out.bed>\n"
		"\n"
		"\t\t-InList      <str>     Input Bam/Sam File List\n"
		"\t\t-InDepthF    <str>     In DepthFA File,bamCoverage OutFile\n"
		"\t\t-OutPut      <str>     OutPut Bed Region File\n"
		"\n"
		"\t\t-LowDepth    <int>     Regard SiteDepth < X as LowDepth[2]\n"
		"\t\t-MinLength   <int>     Filter too short region [1000]\n"
		"\t\t-MinQ        <int>     Ignore too low mapQ read[10]\n"
		"\n"
		"\t\t-help                  Show this help [hewm2008 v1.01]\n"
		"\n";
	return 1;
}

int bamLDR_help01(int argc, char **argv , In3str1v * paraFA04 )
{
	if (argc <=2 ) {bamLDR_help();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InList")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr1=argv[i];
		}
		else if (flag  ==  "InDepthF")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "LowDepth" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "MinQ" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt2=atoi(argv[i]);
		}
		else if (flag  ==  "MinLength" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InF=atof(argv[i]);
		}
		else if (flag == "help")
		{
			bamLDR_help();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ( (paraFA04->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	if  ((paraFA04->InStr1).empty()   &&  (paraFA04->InStr3).empty())
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;		
	}
	(paraFA04->InStr2)=add_Asuffix(paraFA04->InStr2);
	return 1 ;
}


int bamLowDepth_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=2;
	paraFA04->InInt2=10;
	paraFA04->InF=1000;
	if ((bamLDR_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 1 ;
	}

	int MinLength=int(paraFA04->InF);

	ogzstream  OUT ((paraFA04->InStr2).c_str());

	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<(paraFA04->InStr2)<<endl;
		delete  paraFA04 ; return  1;
	}


	if  ((paraFA04->InStr3).empty())
	{

	}
	else
	{

		igzstream FAIN ((paraFA04->InStr3).c_str(),ifstream::in);

		bool S_E=false;
		int Start=0;
		int End=0;
		int Count=0;
		string chrName="";
		int Depth=0;
		getline(FAIN,chrName);
		chrName=chrName.substr(1);

		while(!FAIN.eof())
		{
			string  line ;
			getline(FAIN,line);
			if (line.length()<=0)  { continue  ; }
			if (line[0] == '>')
			{
				if (S_E)
				{
					if ((End-Start)>MinLength)
					{
						OUT<<chrName<<"\t"<<Start<<"\t"<<End<<endl;
					}
				}
				chrName=line.substr(1);
				S_E=false;
				Count=0;
			}
			else
			{
				istringstream isone (line,istringstream::in);
				while(isone>>Depth)
				{				 		
					Count++;
					if ( (!S_E)   &&  ( Depth > (paraFA04->InInt)) )
					{
						continue ;
					}
					else if ( (S_E)   &&  ( Depth  > (paraFA04->InInt))  )
					{
						S_E=false ;
						if ((End-Start)>MinLength)
						{
							OUT<<chrName<<"\t"<<Start<<"\t"<<End<<"\n";
						}
					}
					else if   ( (!S_E)   &&  ( Depth  < (paraFA04->InInt)) )
					{
						S_E=true ;
						Start=Count ;  End=Count;
					}
					else 
					{
						End=Count;
					}
				}
			}


		}

		if (S_E)
		{
			if ((End-Start)>MinLength)
			{
				OUT<<chrName<<"\t"<<Start<<"\t"<<End<<endl;
			}
		}

		OUT.close();
		FAIN.close();
		delete  paraFA04 ; return  0;
	}





	igzstream LISTA ((paraFA04->InStr1).c_str(),ifstream::in); // igzstream
	string  BamPath="";
	getline(LISTA,BamPath);
	LISTA.close();
	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);


	cout<<"begin new the memory ...\n";

	unsigned short int **depth = new unsigned short int*[(header->n_targets)]; //开辟行  
	for(int i = 0; i < (header->n_targets); i++)  
	{
		depth[i] = new unsigned short int [(header->target_len[i])]; //开辟列  
		for (int32_t j =0 ; j< (header->target_len[i]) ; j++)
		{
			depth[i][j]=0;
		}
	}

	cout<<"new the memory done"<<endl;


	igzstream LIST ((paraFA04->InStr1).c_str(),ifstream::in); // igzstream
	if (!LIST.good())
	{
		cerr << "open List error: "<<(paraFA04->InStr1)<<endl;
	}
	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0)  { continue  ; }

		cout <<"Begin Bam :"<<line<<endl;
		bam_hdr_t *headerA;
		bam1_t *aln = bam_init1();

		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);

		if ((header->n_targets)!=(headerA->n_targets))
		{
			cerr<<"Error1 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			for(int i = 0; i <(header->n_targets); i++)
			{
				delete[] depth[i];  
			}
			delete[] depth;
			delete paraFA04 ;
			return 1;
		}
		bool NoSameBam=false;
		for(int i = 0; i < (header->n_targets); i++) 
		{
			if (strcmp(header->target_name[i],headerA->target_name[i])!=0)
			{
				NoSameBam=true;
				cerr<<header->target_name[i]<<"\t"<<headerA->target_name[i]<<endl;
				break ;
			}
		}
		if (NoSameBam)
		{
			cerr<<"Error2 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			for(int i = 0; i <(header->n_targets); i++)
			{
				delete[] depth[i];  
			}
			delete[] depth;
			delete paraFA04 ;
			return 1;
		}

		while (sam_read1(InBam, header, aln) >= 0)
		{
			if ((aln->core).tid < 0) {continue ;}
			if ( (aln->core).qual < (paraFA04->InInt2))
			{
				continue ;
			}
			for (int32_t j=0 ; j< ((aln->core).l_qseq) ; j++)
			{
				int32_t  This=((aln->core).pos)+j;
				depth[(aln->core).tid][This]++;
			}
		}
		sam_close(InBam);
		bam_destroy1(aln);
		bam_hdr_destroy(headerA);
	}
	LIST.close();

	cout <<"ALL Bam Read done"<<endl;

	for(int i = 0; i <(header->n_targets); i++)
	{
		bool S_E=false ;
		int Start=0;
		int End=0;
		for (int j =0 ; j< (header->target_len[i]) ; j++)
		{
			if ( (!S_E)   &&  ( depth[i][j] > (paraFA04->InInt)) )
			{
				continue ;
			}
			else if ( (S_E)   &&  ( depth[i][j] > (paraFA04->InInt))  )
			{
				S_E=false ;
				if ((End-Start)>MinLength)
				{
					OUT<<(header->target_name[i])<<"\t"<<Start<<"\t"<<End<<"\n";
				}
			}
			else if   ( (!S_E)   &&  ( depth[i][j] < (paraFA04->InInt)) )
			{
				S_E=true ;
				Start=j ;  End=j;
			}
			else 
			{
				End=j;
			}
		}
		if (S_E)
		{
			if ((End-Start)>MinLength)
			{
				OUT<<(header->target_name[i])<<"\t"<<Start<<"\t"<<End<<endl;
			}
		}
	}

	OUT.close();

	//释放开辟的资源  
	for(int i = 0; i <(header->n_targets); i++)
	{
		delete[] depth[i];  
	}
	delete[] depth;  


	bam_hdr_destroy(header);
	delete paraFA04 ;
	return 0;
}
#endif // bamLDR_H_  //
///////// swimming in the sky and flying in the sea ////////////





