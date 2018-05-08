
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
#include <algorithm>
using namespace std;
typedef long long llong ;



int  print_Ausage_Xam2fa()
{
	cout <<""
		"\n"
		"\tUsage: bam2fa -InPut <in.bam> -OutPut <out.fa>\n"
		"\n"
		"\t\t-InPut     <str>   InPut sam/bam File\n"
		"\t\t-OutPut    <str>   OutPut fa File\n"
		"\n"
		"\t\t-UnMap             only output UnMapp read [NA]\n"
		"\n"
		"\t\t-help              show this help\n"
		"\n";
	//"\t\t-ShiftQ    <int>   Final phred quality [+31/0/-31]\n"
	return 1;
}

int parse_Acmd_Xam2fa(int argc, char **argv, In3str1v * para_Xam2fa)
{
	if (argc <=2 ) {print_Ausage_Xam2fa();return 0;}

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
			para_Xam2fa->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_Xam2fa->InStr2=argv[i];
		}
		else if (flag  ==  "UnMap")
		{
			para_Xam2fa->TF=false;
		}
		///*////
		else if (flag  == "help")
		{
			print_Ausage_Xam2fa();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_Xam2fa->InStr2).empty() ||   (para_Xam2fa->InStr1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	para_Xam2fa->InStr2=add_Asuffix(para_Xam2fa->InStr2);
	return 1 ;
}



//programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int Bam2fa_main(int argc, char **argv)
{
	In3str1v * para_Xam2fa = new In3str1v;
	if (parse_Acmd_Xam2fa(argc, argv , para_Xam2fa )==0)
	{
		delete para_Xam2fa ;
		return 0 ;
	}

	bam_hdr_t *header;
	bam1_t *aln = bam_init1();
	samFile *in = sam_open((para_Xam2fa->InStr1).c_str(), "r");
	header = sam_hdr_read(in);


	ogzstream OUT(para_Xam2fa->InStr2.c_str()) ;


	int Samp[256]={0};
	Samp['A']='T'; Samp['C']='G';
	Samp['T']='A'; Samp['G']='C';
	Samp['N']='N';

	int  Count=0;
	uint8_t Base[16] = {0,65,67,0,71,0,0,0,84,0,0,0,0,0,0,78};
	if ( para_Xam2fa->TF)
	{
		while  ( (sam_read1(in, header, aln) >= 0)    )
		{
			string readID=bam_get_qname(aln);

			if  ( (((aln)->core.flag&64) != 0))
			{
				readID= "@"+readID + "/1" ;
			}
			else if ((((aln)->core.flag&128) != 0) )
			{
				readID= "@"+readID + "/2" ;
			}
			else
			{
				readID= "@"+readID ;
			}


			uint8_t   *seq=bam_get_seq(aln);
			string BaseSeq="";
			for(int i=0; i < (aln->core).l_qseq ; i++)
			{
				BaseSeq=BaseSeq+(char)(Base[bam_seqi(seq, i)]);
			}




			if  ((((aln)->core.flag&16) != 0))
			{
				reverse(BaseSeq.begin(), BaseSeq.end());
				for (int i=0 ; i<(aln->core).l_qseq ; i++)
				{
					BaseSeq[i]=Samp[BaseSeq[i]];
				}
			}

			OUT<<">"<<readID<<"\n"<<BaseSeq<<"\n";


		}




	}



	else
	{



		while  ( (sam_read1(in, header, aln) >= 0)    )
		{
			if ((aln->core).tid >= 0)
			{
				continue ;
			}

			string readID=bam_get_qname(aln);

			if  ( (((aln)->core.flag&64) != 0))
			{
				readID= "@"+readID + "/1" ;
			}
			else if ((((aln)->core.flag&128) != 0) )
			{
				readID= "@"+readID + "/2" ;
			}
			else
			{
				readID= "@"+readID ;
			}



			uint8_t   *seq=bam_get_seq(aln);
			string BaseSeq="";
			for(int i=0; i < (aln->core).l_qseq ; i++)
			{
				BaseSeq=BaseSeq+Int2Str(Base[bam_seqi(seq, i)]);
			}




			if  ((((aln)->core.flag&16) != 0))
			{
				reverse(BaseSeq.begin(), BaseSeq.end());
				for (int i=0 ; i<(aln->core).l_qseq ; i++)
				{
					BaseSeq[i]=Samp[BaseSeq[i]];
				}
			}

			OUT<<">"<<readID<<"\n"<<BaseSeq<<"\n";


		}









	}
	sam_close(in);
	bam_destroy1(aln);
	bam_hdr_destroy(header);
	OUT.close();
	delete para_Xam2fa ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////


