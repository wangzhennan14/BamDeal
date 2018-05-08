#ifndef BamSubChr_H_
#define BamSubChr_H_

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
#include <htslib/kstring.h>
#include <cstdlib>
#include <sys/select.h>

using namespace std;
typedef long long llong ;

int  print_SubChr_Usage()
{
	cout <<""
		"\n"
		"\tUsage: SubChr -InPut <in.bam> -RmList  <scaf.list> -OutPut <out.bam>\n"
		"\n"
		"\t\t-InPut       <str>    InPut sam/bam File\n"
		"\t\t-OutPut      <str>    OutPut bam\n"
		"\n"
		"\t\t-KeepList    <str>    Only keep some Chr in bam\n"
		"\t\t-RmList      <str>    Remove some chr in bam\n"
		"\n"
		"\t\t-RmUnmap              Remove UnMapping Read\n"
		"\t\t-ReSetHead            Reset Out.bam Header\n"
		"\n"
		"\t\t-help                 show this help\n" 
		"\n";
	return 1;
}

int parse_Acmd_SubChr(int argc, char **argv, In3str1v * para_Xam04)
{
	if (argc <=2 ) {print_SubChr_Usage();return 0;}

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
			para_Xam04->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_Xam04->InStr2=argv[i];
		}
		else if (flag  ==  "KeepList")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam04->InStr3=(argv[i]);
			para_Xam04->InInt=(para_Xam04->InInt) | 0x1;
		}
		else if (flag  ==  "RmList")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam04->InStr3=(argv[i]);
			para_Xam04->InInt=(para_Xam04->InInt) | 0x2;
		}
		else if (flag  == "ReSetHead")
		{
			para_Xam04->TF=false;
		}
		else if (flag  == "RmUnmap")
		{
			para_Xam04->TF2=false;
		}

		else if (flag  == "help")
		{
			print_SubChr_Usage();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ((para_Xam04->InStr1).empty() ||   (para_Xam04->InStr2).empty()  ||  (para_Xam04->InStr3).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;    
	}

	if (para_Xam04->InInt==3)
	{
		cerr<< "Para [ -RmList ]  and [ -KeepList ] can not use together"<<endl;
		return 0;
	}
	else if   (para_Xam04->InInt==0)
	{
		cerr<< "One of Para [ -RmList ]  and [ -KeepList ] must be seted"<<endl;
		return 0;
	}

	string path=(para_Xam04->InStr2);

	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	if (ext != "bam")
	{
		(para_Xam04->InStr2)=path+".bam";
	}

	return 1 ;


}


//programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int bam_SubChr_main(int argc, char **argv)
{
	In3str1v * para_Xam04 = new In3str1v;
	if (parse_Acmd_SubChr(argc, argv , para_Xam04 )==0)
	{
		delete para_Xam04; 
		return 0;
	}

	samFile *outR1 = sam_open(para_Xam04->InStr2.c_str(), "wb");



	bam_hdr_t *header_SS;
	bam1_t *aln_SS = bam_init1();
	samFile *in_SS = sam_open((para_Xam04->InStr1).c_str(), "r");


	header_SS = sam_hdr_read(in_SS);


	igzstream LIST (para_Xam04->InStr3.c_str(),ifstream::in); // igzstream

	int soapfilecout=0;
	if (!LIST.good())
	{
		cerr << "open List error: "<<para_Xam04->InStr3<<endl;
	}

	map <string ,int > ChrInfo;


	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0)  { continue  ; }
		ChrInfo.insert( map <string ,int >  :: value_type (line,-1));	
		soapfilecout++;
	}
	LIST.close();


	string chrID ;
	int tmp=1;

	if ( para_Xam04->TF )
	{
		if (sam_hdr_write(outR1, header_SS) < 0) 
		{
			fprintf(stderr, "Error writing output.\n");
			exit(-1);
		}


		if ( (para_Xam04->InInt)==1 )
		{
			while((sam_read1(in_SS, header_SS, aln_SS) >= 0))
			{
				if ((aln_SS->core).tid >= 0)
				{ // chr
					chrID=header_SS->target_name[(aln_SS->core).tid];
					if (ChrInfo.find(chrID)!=ChrInfo.end())
					{
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
			}
		}
		else
		{
			while((sam_read1(in_SS, header_SS, aln_SS) >= 0))
			{
				if ((aln_SS->core).tid >= 0)
				{ // chr
					chrID=header_SS->target_name[(aln_SS->core).tid];
					if (ChrInfo.find(chrID)==ChrInfo.end())
					{
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
				else
				{
					if (para_Xam04->TF2)
					{
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
			}
		}

	}


	else
	{


		kstring_t str = { 0, 0, NULL };
		map <string , int > :: iterator MapIt ;
		int Flag=0;

		if ( (para_Xam04->InInt)==1)
		{
			for(int i = 0; i < (header_SS->n_targets); i++)
			{
				string chr=header_SS->target_name[i];
				int TR=header_SS->target_len[i];
				MapIt=ChrInfo.find(chr);
				if (MapIt!=ChrInfo.end())
				{
					kputsn("@SQ\tSN:", 7, &str);
					kputs((MapIt->first).c_str(), &str);
					kputsn("\tLN:", 4, &str);
					kputw(TR, &str);
					kputc('\n', &str);
					MapIt->second= Flag ;
					Flag++;
				}
			}



			bam_hdr_t *NewHeader;
			NewHeader= sam_hdr_parse(str.l, str.s);
			NewHeader->l_text = str.l; NewHeader->text = str.s;
			NewHeader=sam_hdr_sanitise(NewHeader);
			tmp=sam_hdr_write(outR1, NewHeader);
			bam_hdr_destroy(NewHeader);

			while((sam_read1(in_SS, header_SS, aln_SS) >= 0))
			{
				if ((aln_SS->core).tid >= 0)
				{ // chr
					chrID=header_SS->target_name[(aln_SS->core).tid];
					MapIt=ChrInfo.find(chrID);
					if (MapIt!=ChrInfo.end())
					{
						(aln_SS->core).tid=MapIt->second;
						if ((aln_SS->core).mtid>=0)
						{				
							chrID=header_SS->target_name[(aln_SS->core).mtid];
							MapIt=ChrInfo.find(chrID);
							if (MapIt!=ChrInfo.end())
							{
							(aln_SS->core).mtid=MapIt->second;
							}
							else
							{
								(aln_SS->core).mtid==-1;
							}
						}
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
			}



		}
		else
		{



			for(int i = 0; i < (header_SS->n_targets); i++)
			{
				string chr=header_SS->target_name[i];
				int TR=header_SS->target_len[i];
				MapIt=ChrInfo.find(chr);
				if (MapIt==ChrInfo.end())
				{
					kputsn("@SQ\tSN:", 7, &str);
					kputs((MapIt->first).c_str(), &str);
					kputsn("\tLN:", 4, &str);
					kputw(TR, &str);
					kputc('\n', &str);
					MapIt->second= Flag ;
					Flag++;
				}
			}


			bam_hdr_t *NewHeader; 

			NewHeader= sam_hdr_parse(str.l, str.s);
			NewHeader->l_text = str.l; NewHeader->text = str.s;
			NewHeader=sam_hdr_sanitise(NewHeader);
			tmp=sam_hdr_write(outR1, NewHeader);
			bam_hdr_destroy(NewHeader);




			while((sam_read1(in_SS, header_SS, aln_SS) >= 0))
			{
				if ((aln_SS->core).tid >= 0)
				{ // chr
					chrID=header_SS->target_name[(aln_SS->core).tid];

					MapIt=ChrInfo.find(chrID);
					if (MapIt==ChrInfo.end())
					{
						(aln_SS->core).tid=MapIt->second;
						chrID=header_SS->target_name[(aln_SS->core).mtid];
						if ((aln_SS->core).mtid>=0)
						{				
							chrID=header_SS->target_name[(aln_SS->core).mtid];
							MapIt=ChrInfo.find(chrID);
							if (MapIt!=ChrInfo.end())
							{
							(aln_SS->core).mtid=MapIt->second;
							}
							else
							{
								(aln_SS->core).mtid==-1;
							}
						}
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
				else
				{
					if (para_Xam04->TF2)
					{
						if ((aln_SS->core).mtid>=0)
						{				
							chrID=header_SS->target_name[(aln_SS->core).mtid];
							MapIt=ChrInfo.find(chrID);
							if (MapIt!=ChrInfo.end())
							{
							(aln_SS->core).mtid=MapIt->second;
							}
							else
							{
								(aln_SS->core).mtid==-1;
							}
						}
						tmp=bam_write1(outR1->fp.bgzf, aln_SS);
					}
				}
			}
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
