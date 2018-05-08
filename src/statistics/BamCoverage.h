#ifndef bamCov_H_
#define bamCov_H_
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


int  bamCov_help()
{
	cout <<""
		"\n"
		"\tUsage: Coverage  -List  <bam.list>  -OutPut  <out.depth.gz>\n"
		"\n"
		"\t\t-List      <str>     Input Bam/Sam File List\n"
		"\t\t-OutPut    <str>     OutPut Depth info Fa file\n"
		"\n"
		"\t\t-Ref       <str>     In Ref.fa If Want Out Depth-GC wig info\n"
		"\t\t-Windows   <int>     Windows size for Depth-GC wig[10000]\n"
		"\t\t-Stat                Stat Coverage,MeanDepth & Depth Dis info\n"
		"\n"
		"\t\t-MinQ      <int>     Filter the read low mapQ[10]\n"
		"\n"
		"\t\t-help                Show this help [hewm2008 v1.10]\n"
		"\n";
	return 1;
}

int bamCov_help01(int argc, char **argv , In3str1v * paraFA04 )
{
	if (argc <=2 ) {bamCov_help();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "List")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "Ref")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  ==  "Windows")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InF=atof(argv[i]);
		}
		else if (flag  ==  "MinQ")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "Stat")
		{
			paraFA04->TF=false ;
		}
		else if (flag == "help")
		{
			bamCov_help();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((paraFA04->InStr1).empty() || (paraFA04->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	(paraFA04->InStr2)=add_Asuffix(paraFA04->InStr2);
	return 1 ;
}


int bamCov_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=10;
	paraFA04->InF=10000;
	if ((bamCov_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 0 ;
	}


	ogzstream  OUT ((paraFA04->InStr2).c_str());

	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<(paraFA04->InStr2)<<endl;
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
			if ( (aln->core).qual < (paraFA04->InInt))
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

	if (!(paraFA04->InStr3).empty())
	{
		string GCDepth=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-3);
		GCDepth=GCDepth+".DepthGC.wig.gz";
		ogzstream  OUTGC (GCDepth.c_str());
		OUTGC<<"##Bigin\tDepth(X)\tGC(%)\n";

		gzFile fp;
		kseq_t *seq;
		int l;
		fp = gzopen((paraFA04->InStr3).c_str(), "r");
		seq = kseq_init(fp);
		int Windows=int(paraFA04->InF);
		int cutOffLen=int(Windows*0.4);
		while ((l = kseq_read(seq)) >= 0)
		{
			int  NumIntSca=0;
			for(NumIntSca=0; NumIntSca <(header->n_targets); NumIntSca++)
			{
				if ( strcmp(header->target_name[NumIntSca],(seq->name.s)) ==0)
				{
					break;
				}
			}
			int  LastDis=(header->target_len[NumIntSca])-(Windows);
			if (LastDis<0) {continue;}
			OUTGC<<"#>"<<(seq->name.s)<<"\tWindowsSite:\t"<<(Windows)<<endl;

			for (int j =0 ; j< LastDis ; j+=(Windows))
			{
				int  Ascii[256] = {0};
				ubit64_t TotalDepth=0;
				for (int kk=0 ; kk<(Windows) ; kk++)
				{
					int site=j+kk;
					Ascii[seq->seq.s[site]]++;
					TotalDepth+=depth[NumIntSca][site];
				}
				int  N_Aleng= Ascii['N']+Ascii['n'] ;
				int  eff_Acout=(Windows)-N_Aleng ;
				if  (eff_Acout<cutOffLen)
				{
					OUTGC<<(j+1)<<"\tNA\tNA"<<endl;
					continue ;
				}
				int  GC_Aleng=Ascii['C']+Ascii['G'] + Ascii['c']+Ascii['g'] ;
				double MeanDepth=TotalDepth*1.0/eff_Acout;
				double GC_rio=GC_Aleng*100.0/eff_Acout;

				OUTGC<<(j+1)<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<MeanDepth<<"\t"<<GC_rio<<endl;

			}
		}

		OUTGC.close();
		kseq_destroy(seq);
		gzclose(fp);

	}


	if (paraFA04->TF)
	{
		for(int i = 0; i <(header->n_targets); i++)
		{
			OUT<<">"<<(header->target_name[i]);		
			for (int j =0 ; j< (header->target_len[i]) ; j++)
			{
				if (j%50000==0)
				{
					OUT<<endl<<depth[i][j];
				}
				else
				{
					OUT<<" "<<depth[i][j];
				}
			}
			OUT<<endl ;
		}
		//		OUT<<endl ;
	}
	else
	{

		string PrefixO=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-3);
		string StatF=PrefixO+".stat";
		string PLot=PrefixO+".plot";

		ofstream  OUTPLOT (PLot.c_str());
		ofstream  OUTSTAT (StatF.c_str());

		OUTSTAT<<"#ID\tLength\tCoveredBase\tTotalDepth\tCoverage%\tMeanDepth"<<endl;
		OUTPLOT<<"#Depth\tNumber"<<endl;

		map <unsigned short int,ubit64_t> MapDepthSum ;
		map <unsigned short int,ubit64_t> :: iterator Th_mapiT ;
		ubit64_t  Sum_chrLen=0;
		ubit64_t Sum_CoverBase=0;
		ubit64_t Sum_Depth=0;

		for(int i = 0; i <(header->n_targets); i++)
		{

			OUT<<">"<<(header->target_name[i]);		
			map <unsigned short int,ubit64_t> MapDepthChr ;
			map <unsigned short int,ubit64_t> :: iterator mapIt_Depth ;

			for (int j =0 ; j< (header->target_len[i]) ; j++)
			{

				mapIt_Depth=MapDepthChr.find(depth[i][j]);
				if  (mapIt_Depth==MapDepthChr.end())
				{
					MapDepthChr.insert( map <unsigned short int,ubit64_t>  :: value_type (depth[i][j],1));
				}
				else
				{
					(mapIt_Depth->second)++;
				}

				if (j%50000==0)
				{
					OUT<<endl<<depth[i][j];
				}
				else
				{
					OUT<<" "<<depth[i][j];
				}
			}
			OUT<<endl ;
			ubit64_t CovBase=0;
			ubit64_t TotalDepth=0;
			for (mapIt_Depth=MapDepthChr.begin(); mapIt_Depth!=MapDepthChr.end(); mapIt_Depth++)
			{
				if ((mapIt_Depth->first)!=0)
				{
					CovBase+=mapIt_Depth->second;
				}
				Th_mapiT=MapDepthSum.find(mapIt_Depth->first);
				TotalDepth+=(mapIt_Depth->second)*(mapIt_Depth->first);
				if (Th_mapiT==MapDepthSum.end())
				{
					MapDepthSum.insert( map <unsigned short int,ubit64_t>  :: value_type (mapIt_Depth->first,mapIt_Depth->second));
				}
				else
				{
					(Th_mapiT->second)+=(mapIt_Depth->second);
				}
			}

			double Coverage=CovBase*100.0/(header->target_len[i]);
			double MeanDepth=TotalDepth*1.0/(header->target_len[i]);
			Sum_chrLen+=(header->target_len[i]);
			Sum_CoverBase+=CovBase ;
			Sum_Depth+=TotalDepth;
			OUTSTAT<<(header->target_name[i])<<"\t"<<(header->target_len[i])<<"\t"<<CovBase<<"\t"<<TotalDepth<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<Coverage<<"\t"<<MeanDepth<<endl;

		}
		double Coverage=Sum_CoverBase*100.0/Sum_chrLen ;
		double MeanDepth=Sum_Depth*1.0/Sum_chrLen ;

		OUTSTAT<<"#Genome\t"<<Sum_chrLen<<"\t"<<Sum_CoverBase<<"\t"<<Sum_Depth<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<Coverage<<"\t"<<MeanDepth<<endl;

		for (Th_mapiT=MapDepthSum.begin(); Th_mapiT!=MapDepthSum.end() ;Th_mapiT++)
		{
			OUTPLOT<<Th_mapiT->first<<"\t"<<Th_mapiT->second<<endl;
		}


		OUTSTAT.close();
		OUTPLOT.close();

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
#endif // bamCov_H_  //
///////// swimming in the sky and flying in the sea ////////////





