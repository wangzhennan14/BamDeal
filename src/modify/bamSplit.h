#ifndef BamSplit_H_
#define BamSplit_H_
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
#include <htslib/kstring.h>
#include <stdio.h>

using namespace std;

int  BamSplit_help()
{
	cout <<""
		"\n"
		"\tUsage: bamSplit  -InList  <bam.list>  -ReSetHead \n"
		"\tUsage: bamSplit  -InFile A.bam  B.bam  \n"
		"\n"
		"\t\t-InFile    <str>     Input Bam/Sam File to split by chr\n"
		"\t\t-InList    <str>     Input Bam/Sam File List\n"
		"\n"
		"\t\t-OutDir    <str>     OutPut Dir for split [./]\n"
		"\t\t-MinQ      <int>     classify low mapQ<X read to unmap.bam[10]\n"
		"\t\t-OutSam              OutPut sam.gz File, not [bam]\n"
		"\t\t-ReSetHead           Reset Out.bam Header only with it's chr\n"
		"\n"
		"\t\t-help                Show this help [hewm2008 v1.02]\n"
		"\n";
	return 1;
}

int BamSplit_help01(int argc, char **argv , In3str1v * paraBamS )
{
	if (argc <=2 ) {BamSplit_help();return 0;}
	int file_count=0;
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
			string A=argv[i];
			file_count+=(ReadList(A ,(paraBamS->List)));
		}
		else if (flag  == "InFile")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			(paraBamS->List).push_back(A);
			file_count++;
			bool RunT=true;
			while(RunT)
			{
				if ( (i+1) < argc  && (argv[i+1][0]!='-'))
				{
					i++;
					A=argv[i];
					(paraBamS->List).push_back(A);
					file_count++;
				}
				else
				{
					RunT=false;
				}
			}
		}
		else if (flag  ==  "OutDir")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InStr2=argv[i];			
		}
		else if (flag  ==  "MinQ")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraBamS->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "ReSetHead")
		{
			paraBamS->TF=true ;
		}
		else if (flag  ==  "OutSam")
		{
			paraBamS->TF2=false ;
		}
		else if (flag == "help")
		{
			BamSplit_help();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((file_count<1) )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	return 1 ;
}


int Split2Sam( In3str1v * paraBamS )
{

	string  BamPath=(paraBamS->List)[0];

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);

	ogzstream *Soap2Chr = new ogzstream[(header->n_targets)];
	ogzstream UnMapOut ;
	string UUUtmp=(paraBamS->InStr2)+"/Unmap.sam.gz";
	UnMapOut.open(UUUtmp.c_str());
	for (int j=0; j<(header->n_targets) ; j++)
	{
		UnMapOut<<"@SQ\tSN:"<<(header->target_name[j])<<"\tLN:"<<(header->target_len[j])<<endl;
	}


	if  (!(paraBamS->TF))
	{

		for (int i=0; i<(header->n_targets) ; i++)
		{
			string aaa=(paraBamS->InStr2)+"/"+(header->target_name[i])+".sam.gz" ;
			Soap2Chr[i].open(aaa.c_str());
			if  (!Soap2Chr[i].good())
			{
				cerr<<"Can't open follow output:\n"<<aaa<<endl;
			}
			for (int j=0; j<(header->n_targets) ; j++)
			{
				Soap2Chr[i]<<"@SQ\tSN:"<<(header->target_name[j])<<"\tLN:"<<(header->target_len[j])<<endl;
			}
		}
	}
	else
	{
		for (int i=0; i<(header->n_targets) ; i++)
		{
			string aaa=(paraBamS->InStr2)+"/"+(header->target_name[i])+".sam.gz" ;
			Soap2Chr[i].open(aaa.c_str());
			if  (!Soap2Chr[i].good())
			{
				cerr<<"Can't open follow output:\n"<<aaa<<endl;
			}
			Soap2Chr[i]<<"@SQ\tSN:"<<(header->target_name[i])<<"\tLN:"<<(header->target_len[i])<<endl;
		}
	}

	int FileNum=(paraBamS->List).size();
	for (int ii=0 ; ii< FileNum ; ii++)
	{
		string  line =(paraBamS->List)[ii];
		if (line.length()<=0)  { continue  ; }		
		cout <<"Begin Bam :"<<line<<endl;
		bam_hdr_t *headerA;
		bam1_t *aln = bam_init1();

		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);

		if ((header->n_targets)!=(headerA->n_targets))
		{
			cerr<<"Error1 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			delete paraBamS ;
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
			delete paraBamS ;
			return 1;
		}

		while (sam_read1(InBam, header, aln) >= 0)
		{
			kstring_t Kstr = { 0, 0, NULL };
			int A=sam_format1(header, aln, &Kstr);

			if ( ((aln->core).tid < 0)  || (aln->core).qual < (paraBamS->InInt) )
			{
				UnMapOut<<Kstr.s<<endl;
			}
			else
			{
				Soap2Chr[(aln->core).tid]<<Kstr.s<<endl;
			}
		}
		sam_close(InBam);
		bam_destroy1(aln);
		bam_hdr_destroy(headerA);
	}

	bam_hdr_destroy(header);


	UnMapOut.close();
	for (int i=0; i<(header->n_targets) ; i++)
	{
		Soap2Chr[i].close() ;
	}

	delete [] Soap2Chr ;

	return 1 ;
}




int Split2Bam( In3str1v * paraBamS )
{

	string  BamPath=(paraBamS->List)[0];

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);

	htsFile **OUTBam = new  htsFile *[(header->n_targets)];
	bam_hdr_t **OUTheader=  new  bam_hdr_t *[(header->n_targets)];

	htsFile *UUM;
	string aaa=(paraBamS->InStr2)+"/UnMap.bam" ;
	UUM=hts_open(aaa.c_str(), "wb");
	int tmp=sam_hdr_write(UUM, header);


	for (int i=0; i<(header->n_targets) ; i++)
	{
		string aaa=(paraBamS->InStr2)+"/"+(header->target_name[i])+".bam" ;
		OUTBam[i]=hts_open(aaa.c_str(), "wb");
		OUTheader[i]=NULL;

		kstring_t str = { 0, 0, NULL };
		kputsn("@SQ\tSN:", 7, &str);
		kputs((header->target_name[i]), &str);
		kputsn("\tLN:", 4, &str);
		kputw(header->target_len[i], &str);
		kputc('\n', &str);

		OUTheader[i]= sam_hdr_parse(str.l, str.s);
		OUTheader[i]->l_text = str.l; OUTheader[i]->text = str.s;
		OUTheader[i]=sam_hdr_sanitise(OUTheader[i]);
		if (sam_hdr_write(OUTBam[i], OUTheader[i]) < 0)
		{
			cerr<<"Can't Write bam Header:\n"<<aaa<<endl;
			exit(1);
		}
	}
	//sam_close(OUTBam[0]);  return 1 ;

	int FileNum=(paraBamS->List).size();
	for (int ii=0 ; ii< FileNum ; ii++)
	{
		string  line =(paraBamS->List)[ii];
		if (line.length()<=0)  { continue  ; }

		cout <<"Begin Bam :"<<line<<endl;
		bam_hdr_t *headerA;
		bam1_t *aln = bam_init1();

		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);

		if ((header->n_targets)!=(headerA->n_targets))
		{
			cerr<<"Error1 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			delete paraBamS ;
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
			delete paraBamS ;
			return 1;
		}

		while (sam_read1(InBam, header, aln) >= 0)
		{
			int32_t i=(aln->core).tid ;
			if ( ((aln->core).qual < (paraBamS->InInt))  || (i < 0 ) )
			{
				int AA=bam_write1(UUM->fp.bgzf, aln);
			}			
			else
			{
				if ((aln->core).tid==(aln->core).mtid)
				{
					(aln->core).mtid=0;
				}
				(aln->core).tid=0;

				int AA=bam_write1(OUTBam[i]->fp.bgzf, aln);
			}
		}
		sam_close(InBam);
		bam_destroy1(aln);

		bam_hdr_destroy(headerA);
	}

	bam_hdr_destroy(header);

	sam_close(UUM);
	for (int i=0; i<(header->n_targets) ; i++)
	{
		sam_close(OUTBam[i]);
		bam_hdr_destroy(OUTheader[i]);
	}

	delete [] OUTheader ;
	delete [] OUTBam ;
	bam_hdr_destroy(header);
	return 1 ;
}





int Split2BamSameHeader( In3str1v * paraBamS )
{

	string  BamPath=(paraBamS->List)[0];

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	header = sam_hdr_read(BamIn);
	sam_close(BamIn);

	htsFile **OUTBam = new  htsFile *[(header->n_targets)];

	htsFile *UUM;
	string aaa=(paraBamS->InStr2)+"/UnMap.bam" ;
	UUM=hts_open(aaa.c_str(), "wb");
	int tmp=sam_hdr_write(UUM, header);

	for (int i=0; i<(header->n_targets) ; i++)
	{
		aaa=(paraBamS->InStr2)+"/"+(header->target_name[i])+".bam" ;
		OUTBam[i]=hts_open(aaa.c_str(), "wb");		
		if (sam_hdr_write(OUTBam[i], header) < 0)
		{
			cerr<<"Can't Write bam Header:\n"<<aaa<<endl;
			exit(1);
		}
	}
	//sam_close(OUTBam[0]);  return 1 ;

	int FileNum=(paraBamS->List).size();
	for (int ii=0 ; ii< FileNum ; ii++)
	{
		string  line =(paraBamS->List)[ii];
		if (line.length()<=0)  { continue  ; }

		cout <<"Begin Bam :"<<line<<endl;
		bam_hdr_t *headerA;
		bam1_t *aln = bam_init1();

		samFile *InBam = sam_open(line.c_str(), "r");
		headerA = sam_hdr_read(InBam);

		if ((header->n_targets)!=(headerA->n_targets))
		{
			cerr<<"Error1 :Bam header should be the same\n"<<BamPath<<"\n"<<line<<endl;
			delete paraBamS ;
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
			delete paraBamS ;
			return 1;
		}

		while (sam_read1(InBam, header, aln) >= 0)
		{
			int32_t i=(aln->core).tid ;
			if (( (aln->core).qual < (paraBamS->InInt)  )   ||  (i < 0 ))
			{
				int A=bam_write1(UUM->fp.bgzf, aln);
			}
			else
			{
				int A=bam_write1(OUTBam[i]->fp.bgzf, aln);
			}
		}
		sam_close(InBam);
		bam_destroy1(aln);

		bam_hdr_destroy(headerA);
	}

	bam_hdr_destroy(header);

	sam_close(UUM);

	for (int i=0; i<(header->n_targets) ; i++)
	{
		sam_close(OUTBam[i]);
	}

	delete [] OUTBam ;
	bam_hdr_destroy(header);
	return 1 ;
}








int bamSplit_main(int argc, char *argv[])
{
	In3str1v *paraBamS = new In3str1v;
	paraBamS->InInt=10;
	paraBamS->TF=false ;
	paraBamS->InStr2="./" ;
	if ((BamSplit_help01(argc, argv, paraBamS)==0))
	{
		delete paraBamS ;
		return 0 ;
	}

	if (paraBamS->TF2)
	{
		if  (paraBamS->TF)
		{
			Split2Bam( paraBamS );
		}
		else
		{
			Split2BamSameHeader( paraBamS );
		}
	}
	else
	{
		Split2Sam( paraBamS );
	}

	delete paraBamS ;
	return 0;
}
#endif // BamSplit_H_  //
///////// swimming in the sky and flying in the sea ////////////





