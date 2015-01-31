#include "HaplotypeSet.h"
#include <algorithm>
#include <boost/algorithm/string.hpp>



void HaplotypeSet::PrintVcfFile(String filename,bool gzip)
{

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                                OUTPUT VCF FILE                                "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;


    std::cout << "\n Printing out to VCF File                            : " <<filename + ".vcf" + (gzip ? ".gz" : "") << endl;


    IFILE vcfOutFile=NULL;

    vcfOutFile = ifopen(filename + ".vcf" + (gzip ? ".gz" : ""), "wb");

    ifprintf(vcfOutFile,"##fileformat=VCFv4.1\n");
    time_t t = time(0);
    struct tm * now = localtime( & t );
    ifprintf(vcfOutFile,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
    ifprintf(vcfOutFile,"##source=Minimac3Convertor\n");
    ifprintf(vcfOutFile,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    ifprintf(vcfOutFile,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for(int hapId=0;hapId<(int)individualName.size();hapId++)
    {
        ifprintf(vcfOutFile,"\t%s",individualName[hapId].c_str());
    }
    int prevBp=0,duplicates=0;
    string prevRef="D",prevAlt="S";
    for (int i =0; i <numMarkers; i++)
   {


        if(prevBp==VariantList[i].bp && prevRef==VariantList[i].refAlleleString
           && prevAlt==VariantList[i].altAlleleString)

        {
//             cout<< VariantList[i].bp<<"\t"<<
//             VariantList[i].refAlleleString.c_str()<<"\t"<<
//             VariantList[i].altAlleleString.c_str()<<endl;
             duplicates++;
             continue;
        }
        else if(prevBp==VariantList[i].bp && prevAlt==VariantList[i].refAlleleString
           && prevRef==VariantList[i].altAlleleString)

        {
//             cout<< VariantList[i].bp<<"\t"<<
//             VariantList[i].refAlleleString.c_str()<<"\t"<<
//             VariantList[i].altAlleleString.c_str()<<endl;
             duplicates++;
               continue;
        }



        ifprintf(vcfOutFile,"\n%s\t%d\t%s\t%s\t%s\t.\tPASS\t.\tGT",
        VariantList[i].chr.c_str(),VariantList[i].bp,
        VariantList[i].name.c_str(),VariantList[i].refAlleleString.c_str(),
        VariantList[i].altAlleleString.c_str());

        prevBp=VariantList[i].bp;
        prevRef=VariantList[i].refAlleleString;
        prevAlt=VariantList[i].altAlleleString;


        for(int hapId=0;hapId<(int)individualName.size();hapId++)
        {
            char a1=haplotypesUnscaffolded[2*hapId][i];
            char a2=haplotypesUnscaffolded[2*hapId+1][i];
            ifprintf(vcfOutFile,"\t%d|%d",a1,a2);
        }
    }

        ifclose(vcfOutFile);


	std::cout << "\n Number of Markers Recorded                          : " << markerName.size() << endl;
	std::cout << " Number of Duplicate CHR:POS:REF:ALT found Recorded  : " << duplicates << "\n";
    std::cout << " Final Number of Markers written to VCF File         : " << markerName.size() - duplicates<< "\n";

    cout<<endl;


}

bool HaplotypeSet::LoadTargetHaplotypes(String &haps, String &snps)
{

 return LoadMachHaplotypes(haps,snps);

}


bool HaplotypeSet::LoadSnpList(String filename)
{

	std::cout << "\n Loading Marker List from File                       : " << filename << endl;
    typedef boost::tokenizer< boost::char_separator<char> > wsTokenizer;

    IFILE ifs = ifopen(filename, "r");
    string line;
//    SkipPositions
    if(ifs)
    {
        while ((ifs->readLine(line))!=-1)
        {
            markerName.push_back(line);
            line.clear();

        }
    }
    else
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl;
        return false;
    }
    ifclose(ifs);

    numMarkers = markerName.size();

    VariantList.resize(numMarkers);
    boost::char_separator<char> sep(":");

    for(int i=0;i<numMarkers;i++)
    {

        wsTokenizer t(markerName[i], sep);
        wsTokenizer::iterator ii=t.begin();

        string temp=(ii->c_str());
        if(temp.substr(0,3)=="chr")
        {
            temp=temp.erase(0,3);
        }
        VariantList[i].chr=temp;
        ++ii;
        VariantList[i].bp=atoi(ii->c_str());
        VariantList[i].name=VariantList[i].chr+":"+boost::lexical_cast<string>(VariantList[i].bp);

    }


	return true;

}


bool HaplotypeSet::LoadMachHaplotypes(String &haps, String &snps)
{
	typedef boost::tokenizer< boost::char_separator<char> > wsTokenizer;
    cout<<"\n Format = MaCH (MArkov Chain based Haplotyper) "<<endl;


    String targetSnpfile=snps;
    String targetHapFile=haps ;
	if(!LoadSnpList(targetSnpfile))
        return false;

	std::cout << "\n Number of Markers in Data                           : " << markerName.size() << endl<<endl;

	vector<int> knownPosition;
//	int counter = 0;


	int n = 0;

	IFILE ifs = ifopen(targetHapFile, "r");


	string line;
	if(ifs)
    {
        while ((ifs->readLine(line))!=-1)
        {
            n++;
            boost::char_separator<char> sep(" \t");
            wsTokenizer t(line, sep);
            wsTokenizer::iterator i;
            int col = 0;
            for (i = t.begin(); i != t.end(); ++i)
            {
                col++;

                if (col == 1)
                {

                    string tempInd;

                    tempInd = i->c_str();
                    if ((n-1) % 2 == 0)
                    {
                        individualName.push_back(tempInd);
                    }
                }
                if (col == 3)
                {
                    string tempHap;
                    tempHap = i->c_str();
                    if(tempHap.length()!=markerName.size())
                    {
                        cerr<<" Error !!! Incorrect number of entries for HAPLO-"<< 1+((n-1) % 2)<<" Sample "<<individualName[(n-1)/2]<<" "<<endl;
                        cerr<<" Number of variants is " <<tempHap.length() <<" while it should be "<<markerName.size()<<endl;
                        cerr<<" Please check file "<<targetHapFile<<" !!! "<<endl<<endl;
                        return false;
                    }




                    vector<char> tempHaps01(markerName.size(),0);
                    if(n==1)
                    {
                        for (int j = 0; j<(int)markerName.size() ;  j++)
                        {
//                            char allele=0;
                            switch (tempHap[j])
                            {

                            case 'A': case 'a': case '1': tempHap[j] = 'A'; break;
                            case 'C': case 'c': case '2': tempHap[j] = 'C'; break;
                            case 'G': case 'g': case '3': tempHap[j] = 'G'; break;
                            case 'T': case 't': case '4': tempHap[j] = 'T'; break;
                            case 'D': case 'd': case '5': tempHap[j] = 'D'; break;
                            case 'I': case 'i': case '6': tempHap[j] = 'I'; break;
                            case 'R': case 'r': case '7': tempHap[j] = 'R'; break;
                            default:
                                if (!allowMissing)
                                {

                                    cout << " Error: Haplotypes can only contain alleles A ('A', 'a' or '1'),\n";
                                    cout << " C ('C', 'c' or 2), G ('G', 'g' or 3), T ('T', 't' or '4'),\n";
                                    cout << " D ('D', 'd' or 5), I ('I', 'i' or 6) and R ('R', 'r' or 7).\n";
                                    cout << "\n\n For Individual : " << individualName.back() << ", Haplotype #";
                                    cout << (n-1) % 2 +1 << " has allele \"" << tempHap[j] << "\" at position : ";
                                    cout << j + 1 << endl;


                                    return false;
                                }
                            }

                            VariantList[j].refAlleleString=tempHap.substr(j,1);
                            VariantList[j].altAlleleString=tempHap.substr(j,1);
//                            cout<<VariantList[j].refAlleleString<<endl;
                        }
                        haplotypesUnscaffolded.push_back(tempHaps01);
                    }
                    else
                    {
                            for (int j = 0; j<(int)markerName.size() ;  j++)
                            {
//                                char allele=0;
                                switch (tempHap[j])
                                {

                                case 'A': case 'a': case '1': tempHap[j] = 'A'; break;
                                case 'C': case 'c': case '2': tempHap[j] = 'C'; break;
                                case 'G': case 'g': case '3': tempHap[j] = 'G'; break;
                                case 'T': case 't': case '4': tempHap[j] = 'T'; break;
                                case 'D': case 'd': case '5': tempHap[j] = 'D'; break;
                                case 'I': case 'i': case '6': tempHap[j] = 'I'; break;
                                case 'R': case 'r': case '7': tempHap[j] = 'R'; break;
                                default:
                                    if (!allowMissing)
                                    {

                                        cout << " Error: Haplotypes can only contain alleles A ('A', 'a' or '1'),\n";
                                        cout << " C ('C', 'c' or 2), G ('G', 'g' or 3), T ('T', 't' or '4'),\n";
                                        cout << " D ('D', 'd' or 5), I ('I', 'i' or 6) and R ('R', 'r' or 7).\n";
                                        cout << "\n\n For Individual : " << individualName.back() << ", Haplotype #";
                                        cout << (n-1) % 2 +1 << " has allele \"" << tempHap[j] << "\" at position : ";
                                        cout << j + 1 << endl;


                                        return false;
                                    }
                                }

                                if(tempHap[j]==VariantList[j].refAlleleString[0])
                                {
                                    tempHaps01[j]=0;
                                }
                                else if(tempHap[j]==VariantList[j].altAlleleString[0])
                                {
                                    tempHaps01[j]=1;
                                }
                                else
                                {
                                    if(VariantList[j].refAlleleString[0]==VariantList[j].altAlleleString[0])
                                    {
                                        tempHaps01[j]=1;
                                        VariantList[j].altAlleleString[0]=tempHap[j];
                                    }
                                    else
                                    {
                                        cerr<<" Error !!! Only bi-allelic variants supported. \n";
                                        cerr<<" "<<VariantList[j].name<<" has more than 2 alleles : "<<VariantList[j].refAlleleString ;
                                        cout<<" , "<< VariantList[j].altAlleleString <<" , "<< tempHap[j]<<endl;
                                        cerr<<" Please check file "<<targetHapFile<<" !!! "<<endl<<endl;
                                        return false;
                                    }


                                }

                            }

                            haplotypesUnscaffolded.push_back(tempHaps01);

                        }




                }
            }
            line.clear();
        }
    }
    else
    {
        cout<<"\n\n Following File File Not Available : "<<targetHapFile<<endl;
        return false;
    }

	std::cout << " Number of Markers Recorded                          : " << markerName.size() << endl;
	std::cout << " Number of Haplotypes Recorded                       : " << (haplotypesUnscaffolded.size()) << "\n";


    numHaplotypes = haplotypesUnscaffolded.size();
    numMarkers = markerName.size();


    if(numHaplotypes==0)
    {
        cout << "\n No haplotypes recorded from MaCH Input File : "<<targetHapFile<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }
    if(numHaplotypes%2==1)
    {
        cout << "\n Odd Number haplotypes recorded from MaCH Input File : "<<targetHapFile<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }

	std::cout << "\n Haplotype Set successfully loaded from MaCH File    : " << haps << endl;
	ifclose(ifs);

	return true;
}






