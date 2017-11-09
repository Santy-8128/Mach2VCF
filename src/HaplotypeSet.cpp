#include "HaplotypeSet.h"



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
            ifprintf(vcfOutFile,"\t%c|%c",a1,a2);
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

    IFILE ifs = ifopen(filename, "r");
    string line;
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
    char *end_str1;


    VariantList.resize(numMarkers);

    for(int i=0;i<numMarkers;i++)
    {

        string temp = (string) strtok_r ((char*)markerName[i].c_str(),":", &end_str1);
        if(temp.substr(0,3)=="chr")
        {
            temp=temp.erase(0,3);
        }
        VariantList[i].chr=temp;
        VariantList[i].bp=atoi(strtok_r (NULL, ":", &end_str1));

        stringstream strs;
        strs<<(VariantList[i].bp);

        VariantList[i].name=VariantList[i].chr+":"+(string)(strs.str());
        VariantList[i].refAlleleString="---";
        VariantList[i].altAlleleString="---";
        
    }
	return true;

}


bool HaplotypeSet::LoadMachHaplotypes(String &haps, String &snps)
{
	cout<<"\n Format = MaCH (Markov Chain based Haplotyper) "<<endl;


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

            char *end_str1;



            string temp = (string) strtok_r ((char*)line.c_str()," \t\n", &end_str1);
            string tempInd;

            tempInd = temp;
            if ((n-1) % 2 == 0)
            {
                individualName.push_back(tempInd);
            }


            string tempHap1;
            tempHap1 = (string) strtok_r (NULL," \t\n", &end_str1);
            string  tempHap = (string) strtok_r (NULL," \t\n", &end_str1);

            if(tempHap.length()!=markerName.size())
            {
                cerr<<" Error !!! Incorrect number of entries for HAPLO-"<< 1+((n-1) % 2)<<" Sample "<<individualName[(n-1)/2]<<" "<<endl;
                cerr<<" Number of variants is " <<tempHap.length() <<" while it should be "<<markerName.size()<<endl;
                cerr<<" Please check file "<<targetHapFile<<" !!! "<<endl<<endl;
                return false;
            }

            vector<char> tempHaps01(markerName.size(),0);
            
            for (int j = 0; j<(int)markerName.size() ;  j++)
            {
                
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
                    tempHap[j] = '0';
                    
                }

                char _a0_=VariantList[j].refAlleleString.c_str()[0];
                char _a1_=VariantList[j].altAlleleString.c_str()[0];
                
                bool setrefallele = true;
                setrefallele = setrefallele && VariantList[j].refAlleleString.compare("---")==0;
                setrefallele = setrefallele && tempHap[j]!='0';
                if (setrefallele)
                    VariantList[j].refAlleleString=tempHap.substr(j,1);
                bool setaltallele = true;
                setaltallele = setaltallele && (VariantList[j].altAlleleString.compare("---")==0 || _a0_==_a1_);
                setaltallele = setaltallele && tempHap[j]!='0';
                if (setaltallele)
                    VariantList[j].altAlleleString=tempHap.substr(j,1);

                if(tempHap[j]==VariantList[j].refAlleleString[0])
                {
                    tempHaps01[j]='0';
                }
                else if(tempHap[j]==VariantList[j].altAlleleString[0])
                {
                    tempHaps01[j]='1';
                }
                else
                {
                    tempHaps01[j]='.';
                }

            }

            haplotypesUnscaffolded.push_back(tempHaps01);

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






