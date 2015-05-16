#ifndef HAPLOTYPESET_H_INCLUDED
#define HAPLOTYPESET_H_INCLUDED

#include "StringBasics.h"
#include "VcfFileReader.h"
#include <fstream>
#include <sstream>
using namespace std;


class variant
{
public:

    string name;
    int bp;
    string chr;
    char refAllele,altAllele;
    string refAlleleString,altAlleleString;
    string MajAlleleString,MinAlleleString;

    variant(){};
    variant(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
    void assignValues(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
     void assignRefAlt(string &refe,string &alt)
    {
        refAlleleString=refe;
        altAlleleString=alt;
    };


};


class HaplotypeSet
{

	public:
		int         numHaplotypes;
		int         numMarkers;
		vector<int>        optEndPoints;
		vector<int>        ScaffoldIndex;
		vector<vector<char> >     haplotypesUnscaffolded;
		vector<vector<double> > alleleFreq;
		vector<vector<double> > Dosage;
		vector<vector<char> > ImputedAlleles;
		int PrintStartIndex,PrintEndIndex;
        bool GT,DS,GP;
        vector<int> SkipPositions;


		vector<string> markerName;
		vector<string> individualName;
		vector<char> refAlleleList,major, minor;
		vector<bool> missing, MarkerIndices;
		vector<variant> VariantList;

		bool allowMissing, vcfType,optmxType,machType;


        HaplotypeSet()
		{
			numHaplotypes = 0;
			numMarkers = 0;
			alleleFreq.clear();
			//haplotypes.clear();
			individualName.clear();
			markerName.clear();
			refAlleleList.clear();
			major.clear();
			minor.clear();
			missing.clear();
			MarkerIndices.clear();
			allowMissing = true;
			vcfType = false;
			optmxType=false;
			machType=false;
			PrintStartIndex=0;
			PrintEndIndex=0;



		}

        char    getScaffoldedHaplotype          (int sample,int marker);
        void    Create                          (vector<char> &tempHaplotype);
        void    calculateFreq                   ();
		bool    LoadSnpList                     (String filename);
		bool    LoadMachHaplotypes              (String &haps, String &snps);



        void    PrintVcfFile                    (String filename,bool gzip);
        char    convertAlleles                  (string markerId, string indivId, const char *alleles, string refAlleleString, string altAlleleString);
		bool    LoadHaplotypes                  (String filename, String snpNames, int maxIndiv, int maxMarker,bool optmx,String CNO,int START,int END);
		bool    FastLoadHaplotypes              (String filename, int maxIndiv, int maxMarker,String CNO,int START,int END,int WINDOW,bool rsid);
		bool    LoadTargetHaplotypes            (String &haps, String &snps);
		bool    LoadVcfTargetHaplotypes         (String filename, String snpNames, vector<string> &refSnpList, vector<variant> &refMarkerList);
        void    convertToReducedStructure       (vector<int> &optStructure);
        void    writeOptmFile                  (String &filename,bool &gzip);
        bool    readOptmFile                   (String optmFile,String CHR,int START,int END,int WINDOW);
        void    reconstructHaplotype            (vector<char> &reHaplotypes,int &index);
        void    SaveDosageForVcfOutput          (int hapID,vector<double> dose,vector<char> impAlleles);
        void    SaveDosageForVcfOutputSampleWise(int SamID,string &SampleName, vector<double> &dose1,vector<double> &dose2,vector<char> &impAlleles1,vector<char> &impAlleles2);
        void    InitializeDosageForVcfOutput    (int NHaps,int NMarkers);
        void    InitializePartialDosageForVcfOutput    (int NHaps,int NMarkers, vector<bool> &Format);
        void    PrintDosageForVcfOutputForID    (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);
        void    PrintPartialDosageForVcfOutputForID    (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);





};









#endif // HAPLOTYPESET_H_INCLUDED
