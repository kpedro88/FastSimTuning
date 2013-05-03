#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TSpectrum.h"

#define maxHBeta 16
#define maxHOeta 15
#define maxHEeta 14
#define maxHFeta 13
#define maxHDdet 4

void get_time_peaks(std::string fname, std::string oname="Etime_python.txt"){
	//open input root file
	TFile *file = TFile::Open(fname.c_str());
	
	//open output file
	std::ofstream output(oname.c_str());
	if (!output) {
		std::cerr << "Cannot open the output file " << oname << "\n";
		return;
	}
	output << std::fixed << std::setprecision(1); 
	
	TH1F* htmp;
	string dname[maxHDdet] = {"HB","HE","HO","HF"};
	std::stringstream hname;
	int maxHDeta[maxHDdet] = {maxHBeta, maxHEeta, maxHOeta, maxHFeta};
	//ieta shift values
	int ietaShift[maxHDdet] = {1,16,1,29};
	
	//get energy-weighted timing histos
	for(int d = 0; d < maxHDdet; d++){
		output << "ietaShift" << dname[d] << " = cms.int32(" << ietaShift[d] << ")," << std::endl;
		output << "timeShift" << dname[d] << " = cms.vdouble(";
		for(int i = 0; i < maxHDeta[d]; i++){
			hname.str(std::string());
			hname.clear();
			hname << "Etime_" << dname[d] << "_ieta" << i;
			htmp = (TH1F*)file->Get((hname.str()).c_str());
			
			//find peak
			TSpectrum *spec = new TSpectrum(5);
			spec->Search(htmp,6,"nodrawroot  goff");
			Float_t* xpos = spec->GetPositionX();
			Float_t* ypos = spec->GetPositionY();
			Double_t p = xpos[0];
			Double_t ph = ypos[0];
			
			if(d==3 && p < htmp->GetMean()) p = xpos[1];
			
			output << p;
			if(i<maxHDeta[d]-1) output << ", ";
			else output << "),";
		}
		output << std::endl;
	}	
	
}