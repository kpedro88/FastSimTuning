#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TPolyMarker.h"
#include "TLine.h"
#include "TStyle.h"
#include "TDirectory.h"

#define maxHDe 16 //energy points for hadrons
#define maxHDeta 50 //eta bins for hadrons
#define nCballPar 7 //#pars for cball fn

using namespace TMath;

//-------------------------------------
//class to store energy resolution data
class energyRes {
	//vars
	private:
		Double_t energy;
		Int_t imip;
		Int_t ieta;
		std::vector<Double_t> stats; //stats are N, mean, RMS
		std::vector<Double_t> stat_errs;
		TF1 *fit;

	public:
		//constructors
		energyRes(Double_t en, Int_t im, Int_t ie) : energy(en), imip(im), ieta(ie) {}
		
		//set members
		void setStats(std::vector<Double_t>& _stats, std::vector<Double_t>& _stat_errs){
			stats = _stats;
			stat_errs = _stat_errs;
		}
		void setFit(TF1* _fit) { fit = _fit; }
		void setEnergy(Double_t en) { energy = en; }
		void setEta(Int_t ie) { ieta = ie; }
		void setMip(Int_t im) { imip = im; }
	
		//access members
		Double_t getEnergy() { return energy; }
		Int_t getEta() { return ieta; }
		Int_t getMip() { return imip; }
		Double_t getStat(int istat) { return stats[istat]; }
		Double_t getStatErr(int istat) { return stat_errs[istat]; }
		TF1* getFit() { return fit; }

};


//------------------------------------------
//Double-sided Crystal Ball function
//parameters:
//N, mu, sigma, aL, nL, aR, nR
//0,  1,     2,  3,  4,  5,  6
Double_t cball(Double_t *x, Double_t *par){
	//ensure sigma > 0 and a > 0
	//Double_t N = 1/(sigma*(n/a*1/(n-1)*Exp(-a*a/2) + Sqrt(Pi()/2)*(1+Erf(a/Sqrt(2))))); //normalized N
	Double_t N = par[0]; //let N float
	Double_t mu = par[1];
	par[2] = Abs(par[2]);
	Double_t sigma = par[2];
	par[3] = Abs(par[3]);
	Double_t aL = par[3];
	par[4] = (par[4]>1) ? par[4] : 1.01; //n>1 required
	Double_t nL = par[4];
	par[5] = Abs(par[5]);
	Double_t aR = par[5];
	par[6] = (par[6]>1) ? par[6] : 1.01; //n>1 required
	Double_t nR = par[6];	
	Double_t arg = (x[0]-mu)/sigma;
	
	//left tail
	if(arg <= -aL){
		return N*Power(nL/aL,nL)*Exp(-Power(aL,2)/2)*Power(nL/aL-aL-arg,-nL);
	}
	//right tail
	else if(arg >= aR){
		return N*Power(nR/aR,nR)*Exp(-Power(aR,2)/2)*Power(nR/aR-aR+arg,-nR);
	}
	//core
	else{
		return N*Exp(-Power(arg,2)/2);
	}

}

//----------------------------------
//function to fit energy resolutions
energyRes* get_res(int num, int imip, int ieta, bool do_fit, bool do_show, bool do_print, bool do_batch=false){
	//energy values
	Double_t energies[] = {1., 2., 3., 5., 9., 11., 15., 20., 30., 50., 100., 150., 225., 300., 1000., 3000. };

	if (num>=maxHDe || num<0) {
		std::cout << "num must be between 0 and " << maxHDe - 1 << std::endl;
		energyRes* theRes = new energyRes(0,0,0);
		return theRes;
	}

	if (ieta>maxHDeta || ieta<1) {
		std::cout << "ieta must be between 1 and " << maxHDeta << std::endl;
		energyRes* theRes = new energyRes(0,0,0);
		return theRes;
	}
	
	if (imip>2 || imip<0) {
		std::cout << "imip must be between 0 and 2" << std::endl;
		energyRes* theRes = new energyRes(0,0,0);
		return theRes;
	}

	//make filenames
	std::stringstream fname, piname;
	fname << "tree_pion_" << energies[num] << "_500k.root";
	piname << "#pi^{-} " << energies[num] << " GeV, " << (ieta-1)*0.1 << "#leq|#eta|<" << ieta*0.1;

	//open file and tree
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TTree* totalTree = (TTree*)_file->Get("Total");

	//make tree drawing expressions
	//define mip as ecal < 1 gev = 1000 mev
	double mipcut = 0.8;
	std::stringstream drawname;
	std::stringstream cutname;
	std::stringstream etacut;
	//default histo settings
	double Emin = 0.1*energies[num]; //lower cut to remove weird peaks near E=zero
	double Emax = 2*energies[num];
	int nbins = 100;
	
	//ecal & hcal energies are already "calibrated" at sim level (sampling factors for HBHEHO, poisson smearing & pe conversion for HF)
	drawname << "(ecal+hcal)>>htemp(" << nbins << "," << Emin << "," << Emax << ")";

	etacut << (ieta-1)*0.1 << " <= abs(eta) && abs(eta) < " << ieta*0.1;
	//0 is mip, 1 is nomip, 2 is total
	if(imip==0) cutname << "ecal < " << mipcut << " && " << etacut.str();
	else if(imip==1) cutname << "ecal > " << mipcut << " && " << etacut.str();
	else cutname << etacut.str();

	TH1F* h_res; //to store histos drawn from tree
	TF1* gfit;
	TF1* gsnL;
	TF1* gsnR;

	//plotting variables
	TCanvas* can;
	TPad* pad;
	TLegend* leg;
	TPaveText* pave;
	TPaveText* pave_par;
	TLine *aLline;
	TLine *aRline;

	//create instance of energyRes object
	energyRes* theRes = new energyRes(energies[num],imip,ieta);

	//draw w/ appropriate cut
	totalTree->Draw((drawname.str()).c_str(),(cutname.str()).c_str(),"hist goff");
	h_res = (TH1F*)gDirectory->Get("htemp");

	//names
	std::string omip, ofit;
	if (imip==0) omip = "mip";
	else if (imip==1) omip = "nomip";
	else omip = "tot";
	if(do_fit) ofit = "cballD";
	else ofit = "nofit";
	std::stringstream oname;
	oname << "pion_response_" << omip << "_" << ofit << "_" << energies[num] << "gev_ieta" << ieta;
	
	//get values from histo
	Double_t m = h_res->GetMean();
	Double_t me = h_res->GetMeanError();
	//Double_t m = h_res->GetBinCenter(h_res->GetMaximumBin()); //peak
	Double_t s = h_res->GetRMS();
	Double_t se = h_res->GetRMSError();
	Int_t N = h_res->GetEntries();
	
	std::vector<Double_t> stats(3,0);
	std::vector<Double_t> stat_err(3,0);
	stats[0] = N;
	stat_err[0] = 0;
	stats[1] = m;
	stat_err[1] = me;
	stats[2] = s;
	stat_err[2] = se;

	//find peak
	TSpectrum *spec = new TSpectrum(5);
	if(nbins < 100) spec->Search(h_res,6,"nobackground nodraw goff"); //turn off background removal when nbins too small
	else spec->Search(h_res,6,"nodraw goff");
	Float_t* xpos = spec->GetPositionX();
	Float_t* ypos = spec->GetPositionY();
	Double_t p = xpos[0];
	Double_t ph = ypos[0];
	if(do_show) std::cout << "peak: " << p << std::endl;
	
	//setup fitting function & do fit
	if (do_fit){
		gfit = new TF1((oname.str()).c_str(),cball,Emin,Emax,nCballPar);
		gfit->SetParameters(ph,p,s,1,1.1,1,1.1);
		
		//limits on parameters: 0 < a < 10, 1 < n < 200
		gfit->SetParLimits(3,0,10);
		gfit->SetParLimits(5,0,10);
		gfit->SetParLimits(4,1.01,200);
		gfit->SetParLimits(6,1.01,200);
		
		//formatting
		gfit->SetLineColor(kRed);
		gfit->SetMarkerColor(kRed);
		gfit->SetLineWidth(2);
		//fit
		h_res->Fit(gfit,"LNQRB");
	}
	
	//store parameters
	theRes->setStats(stats,stat_err);
	if(do_fit) theRes->setFit(gfit);
	
	std::stringstream muname, signame, aLname, nLname, aRname, nRname, Nname, chiname;
	muname.precision(2);
	signame.precision(2);
	aLname.precision(2);
	nLname.precision(2);
	aRname.precision(2);
	nRname.precision(2);
	chiname.precision(5);
	if (do_fit) {
		muname << fixed << "#mu = " << gfit->GetParameter(1) << " #pm " << gfit->GetParError(1);
		signame << fixed << "#sigma = " << gfit->GetParameter(2) << " #pm " << gfit->GetParError(2);
		aLname << fixed << "a_{L} = " << gfit->GetParameter(3) << " #pm " << gfit->GetParError(3);
		nLname << fixed << "n_{L} = " << gfit->GetParameter(4) << " #pm " << gfit->GetParError(4);
		aRname << fixed << "a_{R} = " << gfit->GetParameter(5) << " #pm " << gfit->GetParError(5);
		nRname << fixed << "n_{R} = " << gfit->GetParameter(6) << " #pm " << gfit->GetParError(6);
		chiname << fixed << "#chi^{2}/ndf = " << gfit->GetChisquare()/gfit->GetNDF();
	}
	else {
		muname << fixed << "Mean = " << m << " #pm " << me;
		signame << fixed << "RMS = " << s << " #pm " << se;	
	}
	Nname << "N = " << N; 

	//plotting
	if (do_show){
		can = new TCanvas(omip.c_str(),omip.c_str(),700,500);
		can->cd();
		pad = new TPad("graph","",0,0,1,1);
		pad->SetMargin(0.12,0.05,0.15,0.05);
		pad->Draw();
		pad->cd();
		
		//formatting
		h_res->SetStats(kTRUE);
		gStyle->SetOptStat("mr");
		h_res->SetTitle("");
		h_res->GetXaxis()->SetTitle("Energy [GeV]");
		h_res->SetLineWidth(2);
		h_res->SetLineColor(kBlack);
		h_res->GetYaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetYaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetXaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetXaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetYaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetXaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		
		//plot histo and fit
		h_res->Draw("hist");
		if(do_fit) gfit->Draw("same");	
	
		//determine placing of legend and paves - par pave goes on side with more space
		Double_t xmin;
		if (m/((h_res->GetXaxis()->GetXmax() + h_res->GetXaxis()->GetXmin())/2) < 1) xmin = 0.65;
		else xmin = 0.2;
	
		if(do_fit) { //legend
			leg = new TLegend(xmin,0.78,xmin+0.2,0.88);
			leg->AddEntry(h_res,"CMSSW");
			leg->AddEntry(gfit,"Fit");
			leg->SetFillColor(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.05);
			leg->SetTextFont(42);
			leg->Draw("same");
			
			can->Update();
			
			//left line
			Double_t bndL = gfit->GetParameter(1) - gfit->GetParameter(2)*gfit->GetParameter(3);
			aLline = new TLine(bndL,pad->GetUymin(),bndL,pad->GetUymax());
			aLline->SetLineStyle(2);
			aLline->SetLineWidth(3);
			aLline->SetLineColor(kBlue);
			aLline->Draw("same");
			
			//left gaussian
			gsnL = new TF1("gsn","gaus",Emin,bndL);
			gsnL->SetParameters(gfit->GetParameter(0),gfit->GetParameter(1),gfit->GetParameter(2));
			gsnL->SetLineColor(kRed);
			gsnL->SetMarkerColor(kRed);
			gsnL->SetLineWidth(2);
			gsnL->SetLineStyle(2);
			gsnL->Draw("same");

			//line
			Double_t bndR = gfit->GetParameter(1) + gfit->GetParameter(2)*gfit->GetParameter(5);
			aRline = new TLine(bndR,pad->GetUymin(),bndR,pad->GetUymax());
			aRline->SetLineStyle(2);
			aRline->SetLineWidth(3);
			aRline->SetLineColor(kBlue);
			aRline->Draw("same");
			
			//right gaussian
			gsnR = new TF1("gsn","gaus",bndR,Emax);
			gsnR->SetParameters(gfit->GetParameter(0),gfit->GetParameter(1),gfit->GetParameter(2));
			gsnR->SetLineColor(kRed);
			gsnR->SetMarkerColor(kRed);
			gsnR->SetLineWidth(2);
			gsnR->SetLineStyle(2);
			gsnR->Draw("same");			

		}
		
		//name
		std::stringstream mipname;
		if(imip==0) mipname << "mip, E_{ecal} < " << mipcut << " GeV";
		else if(imip==1) mipname << "nomip, E_{ecal} > " << mipcut << " GeV";
		else mipname << "total";
		
		//pave
		pave = new TPaveText(xmin,0.68,xmin+0.2,0.78,"NDC");
		pave->AddText((piname.str()).c_str());
		pave->AddText((mipname.str()).c_str());
		pave->SetFillColor(0);
		pave->SetBorderSize(0);
		pave->SetTextFont(42);
		pave->SetTextSize(0.05);
		pave->Draw("same");

		//par pave
		Double_t ymin1;
		if(do_fit) ymin1 = 0.26;
		else ymin1 = 0.51;
		pave_par = new TPaveText(xmin,ymin1,xmin+0.2,0.66,"NDC");
		pave_par->AddText((Nname.str()).c_str());
		pave_par->AddText((muname.str()).c_str());
		pave_par->AddText((signame.str()).c_str());
		if(do_fit){
			pave_par->AddText((aLname.str()).c_str());
			pave_par->AddText((nLname.str()).c_str());
			pave_par->AddText((aRname.str()).c_str());
			pave_par->AddText((nRname.str()).c_str());
			pave_par->AddText((chiname.str()).c_str());
		}
		pave_par->SetFillColor(0);
		pave_par->SetBorderSize(0);
		pave_par->SetTextFont(42);
		pave_par->SetTextSize(0.05);
		pave_par->Draw("same");
		
		if (imip==0) std::cout << "mip:" << std::endl;
		else if (imip==1) std::cout << "nomip:" << std::endl;
		else std::cout << "total:" << std::endl;
		if(do_fit){
			std::cout << "N = " << N << std::endl;
			std::cout << "mu = " << gfit->GetParameter(1) << " +/- " << gfit->GetParError(1) << std::endl;
			std::cout << "sigma = " << gfit->GetParameter(2) << " +/- " << gfit->GetParError(2) << std::endl;
			std::cout << "aL = " << gfit->GetParameter(3) << " +/- " << gfit->GetParError(3) << std::endl;
			std::cout << "nL = " << gfit->GetParameter(4) << " +/- " << gfit->GetParError(4) << std::endl;
			std::cout << "aR = " << gfit->GetParameter(5) << " +/- " << gfit->GetParError(5) << std::endl;
			std::cout << "nR = " << gfit->GetParameter(6) << " +/- " << gfit->GetParError(6) << std::endl;
			std::cout << "chi^2/ndf = " << gfit->GetChisquare()/gfit->GetNDF() << std::endl;
		}
		else{
			std::cout << std::endl;
			std::cout << "N = " << N << std::endl;
			std::cout << "Mean = " << m << std::endl;
			std::cout << "RMS = " << s << std::endl;
		}
	
		if(do_print) can->Print((oname.str()+".png").c_str(),"png");
		if(do_batch) _file->Close();
	}
	else { _file->Close(); }
	
	//return data structure with relevant info
	return theRes;
}

//-----------------------------------------
//function to store all fit results in tree
void allthefits(string outname="tree_cballD_lim.root"){
	//energy values
	Double_t energies[] = {1., 2., 3., 5., 9., 11., 15., 20., 30., 50., 100., 150., 225., 300., 1000., 3000. };	

	//tree variables
	double s_energy, s_mu, s_sigma, s_aL, s_nL, s_aR, s_nR, s_chi2, s_mu_e, s_sigma_e, s_aL_e, s_nL_e, s_aR_e, s_nR_e;
	int s_imip, s_ieta, s_N, s_ndf;

	//open file and tree
	TFile* out_file = new TFile(outname.c_str(), "RECREATE");
	TTree* s_tree = new TTree("tree", "pion fit parameters");
	
	//setup tree
	s_tree->Branch("energy",&s_energy,"s_energy/D");
	s_tree->Branch("mu",&s_mu,"s_mu/D");
	s_tree->Branch("sigma",&s_sigma,"s_sigma/D");
	s_tree->Branch("aL",&s_aL,"s_aL/D");
	s_tree->Branch("nL",&s_nL,"s_nL/D");
	s_tree->Branch("aR",&s_aR,"s_aR/D");
	s_tree->Branch("nR",&s_nR,"s_nR/D");
	s_tree->Branch("mu_err",&s_mu_e,"s_mu_e/D");
	s_tree->Branch("sigma_err",&s_sigma_e,"s_sigma_e/D");
	s_tree->Branch("aL_err",&s_aL_e,"s_aL_e/D");
	s_tree->Branch("nL_err",&s_nL_e,"s_nL_e/D");
	s_tree->Branch("aR_err",&s_aR_e,"s_aR_e/D");
	s_tree->Branch("nR_err",&s_nR_e,"s_nR_e/D");
	s_tree->Branch("chi2",&s_chi2,"s_chi2/D");
	s_tree->Branch("imip",&s_imip,"s_imip/I");
	s_tree->Branch("ieta",&s_ieta,"s_ieta/I");
	s_tree->Branch("N",&s_N,"s_N/I");
	s_tree->Branch("ndf",&s_ndf,"s_ndf/I");

	//for storage of output info
	energyRes* res_temp;
	TF1* fit_temp;
	
	for(int ieta = 1; ieta <= maxHDeta; ieta++){ //loop over ieta
		for(int num = 0; num < maxHDe; num++) { //loop over energy
			for(int imip = 0; imip < 3; imip++){ //loop over imip
				if(ieta>30 && energies[num]<15) continue; //skip low energies in HF
				else if(ieta>31 && energies[num]<30 && imip==1) continue; //skip nomip in HF
				else if(ieta>32 && energies[num]<50 && imip==1) continue; //skip nomip in HF
				else if(ieta>=35 && imip==1) continue; //skip nomip in HF			
				res_temp = get_res(num,imip,ieta,1,0,0);
			
				//get values directly from energyRes
				s_energy = res_temp->getEnergy();
				s_ieta = res_temp->getEta();
				s_imip = res_temp->getMip();
				s_N = res_temp->getStat(0);
			
				//get fit parameters
				fit_temp = res_temp->getFit();
				s_mu = fit_temp->GetParameter(1);
				s_sigma = fit_temp->GetParameter(2);
				s_aL = fit_temp->GetParameter(3);
				s_nL = fit_temp->GetParameter(4);
				s_aR = fit_temp->GetParameter(5);
				s_nR = fit_temp->GetParameter(6);
				s_mu_e = fit_temp->GetParError(1);
				s_sigma_e = fit_temp->GetParError(2);
				s_aL_e = fit_temp->GetParError(3);
				s_nL_e = fit_temp->GetParError(4);
				s_aR_e = fit_temp->GetParError(5);
				s_nR_e = fit_temp->GetParError(6);
				s_chi2 = fit_temp->GetChisquare();
				s_ndf = fit_temp->GetNDF();
				
				//fill output
				s_tree->Fill();
			}
		}
	}
	
	//save output tree
	out_file->cd();
	s_tree->Write();
	out_file->Close();
}

//---------------------------
//function to print all plots
void alltheplots(){
	gROOT->SetBatch();
	
	//energy values
	Double_t energies[] = {1., 2., 3., 5., 9., 11., 15., 20., 30., 50., 100., 150., 225., 300., 1000., 3000. };	
	
	for(int ieta = 1; ieta <= maxHDeta; ieta++){ //loop over ieta
		for(int num = 0; num < maxHDe; num++) { //loop over energy
			for(int imip = 0; imip < 3; imip++){ //loop over imip
				get_res(num,imip,ieta,1,1,1,1);
			}
		}
	}

}

//----------------------------
//function to print python cfg
void fs_print_cfg(std::string oname="HcalResponse_new.py", std::string fname="tree_cballD_lim.root"){
	//energy values
	Double_t energies[] = {1., 2., 3., 5., 9., 11., 15., 20., 30., 50., 100., 150., 225., 300., 1000., 3000. };	

	//open file and tree
	TFile* file = TFile::Open(fname.c_str());
	TTree* tree = (TTree*)file->Get("tree");
	
	//open output file
	std::ofstream output(oname.c_str());
	if (!output) {
		std::cerr << "Cannot open the output file " << oname << "\n";
		return;
	}
	output << std::fixed << std::setprecision(5);
	
	int nPar = 6;
	std::string parNames[] = {"mu","sigma","aL","nL","aR","nR"};
	double defaults[] = {0,0,10,200,10,200};
	std::string dets[] = {"_HB","_HE","_HF"};
	int HDeta[] = {0,16,30,50}; //HCALResponse ieta = my ieta - 1
	int maxHDetas[] = {HDeta[1] - HDeta[0], HDeta[2] - HDeta[1], HDeta[3] - HDeta[2]};
	double minHDe[] = {1, 1, 15}; //minimum energies for each det
	std::string mips[] = {"_mip","_nomip",""};
	std::string s4 = "    "; //python indent string
	
	for(int p = 0; p < nPar; p++){ //loop over parameters
		for(int m = 0; m < 3; m++){ //loop over imip
			for(int d = 0; d < 3; d++){ //loop over dets
				output << s4 << s4 << (parNames[p] + dets[d] + mips[m]) << " = cms.vdouble( *[" << std::endl;
			
				for(int i = 0; i < maxHDe; i++){
					if(energies[i] < minHDe[d]) continue;
					
					//select par in eta range for this det, energy, and imip
					std::stringstream cutname;
					cutname << "imip==" << m << " && energy==" << energies[i] << " && ieta>" << HDeta[d] << " && ieta<=" << HDeta[d+1];
					Long64_t Npars = tree->Draw(parNames[p].c_str(),(cutname.str()).c_str(),"goff");
					Double_t *pars = tree->GetV1();
					
					output << s4 << s4 << s4;
					for(int j = 0; j < maxHDetas[d]; j++){
						//put customizations here
						if(j<Npars) {
							if((p==2 || p==4) && pars[j] < 0.00001) output << 0.00001; //min value for aL, aR
							else output << pars[j];
						}
						else output << defaults[p];
						if((j < maxHDetas[d]-1) || (i < maxHDe-1)) output << ", ";
					}
					
					if(i == maxHDe-1) output << " ] ";
					output << "#" << (int)(energies[i]) << std::endl;
				}
				
				output << s4 << s4 << ")," << std::endl;
			}
		}
	}
	
	


}

//--------------------------------------------
//macro to make ECAL-HCAL energy "banana plot"
void fs_banana_plot(int num, int ieta, bool do_print){
	gStyle->SetPalette(1);

	//energy values
	Double_t energies[] = {1., 2., 3., 5., 9., 11., 15., 20., 30., 50., 100., 150., 225., 300., 1000., 3000. };

	if (num>=maxHDe || num<0) {
		std::cout << "num must be between 0 and " << maxHDe - 1 << std::endl;
		return;
	}

	if (ieta>maxHDeta || ieta<1) {
		std::cout << "ieta must be between 1 and " << maxHDeta << std::endl;
		return;
	}

	//make filenames
	std::stringstream fname, piname;
	fname << "tree_pion_" << energies[num] << "_500k.root";
	piname << "#pi^{-} " << energies[num] << " GeV, " << (ieta-1)*0.1 << "#leq|#eta|<" << ieta*0.1;

	//open file and tree
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TTree* totalTree = (TTree*)_file->Get("Total");

	//make tree drawing expressions
	//define mip as ecal < 1 gev = 1000 mev
	double mipcut = 0.8;
	std::stringstream drawname[3];
	std::stringstream cutname[3];
	std::stringstream etacut;
	double hmin = 0.1*energies[num]; //lower cut to remove weird peaks near E=zero
	double hmax = 3*energies[num];
	double emax = 1.25*energies[num];

	//ecal & hcal energies are already "calibrated" at sim level (sampling factors for HBHEHO, poisson smearing & pe conversion for HF)
	drawname[0] << "hcal:ecal>>htemp0(100,0," << mipcut << ",100," << hmin << "," << hmax << ")";
	drawname[1] << "hcal:ecal>>htemp1(100,0," << emax << ",100," << hmin << "," << hmax << ")";
	drawname[2] << "hcal:ecal>>htemp2(100,0," << emax << ",100," << hmin << "," << hmax << ")";

	etacut << (ieta-1)*0.1 << " <= abs(eta) && abs(eta) < " << ieta*0.1;
	cutname[0] << "ecal < " << mipcut << " && " << etacut.str();
	cutname[1] << "ecal > " << mipcut << " && " << etacut.str();
	cutname[2] << etacut.str();

	TH2F* h_res[3]; //to store histos drawn from tree

	//plotting variables
	TCanvas* can[3];
	TPad* pad[3];
	TPaveText* pave[3];
	
	//loop over imip
	for(int i=0;i<3;i++){
		//draw w/ appropriate cut
		totalTree->Draw((drawname[i].str()).c_str(),(cutname[i].str()).c_str(),"hist goff");
		if(i==0) h_res[i] = (TH2F*)gDirectory->Get("htemp0");
		else if(i==1) h_res[i] = (TH2F*)gDirectory->Get("htemp1");
		else h_res[i] = (TH2F*)gDirectory->Get("htemp2");
		
		//names
		std::string omip;
		if (i==0) omip = "mip";
		else if (i==1) omip = "nomip";
		else omip = "tot";
		std::stringstream mipname;
		if(i==0) mipname << "mip, E_{ecal} < " << mipcut << " GeV";
		else if(i==1) mipname << "nomip, E_{ecal} > " << mipcut << " GeV";
		else mipname << "total";

		std::stringstream oname, Nname;
		oname << "pion_banana_plot_" << omip << "_" << energies[num] << "gev_ieta" << ieta << ".png";
		Nname << "N = " << h_res[i]->GetEntries();
		
		//plotting
		can[i] = new TCanvas(omip.c_str(),omip.c_str(),700,500);
		can[i]->cd();	
		pad[i] = new TPad("graph","",0,0,1,1);
		pad[i]->SetMargin(0.15,0.05,0.15,0.05);
		pad[i]->Draw();
		pad[i]->cd();
		
		//formatting
		h_res[i]->SetTitle("");
		h_res[i]->SetLineWidth(2);
		h_res[i]->GetYaxis()->SetTitleOffset(1.0);
		h_res[i]->GetYaxis()->SetTitleSize(32/(pad[i]->GetWh()*pad[i]->GetAbsHNDC()));
		h_res[i]->GetYaxis()->SetLabelSize(28/(pad[i]->GetWh()*pad[i]->GetAbsHNDC()));
		h_res[i]->GetXaxis()->SetTitleSize(32/(pad[i]->GetWh()*pad[i]->GetAbsHNDC()));
		h_res[i]->GetXaxis()->SetLabelSize(28/(pad[i]->GetWh()*pad[i]->GetAbsHNDC()));
		h_res[i]->GetYaxis()->SetTickLength(12/(pad[i]->GetWh()*pad[i]->GetAbsHNDC()));
		h_res[i]->GetXaxis()->SetTickLength(12/(pad[i]->GetWh()*pad[i]->GetAbsHNDC()));
			
		//plot histo
		h_res[i]->SetTitle("");
		h_res[i]->GetXaxis()->SetTitle("ECAL energy [GeV]");
		h_res[i]->GetYaxis()->SetTitle("HCAL energy [GeV]");
		h_res[i]->Draw("hist col");

		//pave
		pave[i] = new TPaveText(0.6,0.6,0.9,0.8,"NDC");
		pave[i]->AddText("CMSSW");
		pave[i]->AddText((piname.str()).c_str());
		pave[i]->AddText((mipname.str()).c_str());
		pave[i]->AddText((Nname.str()).c_str());
		pave[i]->SetFillColor(0);
		pave[i]->SetBorderSize(0);
		pave[i]->SetTextFont(42);
		pave[i]->SetTextSize(0.05);
		pave[i]->Draw("same");
		
		if(do_print){
			can[i]->Print((oname.str()).c_str(),"png");
		}
	}

}
