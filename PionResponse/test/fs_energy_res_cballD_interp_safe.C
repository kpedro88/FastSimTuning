#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
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
#define maxHDeta2 30 //eta bins for hadrons (no HF)
#define maxHDqty 6 //# qtys to plot

//qty: 0 - mu, 1 - sigma, 2 - aL, 3 - nL, 4 - aR, 5 - nR 
double interp(double E, int qty, double *energies, double *pars, int Npars){
	//find nearest energies to E
	int ie = -1;
	for(int i = 0; i < Npars; i++){
		if(E < energies[i]){
			if(i==0) ie = 0; // less than minimal - back extrapolation with the 1st interval
			else ie = i-1;
			break;
		}
	}
	if(ie==-1) ie = Npars - 2; // more than maximum - extrapolation with last interv.
	
	double x1 = energies[ie];
	double x2 = energies[ie+1];
	double y1 = pars[ie];
	double y2 = pars[ie+1];
	
	double y = (y1*(x2-E) + y2*(E-x1))/(x2-x1); //linear interpolation formula

	//do not let mu or sigma get extrapolated below zero for low energies
	//especially important for HF since extrapolation is used for E < 15 GeV
	if((qty==0 || qty==1) && E < x1){
		double tmp = (y1*x2-y2*x1)/(x2-x1); //extrapolate down to e=0
		if(tmp<0) { //require mu,sigma > 0 for E > 0
			y = y1*E/x1;
		}
	}
	//aL,nL,aR,nR have lower bounds - never extrapolate down
	else if((qty==2 || qty==3 || qty==4 || qty==5)){
		if(E < x1 && y1 < y2) y = y1;
		else if(E > x2 && y2 < y1) y = y2;
	}

	//final check for tail parameters
	//if((qty==2 || qty==4) && y<=0) y = 0.0001;
	//else if((qty==3 || qty==5) && y<=1) y = 1.0001;
	
	//return (y>0) ? y : 0;
	return y;
}

//imip: 0 = mip, 1 = nomip, 2 = total
//qty: 0 - response, 1 - resolution, 2 - aL, 3 - nL, 4 - aR, 5 - nR
void fs_plot_res(int qty, int imip, int ieta, bool do_print, std::string fname="tree_cballD_custom.root"){
	if (ieta>maxHDeta || ieta<1) {
		std::cout << "ieta must be between 1 and " << maxHDeta << std::endl;
		return;
	}
	if (qty>=maxHDqty || qty<0) {
		std::cout << "qty must be between 0 and " << maxHDqty-1 << std::endl;
		return;
	}
	if (imip>2 || imip<0) {
		std::cout << "imip must be between 0 and 2" << std::endl;
		return;
	}
	
	//open file and tree
	TFile* file = TFile::Open(fname.c_str());
	TTree* tree = (TTree*)file->Get("tree");

	//tree formulas for quantities and errors
	std::string qtys[] = {"mu/energy","sigma/mu","aL","nL","aR","nR"};
	std::string errs[] = {"mu_err/energy","sigma/mu*sqrt((sigma_err/sigma)^2+(mu_err/mu)^2)","aL_err","nL_err","aR_err","nR_err"};
	
	//get arrays from tree draw
	std::stringstream drawname, cutname;
	drawname << "energy:" << qtys[qty] << ":" << errs[qty];
	cutname << "imip==" << imip << " && ieta==" << ieta;
	Long64_t Npars = tree->Draw((drawname.str()).c_str(),(cutname.str()).c_str(),"goff");
	TGraphErrors* g_exact = new TGraphErrors(Npars,tree->GetV1(),tree->GetV2(),0,tree->GetV3());
	
	//get single quantities
	Double_t *energies, *pars1, *pars2;
	std::string drawname2[] = {"energy:mu","energy:sigma:mu","energy:aL","energy:nL","energy:aR","energy:nR"};
	Npars = tree->Draw(drawname2[qty].c_str(),(cutname.str()).c_str(),"goff");
	energies = tree->GetV1();
	pars1 = tree->GetV2();
	if(qty==1) pars2 = tree->GetV3(); //need mu for sigma
	
	Int_t Npts = 10000; //number of interpolated points to calculate
	Double_t step = 2*energies[Npars-1]/Npts; //large range, E=0..2*E_max
	Double_t *E = new Double_t[Npts];
	Double_t *p = new Double_t[Npts];
	for(int i = 0; i < Npts; i++){
		E[i] = (i+1)*step;
		if(qty==0) p[i] = interp(E[i],qty,energies,pars1,Npars)/E[i]; //mu/E
		else if(qty==1) {
			double m = interp(E[i],qty,energies,pars2,Npars);
			double s = interp(E[i],qty,energies,pars1,Npars);
			if(m>0) p[i] = s/m;//sigma/mu
			else p[i] = 0;
		}
		else p[i] = interp(E[i],qty,energies,pars1,Npars);
	}
	TGraphErrors* g_interp = new TGraphErrors(Npts,E,p,0,0);

	TGraphErrors** graphs = new TGraphErrors*[2];
	graphs[0] = g_exact;
	graphs[1] = g_interp;

	std::string yname[maxHDqty] = {"Response (#mu/E_{true})","Resolution (#sigma/#mu)","a_{L}","n_{L}","a_{R}","n_{R}"};
	std::stringstream etaname;
	etaname << (ieta-1)*0.1 << "#leq|#eta|<" << ieta*0.1;
	Int_t col = kBlack;
	Int_t mrk[2] = {20, 33};

	TCanvas* can;
	TPaveText* pave;

	//formatting
	for(int i = 0; i < 2; i++){
		graphs[i]->GetXaxis()->SetTitle("Energy [GeV]");
		graphs[i]->GetYaxis()->SetTitle(yname[qty].c_str());
		graphs[i]->SetTitle("");
		graphs[i]->SetMarkerStyle(mrk[i]);
		graphs[i]->SetMarkerColor(col);
		graphs[i]->SetMarkerSize(1.5);
		graphs[i]->SetLineColor(col);
		graphs[i]->SetFillColor(0);
	}
	//if(qty==2 || qty==4) graph->GetYaxis()->SetRangeUser(0,10);
	//else if(qty==3 || qty==5) graph->GetYaxis()->SetRangeUser(1,200);

	std::string cname[maxHDqty] = {"mu","sigma","aL","nL","aR","nR"};
	can = new TCanvas(cname[qty].c_str(),cname[qty].c_str(),700,500);
	can->cd();
	can->SetLogx();
	if(qty==3 || qty==5) can->SetLogy();

	graphs[1]->Draw("AL");
	graphs[0]->Draw("PZ same");

	//legend, pave coords
	double y1;
	if(qty==1) y1 = 0.5;
	else y1 = 0.2;	
	
	std::string gname = "CMSSW";
	
	pave = new TPaveText(0.8,y1,0.95,y1+0.15,"NDC");
	pave->AddText(gname.c_str());
	pave->AddText((etaname.str()).c_str());
	if(imip==0) pave->AddText("mip, E_{ecal} < 1 GeV");
	else if(imip==1) pave->AddText("nomip, E_{ecal} > 1 GeV");
	else pave->AddText("total");
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	pave->Draw("same");
	
	if(do_print){

	}
}

void fs_plot_qtys(int imip, int ieta, bool do_print, std::string fname="tree_cballD_custom.root"){
	for(int qty = 0; qty < maxHDqty; qty++){
		fs_plot_res(qty,imip,ieta,do_print,fname);
	}

}

void fs_all_res(int qty, int imip, bool do_print, std::string fname="tree_cballD_custom.root"){
	Double_t energies[] = {0.,1., 2., 3., 5., 9., 11., 15., 20., 30., 50., 100., 150., 225., 300., 1000., 3000.};


	gStyle->SetPalette(1); //rainbow
	//gStyle->SetPalette(54); //blue-yellow

	//open files and trees
	TFile* file;
	file = TFile::Open(fname.c_str());
	TTree* tree;
	tree = (TTree*)file->Get("tree");
	
	//result histo
	Int_t Npts = 50000; //number of interpolated energy points to calculate
	Double_t Emin = 0.1; //minimum energy for extrapolation
	Double_t step = 2*(energies[maxHDe-1]-Emin)/Npts; //large range, E=Emin..2*E_max
	TH2D* graph2d = new TH2D("graph2d","",maxHDeta,0.5,maxHDeta+0.5,Npts,0.5*step,2*energies[maxHDe-1]+0.5*step);
	graph2d->SetContour(1000);
	
	//loop over eta
	for(int ieta=1; ieta<=maxHDeta; ieta++){
		//get single quantities
		Double_t *ens, *pars1, *pars2;
		std::stringstream cutname;
		cutname << "imip==" << imip << " && ieta==" << ieta;
		std::string drawname[] = {"energy:mu","energy:sigma:mu","energy:aL","energy:nL","energy:aR","energy:nR"};
		Long64_t Npars = tree->Draw(drawname[qty].c_str(),(cutname.str()).c_str(),"goff");
		ens = tree->GetV1();
		pars1 = tree->GetV2();
		if(qty==1) pars2 = tree->GetV3(); //need mu for sigma
		
		//do interpolations and fill histo
		for(int i = 0; i < Npts; i++){
			double E = Emin+i*step;
			double p = 0;
			if(qty==0) p = interp(E,qty,ens,pars1,Npars)/E; //mu/E
			else if(qty==1) {
				double m = interp(E,qty,ens,pars2,Npars);
				double s = interp(E,qty,ens,pars1,Npars);
				if(m>0) p = s/m;//sigma/mu
				else p = 0;
			}
			else p = interp(E,qty,ens,pars1,Npars);
			
			int xbin = graph2d->GetXaxis()->FindBin(ieta);
			int ybin = graph2d->GetYaxis()->FindBin(E);
			graph2d->SetBinContent(xbin,ybin,p);
		}
	}
	//fix scale
	double zmax[] = {1.1,1.5,10,200,10,200};
	graph2d->GetZaxis()->SetRangeUser(0,zmax[qty]);
	graph2d->GetYaxis()->SetRangeUser(Emin,2*energies[maxHDe-1]);
	
	//format graphs
	graph2d->GetYaxis()->SetTitle("Energy [GeV]");
	graph2d->GetXaxis()->SetTitle("i#eta");
	graph2d->GetXaxis()->SetTitleOffset(0.9);
	graph2d->GetYaxis()->SetTitleOffset(0.9);
	graph2d->SetTitle("");

	//names
	std::string omip;
	if (imip==0) omip = "mip";
	else if (imip==1) omip = "nomip";
	else omip = "tot";
	std::string qname[maxHDqty] = {"mu","sigma","aL","nL","aR","nR"};
	
	//draw
	TCanvas* can = new TCanvas((qname[qty]+"_2D").c_str(),(qname[qty]+"_2D").c_str(),700,500);
	TPad* pad = new TPad("graph","",0,0,1,1);
	pad->SetMargin(0.12,0.15,0.12,0.08);
	pad->Draw();
	pad->cd();
	pad->SetLogy();
	//graph->Draw("pcol");
	graph2d->Draw("hist colz");

	//pave
	std::string yname[maxHDqty] = {"Response (#mu/E_{true})","Resolution (#sigma/#mu)","a_{L}","n_{L}","a_{R}","n_{R}"};
	std::string zname = "z = " + yname[qty] + " (" + omip + ")";
	TPaveText* pave = new TPaveText(0.1,0.93,0.5,1.0,"NDC");
	pave->AddText(zname.c_str());
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	pave->Draw("same");
	
	if(do_print){
		std::string oname = qname[qty] + "_" + omip + "_interp.png";
		can->Print(oname.c_str(),"png");
	}

}

void fs_all_qtys(int imip, bool do_print, std::string fname="tree_cballD_custom.root"){
	for(int qty = 0; qty < maxHDqty; qty++){
		fs_all_res(qty,imip,do_print,fname);
	}

}