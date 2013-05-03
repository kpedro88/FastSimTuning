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
#define maxHDqty 7 //# qtys to plot

//imip: 0 = mip, 1 = nomip, 2 = total
//qty: 0 - response, 1 - resolution, 2 - aL, 3 - nL, 4 - aR, 5 - nR, 6 - chi2/ndf
TGraphErrors* fs_plot_res(int qty, int imip, int ieta, bool do_print, std::string fname="tree_cballD_custom.root"){
	if (ieta>maxHDeta || ieta<1) {
		std::cout << "ieta must be between 1 and " << maxHDeta << std::endl;
		return 0;
	}
	if (qty>=maxHDqty || qty<0) {
		std::cout << "qty must be between 0 and " << maxHDqty-1 << std::endl;
		return 0;
	}
	if (imip>2 || imip<0) {
		std::cout << "imip must be between 0 and 2" << std::endl;
		return 0;
	}
	
	//open file and tree
	TFile* file = TFile::Open(fname.c_str());
	TTree* tree = (TTree*)file->Get("tree");

	//tree formulas for quantities and errors
	std::string qtys[] = {"mu/energy","sigma/mu","aL","nL","aR","nR","chi2/ndf"};
	std::string errs[] = {"mu_err/energy","sigma/mu*sqrt((sigma_err/sigma)^2+(mu_err/mu)^2)","aL_err","nL_err","aR_err","nR_err","0"};
	
	//get arrays from tree draw
	std::stringstream drawname, cutname;
	drawname << "energy:" << qtys[qty] << ":" << errs[qty];
	cutname << "imip==" << imip << " && ieta==" << ieta;
	Long64_t Npts = tree->Draw((drawname.str()).c_str(),(cutname.str()).c_str(),"goff");
	
	//put together graph
	TGraphErrors* graph = new TGraphErrors(Npts,tree->GetV1(),tree->GetV2(),0,tree->GetV3());

	std::string yname[maxHDqty] = {"Response (#mu/E_{true})","Resolution (#sigma/#mu)","#alpha_{L}","n_{L}","#alpha_{R}","n_{R}","#chi^{2}/ndf"};
	std::stringstream etaname;
	etaname << (ieta-1)*0.1 << "#leq|#eta|<" << ieta*0.1;
	Int_t col, mrk;
	col = kBlack; mrk = 20;

	TCanvas* can;
	TPaveText* pave;

	graph->GetXaxis()->SetTitle("Energy [GeV]");
	graph->GetYaxis()->SetTitle(yname[qty].c_str());
	graph->SetTitle("");
	graph->SetMarkerStyle(mrk);
	graph->SetMarkerColor(col);
	graph->SetMarkerSize(1.5);
	graph->SetLineColor(col);
	graph->SetFillColor(0);
	if(qty==2 || qty==4) graph->GetYaxis()->SetRangeUser(0,10);
	else if(qty==3 || qty==5) graph->GetYaxis()->SetRangeUser(1,200);


	std::string cname[maxHDqty] = {"mu","sigma","aL","nL","aR","nR","chi2"};
	can = new TCanvas(cname[qty].c_str(),cname[qty].c_str(),700,500);
	can->cd();
	can->SetLogx();
	if(qty==3 || qty==5) can->SetLogy();

	graph->Draw("APZ");

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

	return graph;
}

void fs_plot_qtys(int imip, int ieta, bool do_print, std::string fname="tree_cballD_custom.root"){
	for(int qty = 0; qty < maxHDqty; qty++){
		fs_plot_res(qty,imip,ieta,do_print,fname);
	}

}

void fs_all_res(int qty, int imip, bool do_print, std::string fname="tree_cballD_custom.root"){
	gStyle->SetPalette(1); //rainbow
	//gStyle->SetPalette(54); //blue-yellow

	//open files and trees
	TFile* file;
	file = TFile::Open(fname.c_str());
	TTree* tree;
	tree = (TTree*)file->Get("tree");
	
	//tree variables
	int s_imip;
	int s_ieta;
	int s_N, s_ndf;
	double s_energy;
	double s_mu, s_sigma, s_aL, s_nL, s_aR, s_nR, s_chi2;
	int nentries;
	
	//set tree branches
	tree->SetBranchAddress("imip",&(s_imip));
	tree->SetBranchAddress("ieta",&(s_ieta));
	tree->SetBranchAddress("N",&(s_N));
	tree->SetBranchAddress("energy",&(s_energy));
	tree->SetBranchAddress("mu",&(s_mu));
	tree->SetBranchAddress("sigma",&(s_sigma));
	tree->SetBranchAddress("aL",&(s_aL));
	tree->SetBranchAddress("nL",&(s_nL));
	tree->SetBranchAddress("aR",&(s_aR));
	tree->SetBranchAddress("nR",&(s_nR));
	tree->SetBranchAddress("chi2",&(s_chi2));
	tree->SetBranchAddress("ndf",&(s_ndf));
	nentries = tree->GetEntries();

	//result histo
	Double_t energies[] = {0.,1., 2., 3., 5., 9., 11., 15., 20., 30., 50., 100., 150., 225., 300., 1000., 3000., 5000. };
	Double_t *ebins = new Double_t[maxHDe+1];
	for(int i = 0; i < maxHDe+1; i++){
		ebins[i] = energies[i] + (energies[i+1]-energies[i])/2;
	}
	TH2D* graph2d = new TH2D("graph2d","",maxHDeta,0.5,maxHDeta+0.5,maxHDe,ebins);
	graph2d->SetContour(50);
	
	//loop over trees
	for(int j = 0; j < nentries; j++){
		tree->GetEntry(j);
		
		//skip events with wrong imip
		if(s_imip != imip) continue;
		
		//skip low-energy events in HF
		if(s_ieta>30 && s_energy<20) continue;
		
		//skip non-crossover events in HF for nomip
		if(imip == 1 && s_ieta>=33) continue;
	
		int xbin = graph2d->GetXaxis()->FindBin(s_ieta);
		int ybin = graph2d->GetYaxis()->FindBin(s_energy);
		double val = 0;
		if(qty==0 && s_mu>0) val = s_mu/s_energy;
		else if(qty==1 && s_mu>1) val = s_sigma/s_mu;
		else if(qty==2) val = s_aL;
		else if(qty==3) val = s_nL;
		else if(qty==4) val = s_aR;
		else if(qty==5) val = s_nR;
		else if(qty==6) val = s_chi2/s_ndf;//chi2/ndf

		graph2d->SetBinContent(xbin,ybin,val);
	}
	
	//fix scale
	//Double_t zmax = TMath::Abs(graph2d->GetBinContent(graph2d->GetMaximumBin()));
	//Double_t zmin = TMath::Abs(graph2d->GetBinContent(graph2d->GetMinimumBin()));
	//zmax = (zmax>zmin) ? zmax : zmin;
	//graph2d->GetZaxis()->SetRangeUser(-zmax,zmax);
	//graph2d->GetZaxis()->SetRangeUser(-5,5);
	if(qty==6) graph2d->GetZaxis()->SetRangeUser(0,10);
	
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
	std::string qname[maxHDqty] = {"mu","sigma","aL","nL","aR","nR","chi2"};
	
	//draw
	TCanvas* can = new TCanvas(qname[qty].c_str(),qname[qty].c_str(),700,500);
	TPad* pad = new TPad("graph","",0,0,1,1);
	pad->SetMargin(0.12,0.15,0.12,0.08);
	pad->Draw();
	pad->cd();
	pad->SetLogy();
	//graph->Draw("pcol");
	graph2d->Draw("hist colz");

	//pave
	std::string yname[maxHDqty] = {"Response (#mu/E_{true})","Resolution (#sigma/#mu)","#alpha_{L}","n_{L}","#alpha_{R}","n_{R}","#chi^{2}/ndf"};
	std::string zname = "z = " + yname[qty] + " (" + omip + ")";
	TPaveText* pave = new TPaveText(0.1,0.93,0.5,1.0,"NDC");
	pave->AddText(zname.c_str());
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	pave->Draw("same");
	
	if(do_print){
		std::string oname = qname[qty] + "_" + omip + ".png";
		can->Print(oname.c_str(),"png");
	}

}

void fs_all_qtys(int imip, bool do_print, std::string fname="tree_cballD_custom.root"){
	for(int qty = 0; qty < maxHDqty; qty++){
		fs_all_res(qty,imip,do_print,fname);
	}

}