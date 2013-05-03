// Global FWCore classes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// user include files
#include "ForwardCaloUpgrade/FullSim/interface/FullSimEtaTimeAnalyzer.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

//SimHits
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

//Hcal det id
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

//Ecal det id
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

//cell geometry
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

//gen particles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;
using namespace edm;

double FullSimEtaTimeAnalyzer::phi(double x, double y) {
	double phi_ = atan2(y, x);
	return (phi_>=0) ?  phi_ : phi_ + 2*TMath::Pi();
}
double FullSimEtaTimeAnalyzer::DeltaPhi(double phi1, double phi2) {
	double phi1_= phi( cos(phi1), sin(phi1) );
	double phi2_= phi( cos(phi2), sin(phi2) );
	double dphi_= phi1_-phi2_;
	if( dphi_> TMath::Pi() ) dphi_-=2*TMath::Pi();
	if( dphi_<-TMath::Pi() ) dphi_+=2*TMath::Pi();

	return dphi_;
}
double FullSimEtaTimeAnalyzer::DeltaR(double phi1, double eta1, double phi2, double eta2){
	double dphi = DeltaPhi(phi1,phi2);
	double deta = eta2 - eta1;
	double dR2 = dphi*dphi + deta*deta;
	return sqrt(dR2);
}

FullSimEtaTimeAnalyzer::FullSimEtaTimeAnalyzer(const edm::ParameterSet& iConfig) { 
	outname = iConfig.getParameter<string>("fileName");
	egen = iConfig.getParameter<double>("eGen");
	dRcut = iConfig.getParameter<double>("dRcut");
	
}

FullSimEtaTimeAnalyzer::~FullSimEtaTimeAnalyzer() { }

// ------------ method called for each event  ------------
void
FullSimEtaTimeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	iSetup.get<CaloGeometryRecord>().get (geometry);

	double s_eta, s_phi;
	s_eta = s_phi = 0;
	
	//Access to GenParticles
	Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
	const reco::GenParticle & p = (*genParticles)[0];
	s_eta = p.eta();
	s_phi = p.phi();
		
	Handle<PCaloHitContainer> h_PCaloHitsH;
	bool bH = iEvent.getByLabel("g4SimHits","HcalHits", h_PCaloHitsH);
	
	//iterator
	PCaloHitContainer::const_iterator genSH;

	//get HCAL energy for this event
	if(bH){
		for(genSH = h_PCaloHitsH->begin(); genSH != h_PCaloHitsH->end(); genSH++) {
			HcalDetId cell(genSH->id());
			const CaloCellGeometry* cellGeometry = geometry->getSubdetectorGeometry(cell)->getGeometry(cell);
			double h_eta = cellGeometry->getPosition().eta();
			double h_phi = cellGeometry->getPosition().phi();
			double dR = DeltaR(s_phi,s_eta,h_phi,h_eta);
			HcalSubdetector det = cell.subdet();
		
			if(dR < dRcut){
				//cout << "det " << det << ", ieta " << cell.ietaAbs() << endl;
				if(det == HcalBarrel) h_Etime_eta[0][cell.ietaAbs()-1]->Fill(genSH->time(),genSH->energy());
				else if(det == HcalEndcap) h_Etime_eta[1][cell.ietaAbs()-maxHBeta]->Fill(genSH->time(),genSH->energy());
				else if(det == HcalForward) h_Etime_eta[3][cell.ietaAbs()-maxHEeta-maxHBeta+1]->Fill(genSH->time(),genSH->energy());
				else if(det == HcalOuter) h_Etime_eta[2][cell.ietaAbs()-1]->Fill(genSH->time(),genSH->energy());	
			}

		}
	}

}

// ------------ method called once each job just before starting event loop  ------------
void 
FullSimEtaTimeAnalyzer::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void 
FullSimEtaTimeAnalyzer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
FullSimEtaTimeAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
	out_file = new TFile(outname.c_str(), "RECREATE");
	
	string dname[maxHDdet] = {"HB","HE","HO","HF"};
	std::stringstream hname;
	int maxHDeta[maxHDdet] = {maxHBeta, maxHEeta, maxHOeta, maxHFeta};
	//book energy-weighted timing histos
	for(int d = 0; d < maxHDdet; d++){
		h_Etime_eta[d] = new TH1F*[maxHDeta[d]];
		for(int i = 0; i < maxHDeta[d]; i++){
			hname.str(std::string());
			hname.clear();
			hname << "Etime_" << dname[d] << "_ieta" << i;
			h_Etime_eta[d][i] = new TH1F((hname.str()).c_str(),(hname.str()).c_str(),500,0,100);
		}
	}
	
}

// ------------ method called when ending the processing of a run  ------------
void 
FullSimEtaTimeAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) { 
	out_file->cd();

	int maxHDeta[maxHDdet] = {maxHBeta, maxHEeta, maxHOeta, maxHFeta};
	//write energy-weighted timing histos
	for(int d = 0; d < maxHDdet; d++){
		for(int i = 0; i < maxHDeta[d]; i++){
			h_Etime_eta[d][i]->Write();
		}
	}
	
	out_file->Close();
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FullSimEtaTimeAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FullSimEtaTimeAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FullSimEtaTimeAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



