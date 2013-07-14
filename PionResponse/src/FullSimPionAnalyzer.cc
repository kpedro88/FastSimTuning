// Global FWCore classes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// user include files
#include "FastSimTuning/PionResponse/interface/FullSimPionAnalyzer.h"
#include <cmath>
#include <iostream>

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

//random engine
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandPoissonQ.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace edm;

double FullSimPionAnalyzer::phi(double x, double y) {
	double phi_ = atan2(y, x);
	return (phi_>=0) ?  phi_ : phi_ + 2*TMath::Pi();
}
double FullSimPionAnalyzer::DeltaPhi(double phi1, double phi2) {
	double phi1_= phi( cos(phi1), sin(phi1) );
	double phi2_= phi( cos(phi2), sin(phi2) );
	double dphi_= phi1_-phi2_;
	if( dphi_> TMath::Pi() ) dphi_-=2*TMath::Pi();
	if( dphi_<-TMath::Pi() ) dphi_+=2*TMath::Pi();

	return dphi_;
}
double FullSimPionAnalyzer::DeltaR(double phi1, double eta1, double phi2, double eta2){
	double dphi = DeltaPhi(phi1,phi2);
	double deta = eta2 - eta1;
	double dR2 = dphi*dphi + deta*deta;
	return sqrt(dR2);
}

FullSimPionAnalyzer::FullSimPionAnalyzer(const edm::ParameterSet& iConfig) { 
	outname = iConfig.getParameter<string>("fileName");
	
    samplingHBHE = iConfig.getParameter<vector<double> >("samplingHBHE");
    samplingHF = iConfig.getParameter<vector<double> >("samplingHF");
	simHitToPhotoelectrons = iConfig.getParameter<vector<double> >("simHitToPhotoelectrons");
    samplingHO = iConfig.getParameter<vector<double> >("samplingHO");
	
	dRcut = iConfig.getParameter<double>("dRcut");
	
	// Initialize the random number generator service
	edm::Service<edm::RandomNumberGenerator> rng;
	if ( ! rng.isAvailable() ) {
	throw cms::Exception("Configuration")
		<< "FullSimPionAnalyzer requires the RandomGeneratorService\n"
		"which is not present in the configuration file.\n"
		"You must add the service in the configuration file\n"
		"or remove the module that requires it";
	}
	//initialize poisson distribution
	poisson = new CLHEP::RandPoissonQ(rng->getEngine(), 1);
	
}

FullSimPionAnalyzer::~FullSimPionAnalyzer() { }

// ------------ method called for each event  ------------
void
FullSimPionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	iSetup.get<CaloGeometryRecord>().get (geometry);

	double sum_ecal, sum_hcal;
	sum_ecal = sum_hcal = 0;
	double s_eta, s_phi;
	s_eta = s_phi = 0;
	
	//Access to GenParticles
	Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
	const reco::GenParticle & p = (*genParticles)[0];
	s_eta = p.eta();
	s_phi = p.phi();

	//Access to simHits information
	Handle<PCaloHitContainer> h_PCaloHitsEE;
	bool bEE = iEvent.getByLabel("g4SimHits","EcalHitsEE", h_PCaloHitsEE);
	Handle<PCaloHitContainer> h_PCaloHitsEB;
	bool bEB = iEvent.getByLabel("g4SimHits","EcalHitsEB", h_PCaloHitsEB);
		
	Handle<PCaloHitContainer> h_PCaloHitsH;
	bool bH = iEvent.getByLabel("g4SimHits","HcalHits", h_PCaloHitsH);
	
	//iterator
	PCaloHitContainer::const_iterator genSH;

	//get ECAL energy for this event
	if(bEE){
		for(genSH = h_PCaloHitsEE->begin(); genSH != h_PCaloHitsEE->end(); genSH++) {
			EEDetId cell(genSH->id());
			const CaloCellGeometry* cellGeometry = geometry->getSubdetectorGeometry(cell)->getGeometry(cell);
			double h_eta = cellGeometry->getPosition().eta();
			double h_phi = cellGeometry->getPosition().phi();
			double dR = DeltaR(s_phi,s_eta,h_phi,h_eta);
			
			if(dR < dRcut) sum_ecal += genSH->energy();
		}
	}
	if(bEB){
		for(genSH = h_PCaloHitsEB->begin(); genSH != h_PCaloHitsEB->end(); genSH++) {
			EBDetId cell(genSH->id());
			const CaloCellGeometry* cellGeometry = geometry->getSubdetectorGeometry(cell)->getGeometry(cell);
			double h_eta = cellGeometry->getPosition().eta();
			double h_phi = cellGeometry->getPosition().phi();
			double dR = DeltaR(s_phi,s_eta,h_phi,h_eta);
			
			if(dR < dRcut) sum_ecal += genSH->energy();
		}
	}

	//get HCAL energy for this event
	if(bH){
		for(genSH = h_PCaloHitsH->begin(); genSH != h_PCaloHitsH->end(); genSH++) {
			HcalDetId cell(genSH->id());
			const CaloCellGeometry* cellGeometry = geometry->getSubdetectorGeometry(cell)->getGeometry(cell);
			double h_eta = cellGeometry->getPosition().eta();
			double h_phi = cellGeometry->getPosition().phi();
			double dR = DeltaR(s_phi,s_eta,h_phi,h_eta);
			
			if(dR < dRcut) {
				HcalSubdetector det = cell.subdet();
				if(det == HcalBarrel || det == HcalEndcap) sum_hcal += samplingHBHE[cell.ietaAbs()-1]*(genSH->energy());
				else if(det == HcalForward) {
					//HF is difficult
					double depth = cell.depth();
					double mean_pe = genSH->energy();
					double smeared_pe = poisson->fire(mean_pe*simHitToPhotoelectrons[depth-1])/simHitToPhotoelectrons[depth-1];
					sum_hcal += smeared_pe/samplingHF[depth-1];
				}
				else if(det == HcalOuter) sum_hcal += samplingHO[cell.ietaAbs()-1]*(genSH->energy());
			}
		}
	}
	
	e_ecal = sum_ecal;
	e_hcal = sum_hcal;
	gen_eta = s_eta;
	gen_phi = s_phi;
	
	tree_tot->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
FullSimPionAnalyzer::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void 
FullSimPionAnalyzer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
FullSimPionAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
	out_file = new TFile(outname.c_str(), "RECREATE");
	tree_tot = new TTree("Total", "Energy Calorimeter info");
	tree_tot->Branch("ecal",&e_ecal,"e_ecal/D");
	tree_tot->Branch("hcal",&e_hcal,"e_hcal/D");
	tree_tot->Branch("eta",&gen_eta,"gen_eta/D");
	tree_tot->Branch("phi",&gen_phi,"gen_phi/D");
}

// ------------ method called when ending the processing of a run  ------------
void 
FullSimPionAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) { 
	out_file->cd();
	
	tree_tot->Write();
	
	out_file->Close();
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FullSimPionAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FullSimPionAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FullSimPionAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



