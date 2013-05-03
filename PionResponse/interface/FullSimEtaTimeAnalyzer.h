#ifndef FullSimEtaTimeAnalyzer_h
#define FullSimEtaTimeAnalyzer_h

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//SimHits
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

//geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

//random engine
#include "CLHEP/Random/RandPoissonQ.h"

#define maxHBeta 16
#define maxHOeta 15
#define maxHEeta 14
#define maxHFeta 21
#define maxHDdet 4

class TFile;
class TTree;
class TH1F;

class FullSimEtaTimeAnalyzer : public edm::EDAnalyzer {
	public:
		explicit FullSimEtaTimeAnalyzer(const edm::ParameterSet&);
		~FullSimEtaTimeAnalyzer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		
	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
	
		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

		double phi(double x, double y);
		double DeltaPhi(double phi1, double phi2);
		double DeltaR(double phi1, double eta1, double phi2, double eta2);		
		
		//member variables
		TFile* out_file;
		TH1F **h_Etime_eta[maxHDdet];
		std::string outname;
		double egen;
		double dRcut;
		edm::ESHandle<CaloGeometry> geometry;
};

//define this as a plug-in
DEFINE_FWK_MODULE(FullSimEtaTimeAnalyzer);

#endif