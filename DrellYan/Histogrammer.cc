/*******************************************************************************************
 * Histogrammer.cc
 *
 * creates a root file containing histograms
 *
 *
 *******************************************************************************************
**/


#include <ctime>
#include <sstream>

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEfficiency.h"

#include "TRandom.h"
#include "TRandom3.h"

#include <TCanvas.h>
#include <TGraph.h>


//#include "TreeParticles.hpp"
#include "../../cms_sw/CMSSW_7_4_12/src/TreeWriter/TreeWriter/plugins/TreeParticles.hpp"

/*******************************************************************************************
 * TreeParticles.hpp
 *
 *******************************************************************************************

#ifndef TREEPARTICLES_H
#define TREEPARTICLES_H

#include <TLorentzVector.h>
#include <TVector3.h>

namespace tree
{
	struct Particle
	{
		TVector3 p;
		Float_t  someTestFloat=0.; // only for quick peeking
	};

	struct GenParticle: public Particle
	{
		Int_t pdgId=0;
	};

	struct Photon : public Particle
	{
		Float_t sigmaIetaIeta; // full 5x5
		Float_t hOverE;
		Int_t hasPixelSeed;
		Int_t passElectronVeto;
		Float_t r9;

		Float_t isoChargedHadronsEA;
		Float_t isoNeutralHadronsEA;
		Float_t isoPhotonsEA;
		Float_t isoWorstChargedHadrons;

		Int_t isTrue;
		Int_t isTrueAlternative;

		// IDs
		Bool_t  isLoose;
		Bool_t  isMedium;
		Bool_t  isTight;
		Float_t mvaValue;
	};

	struct Jet : public Particle
	{
		bool isLoose;
		float bDiscriminator;
	};

	struct Muon: public Particle
	{
		bool isTight;
	};

	struct Electron: public Particle
	{
		bool isLoose;
		bool isMedium;
		bool isTight;
	};

	struct MET : public Particle
	{
		TVector3 p_raw;
		Float_t  uncertainty;
	};

	inline bool EtGreater(const tree::Particle p1, const tree::Particle p2) {
		return p1.p.Pt() > p2.p.Pt();
	}

} // end namespace definition
#endif // TREEPARTICLES_H 

*******************************************************************************************/

// */
// namespace stuff
using namespace tree;
using namespace std;

/*******************************************************************************************
 * 
 * own functions and stuff
 *
 ************************************************/

//TCanvas* results;

// template to calculate the sign
template <typename T> int sign(T val) {
      return (T(0) < val) - (val < T(0));
}

// convert numer to string
template <typename T> string NumberToString ( T Number ) {
      ostringstream ss;
      ss << Number;
      return ss.str();
  }


// create output file name depending on input tree
string getOutputFilename( string strIn ) {

  // Converts "/path/to/ntuple/QCD_nTuple.root" to "QCD_hists.root"

  auto PosStart = strIn.rfind("/");
  auto PosEnd = strIn.find(".root");
  string outputName;
  
  if( PosEnd != string::npos ) {
    outputName = "../selector_rootFiles/my_selector_results_" + 
    				strIn.substr( PosStart+1, PosEnd-PosStart-1 ) + 
    				"_" + 
    				NumberToString(time(0)) +
    				".root";
  }
  return outputName;

}




/*
 * 
 * function to count Ng and Ne per bin for f=Ng/(Ne+Ng)
 * to calculate fake rate dependency on variables
 * 
 ************************************************/ 

/*
	void Histogrammer::Count_Pt(float Pt, float Ng_, float Ne_){
		
		for(int i=0; i<20;i++){
			if(Pt>=(i*10) && Pt<(i+1)*10){
				
				Ng["Pt"].at(i) += Ng_;
				Ne["Pt"].at(i) += Ne_;
				
			}
		}
		
	}
// */



/*******************************************************************************************
 * 
 * Class Histogrammer
 *
 ************************************************/
class Histogrammer : public TSelector {
	public:
		
		Histogrammer(); // constructor defined below
		
		virtual ~Histogrammer() { }

		virtual void    Init(TTree *tree_);
		virtual void	Begin(TTree *tree_);
		virtual void    SlaveBegin(TTree *tree_);
		virtual Bool_t  Process(Long64_t entry);
		virtual void	SlaveTerminate();
		virtual void    Terminate();
		virtual Int_t   Version() const { return 2; }

		void resetSelection();
		void defaultSelection();

		void initSelection( string const& s );
		void fillSelection( string const& s );

		void initObjects( string const& s );
		void fillObjects( string const& s );

		void initBkgEst( string const& s );
		void fillBkgEst( string const& s );
		
		
		
		

		TTreeReader fReader;
		TTreeReaderArray<tree::Photon> photons;
		TTreeReaderArray<tree::Jet> jets;
		TTreeReaderArray<tree::Electron> electrons;
		TTreeReaderArray<tree::Muon> muons;
		TTreeReaderArray<tree::Particle> genJets;
		TTreeReaderArray<tree::GenParticle> genParticles;

		TTreeReaderValue<tree::MET> met;
		TTreeReaderValue<Float_t> w_pu;
		TTreeReaderValue<Float_t> w_mc;
		
		TTreeReaderValue<Float_t> genHt;

		TTreeReaderValue<Int_t> nGoodVertices;
		TTreeReaderValue<Int_t> nChargedPfCandidates;
		TTreeReaderValue<Int_t> nPV;
		TTreeReaderValue<Int_t> genLeptonsFromW;

		TTreeReaderValue<ULong64_t> eventNo;
		TTreeReaderValue<UInt_t> runNo;
		TTreeReaderValue<UInt_t> lumNo;

		//TTreeReaderValue<UInt_t> size;

		TTreeReaderValue<Bool_t> isRealData;
		TTreeReaderValue<Bool_t> signalTriggerFired;
		TTreeReaderValue<Bool_t> crossTriggerPhotonFired;
		TTreeReaderValue<Bool_t> crossTriggerHtFired;

		vector<tree::Photon*> selPhotons;
		vector<tree::Jet*> selJets;
		vector<tree::Jet*> selBJets;
		vector<tree::Electron*> selElectrons;
		vector<tree::Muon*> selMuons;
		
		// genParticles

		vector<tree::Photon*> selPhotons2;
		vector<tree::Jet*> selJets2;
		vector<tree::Jet*> selBJets2;
		vector<tree::Electron*> selElectrons2;
		vector<tree::Muon*> selMuons2;
		//vector<tree::

		// INIT
		// 
		// /************************************************************************
		
		float m_Z = 91.1876;	// Z mass in GeV

		Float_t selW = 1.;	// weight
 		float selHt=0;
		int Nnum = 0;
		int	Ncut1 = 0,
			Ncut2 = 0,
			Ncut3 = 0,
			Ncut4 = 0,
			Ntest1 = 0,
			Ntest2 = 0,
			Ntest3 = 0,
			Ntest4 = 0;
			

		float Ht_2g = 0;

		int Ntwo = 0;
		
		// some variables
		float ht=0;
		
		vector<string> vs;

		float	sump_gg=0.,
			sump_eg=0.;
		float	pt1,
			pt2,
			eta1,
			eta2,
			phi1,
			phi2;
		float 	Msquare=0.,
			Htjets=0.;
		float 	Etemp1=0.,
			Etemp2=0.;
			

		struct S_mZdiff_m{
			float diff;
			float m;
		};
		
		vector<S_mZdiff_m>	ownVec1,
					ownVec2,
					ownVec3,
					ownVec4;
		
		int numE_tot = 0;		

		bool 	oneHasNoPixelSeed = false;	// set this TRUE if one has NO pixelseed and
							// therefore seems to be a "real" photon
		
		TVector3 	vTemp1,			// temp vector
				vTemp2,			// temp vector
				vTemp3,			//
				vTemp4;			//
		
		TLorentzVector	lvTemp1,		// temp lorentzvector
				lvTemp2,		// temp lorentzvector
				lvTemp3,		//
				lvTemp4;		//

		vector<TLorentzVector> 	vec_lvTemp1,	// vector to save lorentzvectors
					vec_lvTemp2,
					vec_lvTemp3,
					vec_lvTemp4;	//
					
		vector<float>	vecTemp1,
				vecTemp2,
				vecTemp3,
				vecTemp4;			

		vector<float> 	BestZMatch_ElectronPhoton_Denom,	//
				BestZMatch_ElectronPhoton_Num;
		vector<float> M_ElectronPhoton;			//
		
		float 	Ng_1 =0,
			Ne_1 =0;
			
		float	Ne_pt[10],
			Ng_pt[10];
			
		// for dependency on variables
		map<string,vector<float>> Ne, Ng, f, Xval;
		
		//struct
		//map<string,float> Ne, Ng;
		
		// ************************************************************************
		// */
		
		//pair<float,float> mrr2;
		
		map<string,TEfficiency> eff;
		map<string,TH1F> h;
		map<string,TGraph> g;
		
		map<string,TH2F> h2;
		
		//TLorentzVector v;

		//map<string,TLorentzVector> v; // lorentz vectors
		
		/*map<string,TH2F> h2;
		map<string,TProfile> prof;
 		map<string,TEfficiency> eff;
		// */
};

/*******************************************************************************************
 * 
 * constructor: link the branches to the access variables
 *
 ************************************************/
Histogrammer::Histogrammer(): 
	photons( fReader, "photons" ),
	jets( fReader, "jets" ),
	electrons( fReader, "electrons" ),
	muons( fReader, "muons" ),
	genJets( fReader, "genJets" ),
	genParticles( fReader, "genParticles" ),
	met( fReader, "met" ),
	nGoodVertices( fReader, "nGoodVertices" ),
	//nChargedPfCandidates( fReader, "nChargedPfCandidates" ),
	w_pu( fReader, "pu_weight" ),
	
	genHt( fReader, "genHt"),
	//w_mc( fReader, "mc_weight" ),
	
	//isRealData( fReader, "isRealData" ),
	signalTriggerFired( fReader, "HLT_Photon90_CaloIdL_PFHT500_v" ),
	crossTriggerPhotonFired( fReader, "HLT_Photon90_v" ),
	crossTriggerHtFired( fReader, "HLT_PFHT600_v" )
{
	cout << "constructor()" << endl;
} // constructor


/*******************************************************************************************
 * 
 * initialize Selection
 *
 ************************************************/
void Histogrammer::initSelection( string const& s ) {
	
	cout << "initSelection()" << endl;	
	
	eff["Pt"] = TEfficiency(("eff_Pt_"+s).c_str(),"my efficiency;p;fakerate",200,0,200);
	//eff["Pt"] = ->SetDirectory(gDirectory);
	
	// TH1F(const char* name, const char* title, nbins, binmin, binmax)
	h["m_gg"] 	= TH1F(("m_gg_"+s).c_str(),";m_{#gamma#gamma};a.u.", 100, 0, 200);
	h["m_eg"] 	= TH1F(("m_eg_"+s).c_str(),";m_{e#gamma};a.u.", 100, 0, 200);
	h["m_gg_all"]	= TH1F(("m_gg_all_"+s).c_str(),";m_{#gamma#gamma};a.u.", 100, 0, 200);
	h["m_eg_all"]	= TH1F(("m_eg_all_"+s).c_str(),";m_{#gamma#gamma};a.u.", 100, 0, 200);
	h["met_all"]	= TH1F(("met_all_"+s).c_str(),";#slash{E}_{T};counts",100,0,200);
	h["met_eg"]	= TH1F(("met_eg_"+s).c_str(),";#slash{E}_{T} for e#gamma;counts",100,0,200);
	h["met_eg_all"]	= TH1F(("met_eg_all_"+s).c_str(),";#slash{E}_{T} for e#gamma;counts",100,0,200);
	h["pt_eg"]	= TH1F(("pt_eg_"+s).c_str(),";fake photon p_{T} ",100,0,200);
	
	// invariant mass in states with 2 photon like objects with 1 object pT>25 GeV, the other pT>10 (called ptcut25)
	h["m_eg_ptcut25"] = TH1F(("m_eg_ptcut25_"+s).c_str(),";m [GeV];counts",100,0,200);	// all probes, eg
	h["m_gg_ptcut25"] = TH1F(("m_gg_ptcut25_"+s).c_str(),";m [GeV];counts",200,0,200);		// photon candidate probes, gg

	// Ht from jets in states with 2 photon like objects with 1 object pT>25 GeV, the other pT>10 (called ptcut25)
	h["Htjets_eg_ptcut25"] = TH1F(("Htjets_eg_ptcut25_"+s).c_str(),";H_T;counts",100,0,1000);	// all probes, eg
	h["Htjets_gg_ptcut25"] = TH1F(("Htjets_gg_ptcut25_"+s).c_str(),";H_T;counts",200,0,1000); 	// photon candidate probes, gg

	//h["Htjets_gg_ptcut25"] = TH1F(("Htjets_gg_ptcut25_"+s).c_str(),";H_T;counts",200,0,1000); 	// photon candidate probes, gg

	// use electron and photon objects; clean up by comparing momenta if objects are identified in both categories
	// then calculate invariant mass out of all combinations with one electron and one photon object.
	// no photon requirements: denominator, photon with no pixelseed: numerator
	h["m_ptcut25_0_den"] = TH1F(("m_ptcut25_0_den_"+s).c_str(),";m [GeV];counts",100,0,200);
	h["m_ptcut25_0_num"] = TH1F(("m_ptcut25_0_num_"+s).c_str(),";m [GeV];counts",100,0,200);
		
	// use electron and photon objects; clean up by comparing momenta if objects are identified in both categories
	// then find best z match from one electron object and one photon object
	// no photon requirements: denominator, photon with no pixelseed: numerator
	h["m_ptcut25_1_den"] = TH1F(("m_ptcut25_1_den_"+s).c_str(),";m [GeV];counts",100,0,200);
	h["m_ptcut25_1_num"] = TH1F(("m_ptcut25_1_num_"+s).c_str(),";m [GeV];counts",100,0,200);
	
	// simple fake rate calculation f = Ng / (Ng+Ne)
	h["Ng_Pt"] = TH1F(("Ng_Pt_"+s).c_str(),";p_{T} [GeV];fake rate",100,0,200);
	h["Ng+Ne_Pt"] = TH1F(("Ng+Ne_Pt_"+s).c_str(),";p_{T} [GeV];fake rate",100,0,200);
	
	// count photons
	h["number_photons_no_ptcut25"] = TH1F(("number_photons_no_ptcut25_"+s).c_str(),";number;counts",10,0,10);
	h["number_photons_ptcut25"] = TH1F(("number_photons_ptcut25_"+s).c_str(),";number;counts",10,0,10);
	
	// fakerate in dependency of pT
	h["f_Pt_total"] = TH1F(("f_Pt_total_"+s).c_str(),";P_{T} [GeV];counts",200,0,200);
	h["f_Pt_passed"] = TH1F(("f_Pt_passed_"+s).c_str(),";P_{T} [GeV];counts",200,0,200);
	
	// fakerate in dependency of Nvtx
	h["f_Nvtx_total"] = TH1F(("f_Nvtx_total_"+s).c_str(),";N_{vtx};counts",40,0,40);
	h["f_Nvtx_passed"] = TH1F(("f_Nvtx_passed_"+s).c_str(),";N_{vtx};counts",40,0,40);
	
	// fakerate in dependency of number of jets
	h["f_JetSize_total"] = TH1F(("f_JetSize_total_"+s).c_str(),";N_{jets};counts",40,0,40);
	h["f_JetSize_passed"] = TH1F(("f_JetSize_passed_"+s).c_str(),";N_{jets};counts",40,0,40);
	
	// fakerate in dependency of eta
	h["f_Eta_total"] = TH1F(("f_Eta_total_"+s).c_str(),";|#eta|;counts",200,0,5);
	h["f_Eta_passed"] = TH1F(("f_Eta_passed_"+s).c_str(),";|#eta|;counts",200,0,5);

	// fakerate in dependency of gen Ht
	h["f_met_total"] = TH1F(("f_met_total_"+s).c_str(),";#slash{E}_{T} [GeV];counts",100,0,500);
	h["f_met_passed"] = TH1F(("f_met_passed_"+s).c_str(),";#slash{E}_{T} [GeV];counts",100,0,500);

	// fakerate in dependency of HLT photon 90 CaloIdL PFHT500 trigger
	h["f_met_raw_total"] = TH1F(("f_met_raw_total_"+s).c_str(),";raw #slash{E}_{T} [GeV];counts",2,0,2);
	h["f_met_raw_passed"] = TH1F(("f_met_raw_passed_"+s).c_str(),";raw #slash{E}_{T} [GeV]",2,0,2);
	
	h["Nvtx"] = TH1F(("Nvtx_"+s).c_str(),";N_{vtx};a.u.",100,0,100);

	h["m_gg"].Sumw2();
	h["m_eg"].Sumw2();
	h["m_gg_all"].Sumw2();
	h["m_eg_all"].Sumw2();
	h["met_all"].Sumw2();
	h["met_eg"].Sumw2();
	h["met_eg_all"].Sumw2();
	h["pt_eg"].Sumw2();
	h["Nvtx"].Sumw2();
	h["m_eg_ptcut25"].Sumw2();
	h["m_gg_ptcut25"].Sumw2();
	
	h["Ng_Pt"].Sumw2();
	h["Ng+Ne_Pt"].Sumw2();
	
	h["f_Pt_total"].Sumw2();
	h["f_Pt_passed"].Sumw2();
	
	h["f_Nvtx_total"].Sumw2();
	h["f_Nvtx_passed"].Sumw2();
	
	h["f_JetSize_total"].Sumw2();
	h["f_JetSize_passed"].Sumw2();
	
	h["f_Eta_total"].Sumw2();
	h["f_Eta_passed"].Sumw2();

	h["f_met_total"].Sumw2();
	h["f_met_passed"].Sumw2();
	
	h["f_met_raw_total"].Sumw2();
	h["f_met_raw_passed"].Sumw2();
	
	h["delR"] = TH1F(("delR_"+s).c_str(),";#Delta R;a.u.",50,0,10);
	h["delR"].Sumw2();
	
	
	// *********************************************************************************************
	// 
	// 
	h["Zmass_fit_total_sig_noBin"] = TH1F(("Zmass_fit_total_sig_noBin_"+s).c_str(),";m [GeV];counts",400,0,200);
	h["Zmass_fit_total_bkg_noBin"] = TH1F(("Zmass_fit_total_bkg_noBin_"+s).c_str(),";m [GeV];counts",400,0,200);
	h["Zmass_fit_total_bkg2_noBin"] = TH1F(("Zmass_fit_total_bkg2_noBin_"+s).c_str(),";m [GeV];counts",400,0,200);
	h["Zmass_fit_total_bkg3_noBin"] = TH1F(("Zmass_fit_total_bkg3_noBin_"+s).c_str(),";m [GeV];counts",400,0,200);
	h["Zmass_fit_total_sig_noBin"].Sumw2();
	h["Zmass_fit_total_bkg_noBin"].Sumw2();
	h["Zmass_fit_total_bkg2_noBin"].Sumw2();
	h["Zmass_fit_total_bkg3_noBin"].Sumw2();
	
	
	h["Zmass_fit_passed_sig_noBin"] = TH1F(("Zmass_fit_passed_sig_noBin_"+s).c_str(),";m [GeV];counts",400,0,200);
	h["Zmass_fit_passed_bkg_noBin"] = TH1F(("Zmass_fit_passed_bkg_noBin_"+s).c_str(),";m [GeV];counts",400,0,200);
	h["Zmass_fit_passed_sig_noBin"].Sumw2();
	h["Zmass_fit_passed_bkg_noBin"].Sumw2();
	
	
	// TH2F (const char *name, const char *title, 
	//	Int_t nbinsx, Double_t xlow, Double_t xup, 
	//	Int_t nbinsy, Double_t ylow, Double_t yup)
	h2["delR_plane"] = TH2F(("delR_plane_"+s).c_str(),";#Delta #phi;#Delta #eta",
					400,-7.,7.,
					800,-10,10);
					
	h2["PhiEtaPlane_total"] = TH2F(("PhiEtaPlane_total_"+s).c_str(),";#Delta #phi;#Delta #eta",
					400,-7.,7.,
					800,-10,10);
	
	h2["PhiEtaPlane_passed"] = TH2F(("PhiEtaPlane_passed_"+s).c_str(),";#Delta #phi;#Delta #eta",
					400,-7.,7.,
					800,-10,10);
	
	h2["PhiEtaPlaneZoom_total"] = TH2F(("PhiEtaPlaneZoom_total_"+s).c_str(),";#Delta #phi;#Delta #eta",
					200,-1.,1.,
					200,-1,1);
	
	h2["PhiEtaPlaneZoom_passed"] = TH2F(("PhiEtaPlaneZoom_passed_"+s).c_str(),";#Delta #phi;#Delta #eta",
					200,-1.,1.,
					200,-1,1);
					
	
	h2["f_NtrkPt_total"] = TH2F(("f_NtrkPt_total_"+s).c_str(),";p_{T} [GeV];N_{trk}",
					200,0,200,
					250,0,250);

	h2["f_NtrkPt_passed"] = TH2F(("f_NtrkPt_passed_"+s).c_str(),";p_{T} [GeV];N_{trk}",
					200,0,200,
					250,0,250);
	
	
	fOutput->AddAll(gDirectory->GetList());
	
} // init selection


/*******************************************************************************************
 * 
 * initialize the tree reader
 *
 ************************************************/
void Histogrammer::Init(TTree *tree_)
{
	cout << "Init()" << endl;

	fReader.SetTree(tree_);
	
	cout << "end of Init()" << endl;
	
}



/*******************************************************************************************
 * 
 * reset the selections
 *
 ************************************************/
void Histogrammer::resetSelection() {
	
	selHt = 0;
	selPhotons.clear();
	selJets.clear();

	selBJets.clear();
	selElectrons.clear();
	selMuons.clear();
	
	ht = 0;


		selPhotons2.clear();
		selJets2.clear();
		selBJets2.clear();
		selElectrons2.clear();
		selMuons2.clear();	
	
	// INIT
		// DONE (move this section to an initial function outside process()!) DONE 
		// /************************************************************************
		sump_gg=0.;
		sump_eg=0.;
		pt1=0.;
		pt2=0.;
		eta1=0.;
		eta2=0.;
		phi1=0.;
		phi2=0.;
		Msquare=0.;
		Htjets=0.;
		Etemp1=0.;
		Etemp2=0.;			

		oneHasNoPixelSeed = false;	// set this TRUE if one has NO pixelseed and
						// therefore seems to be a "real" photon

		//TagElectron = false;
		
		vTemp1.SetXYZ(0.,0.,0.);			// temp vector
		vTemp2.SetXYZ(0.,0.,0.);			// temp vector
		vTemp3.SetXYZ(0.,0.,0.);			// temp vector
		vTemp4.SetXYZ(0.,0.,0.);			// temp vector		

		lvTemp1.SetPxPyPzE(0.,0.,0.,0.);		// temp lorentzvector
		lvTemp2.SetPxPyPzE(0.,0.,0.,0.);		// temp lorentzvector
		lvTemp3.SetPxPyPzE(0.,0.,0.,0.);		// temp lorentzvector
		lvTemp4.SetPxPyPzE(0.,0.,0.,0.);		// temp lorentzvector

		vec_lvTemp1.clear();	// clear lvectors of lorentzvectors
		vec_lvTemp2.clear();
		vec_lvTemp3.clear();
		vec_lvTemp4.clear();
		
		vecTemp1.clear();
		vecTemp2.clear();
		vecTemp3.clear();
		vecTemp4.clear();
		
		BestZMatch_ElectronPhoton_Denom.clear();
		BestZMatch_ElectronPhoton_Num.clear();
		
		M_ElectronPhoton.clear();
		
		// ************************************************************************
		// */		


	
	//mrr2.first = 0;
	//mrr2.second = 0;

	//selJets.clear();
	//selBJets.clear();
	//selElectrons.clear();
	//selMuons.clear();

	
} // reset selection

/*******************************************************************************************
 * 
 * begin
 *
 ************************************************/
void Histogrammer::Begin(TTree * /*tree*/)
{
	// The Begin() function is called at the start of the query.
	// When running with PROOF Begin() is only called on the client.
	// The tree argument is deprecated (on PROOF 0 is passed).

	cout << "Begin()" << endl;

	//TString option = GetOption();
	//cout << "Option: " << option << endl;
	//Info("Begin", "starting a simple exercise with process option: %s", option.Data());
	//printf("\n");
	

}

/*******************************************************************************************
 * 
 * slave begin
 *
 ************************************************/
void Histogrammer::SlaveBegin(TTree *tree_)
{
	cout << "SlaveBegin()" << endl;
	
	
	
	// initalize the selection histograms	
	//initSelection("nocut");
	
	initSelection("test");
	
} // slave begin



/*******************************************************************************************
 * 
 * process.
 *
 ************************************************/
Bool_t Histogrammer::Process(Long64_t entry)
{
	if(Nnum > 100000) return kTRUE;
			
	
	//if(Nnum<10) cout << Nnum << " Process()-"; // << endl;
	
	resetSelection();
	
	fReader.SetEntry(entry);	// fReader on current entry
	
	
	//cout << "-"<<entry;		

	// set weight
	//selW = sign(*w_mc);
	
	//float Nvertices; // number of vertices
	//Nvertices = nGoodVertices;	
/*		
	// scan photon objects per event		
	for(auto& p: photons){
					
		selPhotons.push_back(&p);
		vTemp1 = vTemp1 + p.p;	// summed photon momentum
		Etemp1 += p.p.Mag();	// summed photon energy
								
	}

	// scan jets per event
	for(auto& j: jets){
		Htjets += j.p.Pt();

	}
	
	//scan electron objects per event
	for(auto& e: electrons){
		
		selElectrons.push_back(&e);

	} //

	// scan muon objects per event
	for(auto& m: muons){
		
		selMuons.push_back(&m);
	
	}
	
	// set lorentzvector from the sum of momentum three vectors and energy, respectively
	lvTemp1.SetVect(vTemp1);// summed momentum
	lvTemp1.SetE(Etemp1);	// summed energy
		
	if(lvTemp1.M()>1) h["m_gg_all"].Fill(lvTemp1.M());	// photon candidate probes, "gg" = eg
	if(lvTemp2.M()>1) h["m_eg_all"].Fill(lvTemp2.M());	// all probes,		  "eg" = eg + ee
	
	// only final states with 2 photon like objects
	if(selPhotons.size()==2){
					
		h["m_gg"].Fill(lvTemp1.M());
		h["met_eg_all"].Fill(met->p_raw.Mag());
		
		// either one or the other has NO pixel seed
		if( 	(selPhotons[0]->hasPixelSeed && !selPhotons[1]->hasPixelSeed) || 	
			(selPhotons[1]->hasPixelSeed && !selPhotons[0]->hasPixelSeed)
			){

			oneHasNoPixelSeed = true;	// control variable for pixelseed veto
							
			h["m_eg"].Fill(lvTemp1.M());				
			h["met_eg"].Fill(met->p_raw.Mag());
			
			// pT of the electron proxy objects (
			if(selPhotons[0]->hasPixelSeed) h["pt_eg"].Fill(selPhotons[0]->p.Pt());
			else if(selPhotons[1]->hasPixelSeed) h["pt_eg"].Fill(selPhotons[1]->p.Pt());

		}
		
		// ptcut25 for at least one object, ptcut10 for the other
		if(	(selPhotons[0]->p.Pt()>25. && selPhotons[1]->p.Pt()>10.) ||
			(selPhotons[0]->p.Pt()>10. && selPhotons[1]->p.Pt()>25.)
			){
			
			// invariant mass
			h["m_eg_ptcut25"].Fill(lvTemp1.M());				// all probes, "eg" = eg + ee
			
			
			// Ht from jets
			h["Htjets_eg_ptcut25"].Fill(Htjets);				// all probes, "eg" = eg + ee

			if(oneHasNoPixelSeed) {
				h["m_gg_ptcut25"].Fill(lvTemp1.M());	// photon candidate probes, "gg" = eg
				h["Htjets_gg_ptcut25"].Fill(Htjets);	// photon candidate probes, "gg" = eg
			}
			
		}
					
	} // two photon objects
			
	// met
	h["met_all"].Fill(met->p_raw.Mag());
	
	//wpu = *w_pu;

	// number of "good" vertices per event
	h["Nvtx"].Fill(*nGoodVertices);
	
	//h["Nvtx"].Fill(w_pu);
	//
*/
	// *******************************************************************************
	// *******************************************************************************
	// tag and probe in the double photon data set ***********************************
	//
	
	// RESET
	resetSelection();
	
	bool	TagElectron = false,
		ProbePhoton = false;
		
	int	Nphotons_nocut = 0,
		Nphotons_cut = 0;
		
	int	Ng_temp=0,
		Ne_temp=0;
	
	float 	epsilonP = 0.1; //maximum difference between two absolute momenta
				 // to distinguisch originally same objects
				 // (in GeV)
				 
	float	delR_p_gene = 0,	// dR photon object <-> generated electron
		delR_e_gene = 0;
	
	bool	oneDelRmatch = false;
	 
	float 	epPt_temp=0;
	
	S_mZdiff_m mz_m1, mz_m2;
	
	float 	Ptcut_El = 25.; // 25GeV as cut
	
/*		deltaR	 = sqrt(	(genp.p.Phi()-p.p.Phi())*(genp.p.Phi()-p.p.Phi()) + 
		    		(genp.p.Eta()-p.p.Eta())*(genp.p.Eta()-p.p.Eta()) 
		    		);
				
*/
	
	// *******************************************************************************
	// simple fake rate: f=Ng/(Ne+Ng)
	// count Ng and Ne
	bool matchToPF = false;
	bool gotoNextPhoton = false;
	
	for(auto& p: photons){ // loop photon objects
		
		vTemp1 = vTemp1 + p.p;	// summed photon momentum
		Etemp1 += p.p.Mag();	// summed photon energy
	
		gotoNextPhoton = false;
		
		for(auto& genp: genParticles){ // per photon object loop gen particles
			if( (genp.pdgId == 11) || (genp.pdgId == -11)){ // if electron or positron
				Ntest2 += 1;
				
				// TH2F.Fill (Double_t x, Double_t y)
				h2["delR_plane"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
				h["delR"].Fill(genp.p.DeltaR(p.p));
				
				if(p.p.Pt() > 30. && fabs(p.p.Eta()) < 1.4442){ // fake rate in the plane
					if(!p.hasPixelSeed){
						h2["PhiEtaPlane_total"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
						h2["PhiEtaPlane_passed"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
						
						h2["PhiEtaPlaneZoom_total"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
						h2["PhiEtaPlaneZoom_passed"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
						
					}
					if(p.hasPixelSeed){
						h2["PhiEtaPlane_total"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
						h2["PhiEtaPlaneZoom_total"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
					}
				}
				
				if(genp.p.DeltaR(p.p) <0.2 && p.p.Pt() >30.){
					if(!p.hasPixelSeed){
						h["f_Eta_total"].Fill(fabs(p.p.Eta()));
						h["f_Eta_passed"].Fill(fabs(p.p.Eta()));
					}
					if(p.hasPixelSeed){
						h["f_Eta_total"].Fill(fabs(p.p.Eta()));
					}
				}
				
				if(!gotoNextPhoton){ // if, for any chance, the delR requirement is fullfilled
						// more than one time, this makes sure not to count
						// photon objects twice
					if(genp.p.DeltaR(p.p) <0.2 && p.p.Pt() > 30. && fabs(p.p.Eta()) < 1.4442 ){
						if(!p.hasPixelSeed  ){	
							Ng_1 +=1;
										
							h["f_Pt_total"].Fill(p.p.Pt());
							h["f_Pt_passed"].Fill(p.p.Pt());
							
							h["f_Nvtx_total"].Fill(*nGoodVertices);
							h["f_Nvtx_passed"].Fill(*nGoodVertices);
							
							h["f_JetSize_total"].Fill(jets.GetSize());
							h["f_JetSize_passed"].Fill(jets.GetSize());
							
							h["f_met_total"].Fill(met->p.Pt());
							h["f_met_passed"].Fill(met->p.Pt());
							
							//h2["f_NtrkPt_total"].Fill(p.p.Pt(),*nChargedPfCandidates);
							//h2["f_NtrkPt_passed"].Fill(p.p.Pt(),*nChargedPfCandidates);
							
							//h2["PhiEtaPlane_total"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
							//h2["PhiEtaPlane_passed"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
							
						}
						if(p.hasPixelSeed){
							Ne_1 +=1;
								
							h["f_Pt_total"].Fill(p.p.Pt());
							h["f_Nvtx_total"].Fill(*nGoodVertices);
							h["f_JetSize_total"].Fill(jets.GetSize());
							
							h["f_genHt_total"].Fill(*genHt);
							h["f_met_total"].Fill(met->p.Mag());
							
							//h2["f_NtrkPt_total"].Fill(p.p.Pt(),*nChargedPfCandidates);		
							//h2["PhiEtaPlane_total"].Fill(genp.p.Phi()-p.p.Phi(),genp.p.Eta()-p.p.Eta());
						}
						//goto nextPhoton;
						gotoNextPhoton = true;
					}
				}
			}
		} // genparticle loop
		// nextPhoton:
		
		
		
		
		
	} // photon loop
	
	// *********************************************************************************************
	// try to seperate signal and background by simple assumptions:
	//	use all events with 2 photon objects
	//	- both photons match to gen particles: "signal"
	//	- else: "background"
	
	// set lorentzvector from the sum of momentum three vectors and energy, respectively
	lvTemp1.SetVect(vTemp1);// summed momentum
	lvTemp1.SetE(Etemp1);	// summed energy
	
	int nPhoton = 0;
	bool 	photon1 = false, photon2 = false, // match variables
			hasPixelSeed1 = false, hasPixelSeed2 = false; // pixelseed
	if(photons.GetSize() == 2){
		for(auto& p: photons){
			nPhoton += 1;
			for(auto& genp: genParticles){
				if(genp.p.DeltaR(p.p) <0.2 && p.p.Pt() > 30.){
					
					if(nPhoton == 1){ 
						if(p.hasPixelSeed) hasPixelSeed1 = true;
						photon1 = true;
						}
					if(nPhoton == 2){
						if(p.hasPixelSeed) hasPixelSeed2 = true;
						photon2 = true;
						}
						
					goto nextPhoton; // do not count photons twice
					
				}
			}
		nextPhoton:
		Ntest1 += 1;
		}
	}
	
	if(photon1 && photon2){
		h["Zmass_fit_total_sig_noBin"].Fill(lvTemp1.M());
		if( (hasPixelSeed1 && !hasPixelSeed2) || (!hasPixelSeed1 && hasPixelSeed2) ){
			h["Zmass_fit_passed_sig_noBin"].Fill(lvTemp1.M());
			}
		Ntest2 += 1; 
	}
	
	if( (photon1  && !photon2) || (!photon1 && photon2) ){ // || (!photon1 && !photon2) ){
		
		if(photon1 && !photon2) h["Zmass_fit_total_bkg_noBin"].Fill(lvTemp1.M());
		if(!photon1 && photon2) h["Zmass_fit_total_bkg2_noBin"].Fill(lvTemp1.M());
		if( (hasPixelSeed1 && !hasPixelSeed2) || (!hasPixelSeed1 && hasPixelSeed2) ){
			h["Zmass_fit_passed_bkg_noBin"].Fill(lvTemp1.M());
			}
		Ntest3 += 1;
	}
	
	if(!photon1 && !photon2){
		h["Zmass_fit_total_bkg3_noBin"].Fill(lvTemp1.M());
		Ntest4 += 1;
	}
	
/*	
	if(Nnum==4){
		for(auto& p: photons) cout << "p.Pt() = " << p.p.Pt() << " and p.Eta() = " << p.p.Eta() << endl;
		for(auto& g: genParticles){
			if(fabs(g.pdgId) == 11) cout << "gen.Pt() = " << g.p.Pt() << endl;
		}
	} 
*/		
	//cout << Nnum << "\t" << Ng_1 << "\t" << Ne_1 << endl;
	




	// comment this out to get whole tree analyzed:
	Nnum +=1 ;

	

	
	
	// */
	return kTRUE;
	
	

	// comment catcher */
	
} // process


/*******************************************************************************************
 * 
 * terminate slave function - 
 *
 ************************************************/
void Histogrammer::SlaveTerminate()
{
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.
	
	cout << "\nSlaveTerminate()" << endl;
	
	/*
	h["m_gg"].Scale(1/h["m_gg"].GetEntries());
	h["m_ge"].Scale(1/h["m_ge"].GetEntries());
	h["m_gg_all"].Scale(1/h["m_gg_all"].GetEntries());
	h["m_eg_all"].Scale(1/h["m_eg_all"].GetEntries());
	// */	

}


/*******************************************************************************************
 * 
 * terminate - last function to call
 *
 ************************************************/
void Histogrammer::Terminate()
{
	cout << "Terminate()" << endl;
	
	cout << "total: " << Nnum << endl;

	cout << "Ntest1: " << Ntest1 << endl;
	cout << "Ntest2: " << Ntest2 << endl;
	cout << "Ntest3: " << Ntest3 << endl;
	cout << "Ntest4: " << Ntest4 << endl;
	
	cout << "Ng: " << Ng_1 << endl;
	cout << "Ne: " << Ne_1 << endl;
	cout << "f = " << Ng_1/(Ng_1+Ne_1) << endl;
	
		
	//fOutput->Add(&g["Pt"]);*/
	
	//fOutput->Add(eff["Pt"]);
	
	//eff["Pt"].SetDirectory(gDirectory);
	//fOutput->AddAll(gDirectory->GetList());
	
	auto outputName = getOutputFilename( fReader.GetTree()->GetCurrentFile()->GetName() );
	
	cout << "results written in: " << outputName << endl;
	
	TFile file( outputName.c_str(), "RECREATE");
	
	//TFile hfile("../rootFiles/my_selector_results_mcDY_1446802738.root","RECREATE");
	
	
	//TFile hfile("../rootFiles/my_selector_results_DY_delR_cut.root","RECREATE");
	//TFile hfile("my_selector_results_DY_FullSim.root","RECREATE");
	
	
// */
	fOutput->Write(); // write all objects; used by PROOF	
} // terminate




























