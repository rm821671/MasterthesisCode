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
#include <iostream>
#include <fstream>

#include <algorithm>
#include <string>

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
#include "../../cms_sw/CMSSW_7_4_14/src/TreeWriter/TreeWriter/plugins/TreeParticles.hpp"

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
	
	struct PFCandidate : public Particle
   {
      Int_t pdgId=0;      
      Int_t charge;
      
      // 0: no match
      // 1: electron/positron id-charge match
      // 2: pi+/pi- id-charge match
      //Int_t IdChargeMatch;
      
      // primary vertex information
      //
      // enum PVAssociationQuality {
      // 	NotReconstructedPrimary=0,
      // 	OtherDeltaZ=1,
      // 	CompatibilityBTag=4,
      // 	CompatibilityDz=5,
      // 	UsedInFitLoose=6,
      // 	UsedInFitTight=7	
      // };
      Int_t pvAssociationQuality;
      
      
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
const int sizeNcounts = 20;
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

  // Converts "/path/to/ntuple/XYZ.root" 
  // to "../selector_rootFiles/my_selector_results_XYZ_(timestamp).root"



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

void replaceInString(string &s) {
	std::replace( s.begin(), s.end(), '.', 'd'); // replace all '.' to 'd'
}


// create histograms
//
// --- --- pass




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
		TTreeReaderArray<tree::PFCandidate> pfCandidates;

		TTreeReaderValue<tree::MET> met;
		
		TTreeReaderValue<Float_t> pu_weight;
		TTreeReaderValue<Char_t> mc_weight;
		
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
		
		TTreeReaderValue<Bool_t> TriggerR9Id90;
		TTreeReaderValue<Bool_t> TriggerR9Id85;

		vector<tree::Photon*> selPhotons;
		vector<tree::Jet*> selJets;
		vector<tree::Jet*> selBJets;
		vector<tree::Electron*> selElectrons;
		vector<tree::Muon*> selMuons;
		vector<tree::PFCandidate*> selPFCandidates;
		
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
		
		
		
		int nTotalEvents = 0;
		
		
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
		
		//const int sizeNcounts = 20;
		
		int Ncounts[sizeNcounts] = {};
		
		
		
		int N_pfc = 0;
		

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
		//map<string,vector<float>> Ne, Ng, f, Xval;
		
		//struct
		//map<string,float> Ne, Ng;
		
		// ************************************************************************
		// */
		
		
		bool isData;
		
		
		//pair<float,float> mrr2;
		
		map<string,TEfficiency> eff;	// efficiencies
		map<string,TH1F> h;				// 1d histograms
		map<string,TGraph> g;			// graphs , NOT USED
		
		map<string,TH2F> h2;			// 2d histograms
		
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
	pfCandidates( fReader, "PFCandidates" ),
	met( fReader, "met" ),
	nGoodVertices( fReader, "nGoodVertices" ),
	nChargedPfCandidates( fReader, "nChargedPfCandidates" ),
	mc_weight( fReader, "mc_weight" ),
	pu_weight( fReader, "pu_weight" ),
	
	genHt( fReader, "genHt"),
	
	signalTriggerFired( fReader, "HLT_Photon90_CaloIdL_PFHT500_v" ),
	TriggerR9Id90( fReader, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_v" ),
	TriggerR9Id85( fReader, "HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15_v" )
	
	//isRealData( fReader, "isRealData" ),
	//signalTriggerFired( fReader, "HLT_Photon90_CaloIdL_PFHT500_v" ),
	//crossTriggerPhotonFired( fReader, "HLT_Photon90_v" ),
	//crossTriggerHtFired( fReader, "HLT_PFHT600_v" )
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
	
	
	//eff["Pt"] = TEfficiency(("eff_Pt_"+s).c_str(),"my efficiency;p;fakerate",200,0,200);
	//eff["Pt"] = ->SetDirectory(gDirectory);
	
	// TH1F(const char* name, const char* title, nbins, binmin, binmax)
	
	
	//h[""] 	= TH1F((""+s).c_str(),";;",nbins,min,max);
/*	
	string hname;
	string numstr;
	for(float ptbin = 0.; ptbin <= 10.; ptbin += 0.05){
		numstr = NumberToString(ptbin);
		replaceInString(numstr);
		hname = "PFC_n_Ptcut" + numstr;
		h[hname] = TH1F((hname + "_" + s).c_str(),";Number;counts",300,0,300);
		h[hname].Sumw2();
	}
*/	
	
	// bins of Pt
	
	
/*
	string hname;
	string numstr;
	for(int ptbin = 20; ptbin < 200; ptbin += 5){ // goes to 195
		numstr = NumberToString(ptbin) + "-" + NumberToString(ptbin + 5);
		hname = "fZmass_tot_sig_Ptbin_" + numstr;
		h[hname] = TH1F((hname + "_" + s).c_str(),";Number;counts",200,0,200);
		h[hname].Sumw2();
		hname = "fZmass_pas_sig_Ptbin_" + numstr;
		h[hname] = TH1F((hname + "_" + s).c_str(),";Number;counts",200,0,200);
		h[hname].Sumw2();
		
	}
	
	// Vertex dependency
	for(int Nvtx = 1; Nvtx < 40; Nvtx += 1){ // 39
		numstr = NumberToString(Nvtx);
		hname = "fZmass_tot_sig_Nvtx_" + numstr;
		h[hname] = TH1F((hname + "_" + s).c_str(),";Number;counts",200,0,200);
		h[hname].Sumw2();
		hname = "fZmass_pas_sig_Nvtx_" + numstr;
		h[hname] = TH1F((hname + "_" + s).c_str(),";Number;counts",200,0,200);
		h[hname].Sumw2();
		
	}
	
	// number of jets dependency
	for(int Njet = 1; Njet < 30; Njet += 1){ // 29
		numstr = NumberToString(Njet);
		hname = "fZmass_tot_sig_Njet_" + numstr;
		h[hname] = TH1F((hname + "_" + s).c_str(),";Number;counts",200,0,200);
		h[hname].Sumw2();
		hname = "fZmass_pas_sig_Njet_" + numstr;
		h[hname] = TH1F((hname + "_" + s).c_str(),";Number;counts",200,0,200);
		h[hname].Sumw2();
		
	}
*/
	
	
	
	// h["PFC_n_Ptcut5"] 	= TH1F(("PFC_n_Ptcut5_"+s).c_str(),";Number;counts",300,0,300);
	// h["PFC_n_Ptcut5"].Sumw2();
	
	string hname;
	
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
	
	
	hname = "rawdata_photon_pt";
	h[hname] = TH1F((hname+"_"+s).c_str(),";pt [GeV];counts",300,0,300);
	
	hname = "rawdata_electron_pt";
	h[hname] = TH1F((hname+"_"+s).c_str(),";pt [GeV];counts",300,0,300);
	
	hname = "rawdata_jet_pt";
	h[hname] = TH1F((hname+"_"+s).c_str(),";pt [GeV];counts",300,0,300);
	
	hname = "diphotonraw_photon_pt";
	h[hname] = TH1F((hname+"_"+s).c_str(),";pt [GeV];counts",300,0,300);
	
	hname = "diphotonraw_electron_pt";
	h[hname] = TH1F((hname+"_"+s).c_str(),";pt [GeV];counts",300,0,300);
	
	hname = "diphotonraw_jet_pt";
	h[hname] = TH1F((hname+"_"+s).c_str(),";pt [GeV];counts",300,0,300);
	
	hname = "diphotonZpeak_jet_pt";
	h[hname] = TH1F((hname+"_"+s).c_str(),";m [GeV];counts",200,0,200);
	
	
	
	hname = "diphoton_excess";
	h[hname] = TH1F((hname+"_"+s).c_str(),";m [GeV];counts",2000,0,2000);
	
	
	
	// TH2F (const char *name, const char *title, 
	//	Int_t nbinsx, Double_t xlow, Double_t xup, 
	//	Int_t nbinsy, Double_t ylow, Double_t yup)
	string h2name;
	
	h2name = "Ntrk_Pt_plane";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";Number;Pt",
						300,0,300,
						100,0,10);
	
	h2name = "fZmass_h2_tot_pt";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];pt [GeV]",
						200,0,200,
						250,0,250);
	
	h2name = "fZmass_h2_pas_pt";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];pt [GeV]",
						200,0,200,
						250,0,250);
	
	h2name = "fZmass_h2_tot_nvtx";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];N_{vtx}",
						200,0,200,
						60,0,60);

	h2name = "fZmass_h2_pas_nvtx";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];N_{vtx}",
						200,0,200,
						60,0,60);

	h2name = "fZmass_h2_tot_njet";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];N_{jet}",
						200,0,200,
						15,0,15);

	h2name = "fZmass_h2_pas_njet";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];N_{jet}",
						200,0,200,
						15,0,15);
	
	h2name = "h2_2e2g_delR";
	
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";#Delta #phi;#Delta #eta",
						320,-3.2,3.2,
						320,-8,8);
						
	hname = "h_2e2g_delR";
	h[hname] = TH1F((hname+"_"+s).c_str(),";del R;counts",100,0,5);
	
	// for electron tag
	h2name = "fZmass_h2_tot_pt_etag";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];pt [GeV]",
						200,0,200,
						250,0,250);
	
	h2name = "fZmass_h2_pas_pt_etag";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];pt [GeV]",
						200,0,200,
						250,0,250);
	
	h2name = "fZmass_h2_tot_nvtx_etag";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];N_{vtx}",
						200,0,200,
						60,0,60);

	h2name = "fZmass_h2_pas_nvtx_etag";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];N_{vtx}",
						200,0,200,
						60,0,60);

	h2name = "fZmass_h2_tot_njet_etag";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];N_{jet}",
						200,0,200,
						15,0,15);

	h2name = "fZmass_h2_pas_njet_etag";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";m [GeV];N_{jet}",
						200,0,200,
						15,0,15);
	
	
	
	
	
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
	string thisFileName = fReader.GetTree()->GetCurrentFile()->GetName();
	cout << "opened file: " << thisFileName << endl;
	
	// if mc, set true, if data set false
	isData = 	thisFileName.find("DYJetsToLL") != string::npos
			||	thisFileName.find("TTGJets") != string::npos
			||	thisFileName.find("WGToLNuG") != string::npos
			||	thisFileName.find("ZGTo2LG") != string::npos;
	
	
	TFile thisFile(thisFileName.c_str());
	if(!thisFile.IsZombie()){
		TH1F* hcut = (TH1F*)thisFile.Get("TreeWriter/hCutFlow");
		nTotalEvents = hcut->GetBinContent(1);
	}
	
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
	
		selPFCandidates.clear();
	
		ht = 0;
		
		N_pfc = 0;
		


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
		
		vTemp1.SetXYZ(0.,0.,0.);				// temp vector
		vTemp2.SetXYZ(0.,0.,0.);				// temp vector
		vTemp3.SetXYZ(0.,0.,0.);				// temp vector
		vTemp4.SetXYZ(0.,0.,0.);				// temp vector
		
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
	// comment this out to get whole tree analyzed:
	//if(Nnum > 1000000) return kTRUE;
	
	string numstr;
	string hname;
	string h2name;
	sdf
	resetSelection();	
	fReader.SetEntry(entry);	// fReader points on current entry
	
	// set weight
	selW = *mc_weight * *pu_weight;
	
	// *******************************************************************************
	// *******************************************************************************
	// try to construct number of tracks
		
	// loop the pf candidates for pt cuts between 0 and 50
/*
	string h2name = "Ntrk_Pt_plane";
	
	string hname;
	string numstr;
	for(float cut = 0.; cut <= 10.; cut+=0.05){	// make histograms with cuts
		resetSelection();
		numstr = NumberToString(cut);
		replaceInString(numstr);
		hname = "PFC_n_Ptcut" + numstr;
		
		for(auto& pfc : pfCandidates){
			if((pfc.pvAssociationQuality == 6 || pfc.pvAssociationQuality == 7) &&
				pfc.p.Pt() > float(cut)	)
				{
				//
				// only UsedInFitLoose=6,
				// 		UsedInFitTight=7
				selPFCandidates.push_back(&pfc);
				Ntest1 +=1;
			}
		}
		//h[hname].Fill(selPFCandidates.size());
		
		h2[h2name].Fill(selPFCandidates.size(),cut);
		Ntest2 +=1;
	}
*/

// comment catcher */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// raw data
	
	
	
	int numberPhoton40 = 0;
	
	hname = "rawdata_photon_pt";
	for(auto& p: photons){
		h[hname].Fill(p.p.Pt(), selW);
		
		if(p.p.Pt() > 40) numberPhoton40++;
	}
	
	hname = "rawdata_electron_pt";
	for(auto& e: electrons){
		h[hname].Fill(e.p.Pt(), selW);
	}
	
	hname = "rawdata_jet_pt";
	for(auto& j: jets){
		h[hname].Fill(j.p.Pt(), selW);
	}
	
	
	
	if(photons.GetSize() == 2 && numberPhoton40 >= 2 && *TriggerR9Id85){
		hname = "diphotonraw_photon_pt";
		for(auto& p: photons){
			
			h[hname].Fill(p.p.Pt(), selW);
			
			vTemp1 = vTemp1 + p.p;	// summed photon momentum
			Etemp1 += p.p.Mag();	// summed photon energy
		}
		
		hname = "diphotonraw_electron_pt";
		for(auto& e: electrons){
			h[hname].Fill(e.p.Pt(), selW);
		}
	
		hname = "diphotonraw_jet_pt";
		for(auto& j: jets){
			h[hname].Fill(j.p.Pt(), selW);
		}
		
		lvTemp1.SetVect(vTemp1);// summed momentum
		lvTemp1.SetE(Etemp1);	// summed energy
		hname = "diphotonZpeak_jet_pt"; // this is the photon pt!
		h[hname].Fill(lvTemp1.M(), selW);
	}
	
	





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// use 2d histograms

/*
	// booleans
	nPhoton = 0;
	photon1 = false;
	photon2 = false; // match variables
	hasPixelSeed1 = false;
	hasPixelSeed2 = false; // pixelseed
	
	tag1 = false; tag2 = false;
	probe1 = false; probe2 = false;

	string h2name;
*/
	int 	nPhoton = 0,
			nElectron = 0;
	bool 	photon1 = false, // <-- TAG
			photon2 = false, // <-- PROBE
							 // match variables
			hasPixelSeed1 = false, 
			hasPixelSeed2 = false; // pixelseed
	
	bool	tag1 = false, probe1 = false,
			tag2 = false, probe2 = false;
	
	
	resetSelection();
	for(auto& p: photons){ // loop photon objects
		vTemp1 = vTemp1 + p.p;	// summed photon momentum
		Etemp1 += p.p.Mag();	// summed photon energy
	}
	
	// set lorentzvector from the sum of momentum three vectors and energy, respectively
	lvTemp1.SetVect(vTemp1);// summed momentum
	lvTemp1.SetE(Etemp1);	// summed energy
	
	// use events with two photons
	if(photons.GetSize() == 2 && *TriggerR9Id85){
	Ncounts[6] ++;
	
	
	hname = "diphoton_excess";
	h[hname].Fill(lvTemp1.M(),selW);
	
	
	
		if(electrons.GetSize() == 2){
			Ncounts[5]++;
			
			// tag loop
			for(auto& p: photons){
				nPhoton += 1;
					cd
					for(auto& e: electrons){
						h2name = "h2_2e2g_delR";
						hname = "h_2e2g_delR";
						
						h[hname].Fill(p.p.DeltaR(e.p), selW);
						h2[h2name].Fill(p.p.DeltaPhi(e.p), p.p.Eta()-e.p.Eta(), selW);
						
						// tag requirements for an electron as tag
						if(e.p.Pt() > 40 && fabs(e.p.Eta() < 1.4442)){
							
							
						}
							
						
						if(p.p.DeltaR(e.p) > 0.1){
							// not the same object
							
						}
						
						
					}
				
					// these are the tag requirements with a photon as tag (pixelseed)
					if(	p.p.Pt() > 25.	&&	p.hasPixelSeed 	&&	fabs(p.p.Eta()) < 1.4442
						){
					
						if(nPhoton == 1){
							Ncounts[2] ++;
							
							tag1 = true;
							}
						if(nPhoton == 2){
							Ncounts[3] ++;
							
							tag2 = true;
							}				
						goto nextPhoton; // do not count both photons as tag
					}
			
				nextPhoton:
				Ncounts[0] += 1;
			}
			
			
			//h2name = "fZmass_h2_tot_pt";
			//h2name = "fZmass_h2_pas_pt";
			
			// probe
			if(tag1 || tag2){
			Ncounts[7] ++;
			
			
			
				if(tag1 && photons[1].p.Pt() > 10. && fabs(photons[1].p.Eta())<1.4442){				
					Ncounts[1] += 1;
					
					h2name = "fZmass_h2_tot_pt";
					h2[h2name].Fill(lvTemp1.M(), photons[1].p.Pt(), selW);
					
					h2name = "fZmass_h2_tot_nvtx";
					h2[h2name].Fill(lvTemp1.M(), *nGoodVertices, selW);
					
					h2name = "fZmass_h2_tot_njet";
					h2[h2name].Fill(lvTemp1.M(), jets.GetSize(), selW);
					
					
					
					if(!photons[1].hasPixelSeed){
						Ncounts[8]++;
						h2name = "fZmass_h2_pas_pt";
						h2[h2name].Fill(lvTemp1.M(), photons[1].p.Pt(), selW);
					
						h2name = "fZmass_h2_pas_nvtx";
						h2[h2name].Fill(lvTemp1.M(), *nGoodVertices, selW);
					
						h2name = "fZmass_h2_pas_njet";
						h2[h2name].Fill(lvTemp1.M(), jets.GetSize(), selW);
					}
				
				}
				
				if(tag2 && photons[0].p.Pt() > 10. && fabs(photons[0].p.Eta())<1.4442){
					Ncounts[1]+=1;
					
					h2name = "fZmass_h2_tot_pt";
					h2[h2name].Fill(lvTemp1.M(), photons[1].p.Pt(), selW);
					
					h2name = "fZmass_h2_tot_nvtx";
					h2[h2name].Fill(lvTemp1.M(), *nGoodVertices, selW);
					
					h2name = "fZmass_h2_tot_njet";
					h2[h2name].Fill(lvTemp1.M(), jets.GetSize(), selW);
					
					if(!photons[0].hasPixelSeed){
						h2name = "fZmass_h2_pas_pt";
						h2[h2name].Fill(lvTemp1.M(), photons[0].p.Pt(), selW);
					
						h2name = "fZmass_h2_pas_nvtx";
						h2[h2name].Fill(lvTemp1.M(), *nGoodVertices, selW);
					
						h2name = "fZmass_h2_pas_njet";
						h2[h2name].Fill(lvTemp1.M(), jets.GetSize(), selW);
					}
				}
			}
		
		
		
		
		
		} // if(electrons.GetSize == 2)
	} // if(photons.GetSize() == 2)


// comment catcher */

	Nnum +=1 ;

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
	
	
}


/*******************************************************************************************
 * 
 * terminate - last function to call
 *
 ************************************************/
void Histogrammer::Terminate()
{
	auto openedDataFile = fReader.GetTree()->GetCurrentFile()->GetName();

	cout << "Terminate()" << endl;
	
	cout << "opened data file: " << openedDataFile << endl;
	cout << isData << endl;
	cout << "total counted: " << Nnum << endl;
	cout << "total from cutflow: " << nTotalEvents << endl;
	
	cout << "Ntest1: " << Ntest1 << endl;
	cout << "Ntest2: " << Ntest2 << endl;
	cout << "Ntest3: " << Ntest3 << endl;
	cout << "Ntest4: " << Ntest4 << endl;
	
	cout << "Number of real photons: Ng: " << Ng_1 << endl;
	cout << "Number of fake photons: Ne: " << Ne_1 << endl;
	
	float f;
	if(Ng_1<Ne_1) f = Ng_1/(Ng_1+Ne_1);
	else f = 0.;
	
	cout << "Fakerate f = Ng/(Ng+Ne) = " << f << endl;
	
	for(int ll=0;ll<sizeNcounts;ll++){
		cout << "Ncounts[" << ll << "]: " << Ncounts[ll] << "\t";
	}
	cout << endl;
	
	
	
	string comment;
	
	
	comment = "";
	
	
	
			
	//fOutput->Add(&g["Pt"]);*/
	
	//fOutput->Add(eff["Pt"]);
	
	//eff["Pt"].SetDirectory(gDirectory);
	//fOutput->AddAll(gDirectory->GetList());
	
	auto outputName = getOutputFilename( fReader.GetTree()->GetCurrentFile()->GetName() );
	
	
	
	
	// some time stuff
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	
	
	cout << "results written in: " << outputName << endl;
	
	TFile file( outputName.c_str(), "RECREATE");
	
// */
	fOutput->Write(); // write all objects merged by fOutput (used by PROOF)
	
	// publish protocoll
	ofstream prot("protocol_hist.txt", ios::app); //append
	if(prot.is_open()){
		prot << "--------------------------------------------------------------------------------";
		prot << endl;
		prot << "opened data file: " << openedDataFile << endl;
		prot << "finished at: " << asctime (timeinfo) << endl; 	// readable timestamp
		prot << "total counted: " << Nnum << endl;
		prot << "total from cutflow: " << nTotalEvents << endl;
		prot << "Ntest1: " << Ntest1 << endl;
		prot << "Ntest2: " << Ntest2 << endl;
		prot << "Ntest3: " << Ntest3 << endl;
		prot << "Ntest4: " << Ntest4 << endl;		
		prot << "Number of real photons: Ng: " << Ng_1 << endl;
		prot << "Number of fake photons: Ne: " << Ne_1 << endl;
		prot << "Fakerate f = Ng/(Ng+Ne) = " << f << endl;
		for(int ll=0;ll<sizeNcounts;ll++){
			prot << "Ncounts[" << ll << "]: " << Ncounts[ll] << "\t";
		}
		prot << endl;
		prot << "results written in: " << outputName << endl;
		prot << "comment (in code): " << comment << endl;
		prot.close();
		
		cout << "protocol published!" << endl;
	}
	
	// publish filename
	ofstream myFile("outputfiles.txt", ios::app); // append
	if(myFile.is_open()){
		myFile << outputName << endl;
		myFile.close();
		cout << "filename published!" << endl;
		
		
	}
	
	
} // terminate































