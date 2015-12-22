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

// */
// namespace stuff
using namespace tree;
using namespace std;

/*******************************************************************************************
 * 
 * own functions and stuff
 *
 *******************************************************************************************/

// number of test count variables
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
    outputName = "exo_selector_results_" + 
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
		
		int nTotalEvents = 0; // weighted number of monte carlo events
		int nTotalEvents_bin = 2; // bin of the cutflow diagram
		Float_t selW = 1.;	// weight
		int Nnum = 0; // total counter, added up in every event
		int Ncounts[sizeNcounts] = {}; // some test counter int
		TVector3 	vTemp[sizeNcounts]; // array of vectors
		TLorentzVector	lvTemp[sizeNcounts]; // array of lorentz vectors
		
		float 	Ng_1 =0,
				Ne_1 =0;
		float	Ne_pt[10],
				Ng_pt[10];
		
		bool isData;
		
		map<string,TEfficiency> eff;	// efficiencies
		map<string,TH1F> h;				// 1d histograms
		map<string,TGraph> g;			// graphs , NOT USED
		map<string,TH2F> h2;			// 2d histograms
		
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
	
	
	signalTrigger( fReader, "HLT_DoublePhoton60_v" ), /* this is the trigger used for the analysis */
	
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
	
	string hname;
	string h2name;
	
	hname = "diphoton_excess";
	h[hname] = TH1F((hname+"_"+s).c_str(),";m [GeV];counts",2000,0,2000);
	
	// TH2F (const char *name, const char *title, 
	//	Int_t nbinsx, Double_t xlow, Double_t xup, 
	//	Int_t nbinsy, Double_t ylow, Double_t yup)
	h2name = "Ntrk_Pt_plane";
	h2[h2name] = TH2F((h2name + "_" + s).c_str(),";Number;Pt",
						300,0,300,
						100,0,10);
	
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
		nTotalEvents = hcut->GetBinContent(nTotalEvents_bin);
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
		selPhotons2.clear();
		selJets2.clear();
		selBJets2.clear();
		selElectrons2.clear();
		selMuons2.clear();	
		
		//TagElectron = false;
		for(int ii = 0; ii < sizeNcounts; ii++){
			vTemp[ii].SetXYZ(0.,0.,0.);			// reset temp vectors
			lvTemp[ii].SetPxPyPzE(0.,0.,0.,0.);	// reset temp lorentzvectors
		
		}
		
		
		
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
	
	resetSelection();	
	fReader.SetEntry(entry);	// fReader points on current entry
	
	// set weight
	selW = *mc_weight * *pu_weight;
	
	
	
	
	
	
	
	
	
	
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































