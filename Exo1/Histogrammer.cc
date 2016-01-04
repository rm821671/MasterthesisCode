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
const int sizeNcounts = 100;

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



// calculate parameters for corrected iso gamma
void fIsoGammaCorr(float eta, float& alpha, float& A, float& kappa){
	
	alpha = 2.5;		// always
	A = 0.;				// default
	kappa = 4.5e-3;		// only change in specific cases, see below
	
	// AN Table 8:
	if(eta < 0.9) A = 0.17;
	else if(eta > 0.9 && eta < 1.4442) A = 0.14;
	else if(eta > 1.566 && eta < 2.0) A = 0.11;
	else if(eta > 2.0 && eta < 2.2){
		A = 0.14;
		kappa = 3.0e-3;
	}
	else if(eta > 2.2 && eta < 2.5){
		A = 0.22;
		kappa = 3.0e-3;
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
		
		TTreeReaderValue<Float_t> rho;
		
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
		
		TTreeReaderValue<Bool_t> signalTrigger; /* this is the used trigger */
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
		
		float Etemp[sizeNcounts];
		float PSum[sizeNcounts];
		
		float isoGammaCorr;
		
		
		float 	Ng_1 =0,
				Ne_1 =0;
		float	Ne_pt[10],
				Ng_pt[10];
		
		bool isData;
		
		// some name strings
		string 	hname,
				h2name,
				numstr;
		
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
	rho ( fReader, "rho" ),
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
	
	
	hname = "diphoton_EBEB";
	h[hname] = TH1F((hname+"_"+s).c_str(),";m [GeV];counts",2000,0,2000);
	
	hname = "diphoton_EBEE";
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
		
		for(int i = 0; i < sizeNcounts; i++){
			vTemp[i].SetXYZ(0.,0.,0.);			// reset temp vectors
			lvTemp[i].SetPxPyPzE(0.,0.,0.,0.);	// reset temp lorentzvectors
			Etemp[i] = 0.;
			PSum[i] = 0.;
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
	// comment this out to analyze the whole tree:
	//if(Nnum > 1000000) return kTRUE;
	
	resetSelection();
	
	fReader.SetEntry(entry);	// fReader points on current entry
	
	// set weight
	selW = *mc_weight * *pu_weight;
	
	// counter for EBEB and EBEE
	int		nEBEB = 0,
			nEBEE = 0;
	
	
	
	// counter for the combinations:
	// maximum photon number is 14, thus 91 combinations are possible
	int 	comb = 0,
			selComb = 0; // counts selected combinations
	
	// trigger:
	if(*signalTrigger){
		
		// loop photons
		// and count kinematic criteria
		// maximum number of photons: 14
		for(int ip=0,np=photons.GetSize(); ip < np;ip++){
			
			// calculate all diphoton objects per event
			for(int iq=ip+1; iq<np; iq++){
				// check for criteria
				if(	(photons[ip].p.Pt() > 75. && photons[iq].p.Pt() > 75.) &&
					(fabs(photons[ip].p.Eta()) < 2.5 && fabs(photons[iq].p.Eta()) < 2.5) &&
					(fabs(photons[ip].p.Eta()) < 1.4442 || fabs(photons[iq].p.Eta()) < 1.4442) 
					){
					
					lvTemp[comb].SetVect(photons[ip].p + photons[iq].p);			// summed momentum
					lvTemp[comb].SetE(photons[ip].p.Mag() + photons[iq].p.Mag());	// summed energy
					
					// EBEB
					if( ((fabs(photons[ip].p.Eta()) < 1.4442 && fabs(photons[iq].p.Eta()) < 1.4442) &&
						lvTemp[comb].M() > 230.) ){
							nEBEB++;
							selPhotons.push_back(&photons[ip]);
							selPhotons2.push_back(&photons[iq]);
							
							PSum[selComb] = photons[ip].p.Mag() + photons[iq].p.Mag();
							
							selComb++;
							Ncounts[0] ++;
						
					}
					
					// EBEE
					if( ((fabs(photons[ip].p.Eta()) >= 1.4442 || fabs(photons[iq].p.Eta()) >= 1.4442) &&
						lvTemp[comb].M() > 320.) ){
							nEBEE++;
							
							selPhotons.push_back(&photons[ip]);
							selPhotons2.push_back(&photons[iq]);
							
							PSum[selComb] = photons[ip].p.Mag() + photons[iq].p.Mag();
							
							selComb++;
							Ncounts[1] ++;
							
						}
					
				}
				//comb++; // all combinations (no need to count this)
			}
		}
		
		if(nEBEB && nEBEE){
			Ncounts[2]++;
			// this is rougly 1% of all events, 
			// where more than one diphoton object fullfills the criteria
			
			//cout << Ncounts[2] << "\t nEBEB: " << nEBEB << " - nEBEE: " << nEBEE << endl;
		}
		
		// at least one diphoton object passed the selection:
		if(nEBEB>0 || nEBEE>0){
			
			bool	p1 = false,
					p2 = false;
			
			// the index of the combination with the highest scalar sum of momenta:
			int maxComb = distance(PSum, max_element(PSum, PSum+sizeNcounts));
			
			// the criteria have to be fullfilled for both photons individually
			float alpha, A, kappa;
			
			// selection cuts of AN section 6.1
			// photon 1
			fIsoGammaCorr(fabs(selPhotons[maxComb]->p.Eta()), alpha, A, kappa);
			isoGammaCorr = alpha + selPhotons[maxComb]->isoPhotonsEA - (*rho) * A - kappa * selPhotons[maxComb]->p.Pt();
			if(	(fabs(selPhotons[maxComb]->p.Eta()) < 1.4442 &&
				selPhotons[maxComb]->isoChargedHadronsEA > 5 &&
				isoGammaCorr > 2.75 &&
				selPhotons[maxComb]->hOverE > 5.0e-2 &&
				selPhotons[maxComb]->sigmaIetaIeta > 0.0105)
				||
				(fabs(selPhotons[maxComb]->p.Eta()) > 1.566 &&
				selPhotons[maxComb]->isoChargedHadronsEA > 5 &&
				isoGammaCorr > 2.0 &&
				selPhotons[maxComb]->hOverE > 5.0e-2 &&
				selPhotons[maxComb]->sigmaIetaIeta > 0.028)
				){
					p1 = true;
			}
			// photon 2
			fIsoGammaCorr(fabs(selPhotons2[maxComb]->p.Eta()), alpha, A, kappa);
			isoGammaCorr = alpha + selPhotons2[maxComb]->isoPhotonsEA - (*rho) * A - kappa * selPhotons2[maxComb]->p.Pt();
			if(	(fabs(selPhotons2[maxComb]->p.Eta()) < 1.4442 &&
				selPhotons2[maxComb]->isoChargedHadronsEA > 5 &&
				isoGammaCorr > 2.75 &&
				selPhotons2[maxComb]->hOverE > 5.0e-2 &&
				selPhotons2[maxComb]->sigmaIetaIeta > 0.0105)
				||
				(fabs(selPhotons2[maxComb]->p.Eta()) > 1.566 &&
				selPhotons2[maxComb]->isoChargedHadronsEA > 5 &&
				isoGammaCorr > 2.0 &&
				selPhotons2[maxComb]->hOverE > 5.0e-2 &&
				selPhotons2[maxComb]->sigmaIetaIeta > 0.028)
				){
					p2 = true;
			}
			
			// now fill (if both pass the cuts)
			if(p1 && p2){
				lvTemp[comb].SetVect(	selPhotons[maxComb]->p + 		selPhotons2[maxComb]->p);		// summed momentum
				lvTemp[comb].SetE(		selPhotons[maxComb]->p.Mag() + 	selPhotons2[maxComb]->p.Mag());	// summed energy
				
				// EBEB or EBEE?
				if(	fabs(selPhotons[maxComb]->p.Eta()) < 1.4442 && 
					fabs(selPhotons2[maxComb]->p.Eta()) < 1.4442 ){ //EBEB
						hname = "diphoton_EBEB";
						h[hname].Fill(lvTemp[comb].M());
					}
					
				if(	fabs(selPhotons[maxComb]->p.Eta()) >= 1.4442 || 
					fabs(selPhotons2[maxComb]->p.Eta()) >= 1.4442 ){ //EBEE
						hname = "diphoton_EBEE";
						h[hname].Fill(lvTemp[comb].M());
					}
				
				
				
			}
		}
		
		
		
		
		
		

		
		
		
	
	} // if signalTrigger
	
	
	
	
	
	
	
	
	
	
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































