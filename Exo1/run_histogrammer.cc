//
//


void run_histogrammer(string comment, bool allData = false){
	
	bool run = true;
	
	//bool allData = false;
	
	// tree file
	string path = "/user/rmeyer/mergedNTuples/";
	
	//string file= path+"DY_v1.root"; // mc sample with particle flow candidates
	//string file= path+"DYJetsToLL_M-50.root";
	
	int kk = 0;
	int Ndatasets = 2; // number of datasets
	
	string datasets[] = {
						/*
							"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v01.root",		// 0
							"TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_v01.root", // 1
							"WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v01.root",	// 2
							"ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v01.root", // 3
							"DoubleEG_Run2015D-05Oct2015-v1.root",	// data  4
							"DoubleEG_Run2015D-PromptReco-v4.root"	// data  5
						*/
							"DoubleEG_Run2015D-05Oct2015-v1_Exo_v2.root",
							"DoubleEG_Run2015D-PromptReco-v4_Exo_v2.root"
						};
	
	string file;
	
	cout << datasets[0] << endl;
	cout << datasets[1] << endl;
	
	//cout << "printed all" << endl;
	//cout << *(datasets+5) << endl;
	
	//string file="../../cms_sw/NTupel/DYJetsToLL_M-50_nTuple.root";
	//string file="../../NTupel/DY_FullSim.root";
	//string file="../../cms_sw/NTupel/myTestTree_unpacked_6_pt095charged.root";
	//string file="../../cms_sw/myTrees/mcDY_1446802738.root";
	
	if(run){		// if boolean is true
		
		const char* selector="Histogrammer.cc+";
		
		if(allData){		// analyse all data sets
		
			cout << "allData" << endl;
			for(; kk < Ndatasets; kk++){		//loop datasets
				cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
				cout << "++++++++++++++++++++++ "<<kk<<" ++++++++++++++++++++++++++++++++++++++" << endl;
				cout << datasets[kk] << endl;
				
				file = path + datasets[kk];
				
				ifstream f(file.c_str());
				if(!f.good()) {
					cout << "No file found, just compile." << endl;
					TSelector::GetSelector(selector);
					return;
				}
				
				cout << file.c_str() << endl;
				
				gSystem->Load("pluginTreeWriterTreeWriterAuto.so");
				cout << "successfull" << endl;
				
				// create chain
				TChain ch("TreeWriter/eventTree");
				// 
				ch.AddFile( file.c_str() );
				
				double start_time = time(NULL); // measure running time
				
				cout << "processing with file: " << file << endl;
				cout << "... processing .." << endl;
				
				ch.Process(selector);
				
				double end_time = 1.*( time(NULL));
				cout << "Runtime ~" << (end_time - start_time)/1 << " sec." << endl;
				
				
				ofstream prot("protocol_hist.txt", ios::app); //append
				if(prot.is_open()){
					prot << "comment (call): " << comment << endl;
					prot << "Runtime ~" << (end_time - start_time)/1 << " sec." << endl;
					prot.close();
				}
				
				
			} // for(NDatasets) loop datasets
		} // if allData
		
	}
	
	
} // void run_histogrammer.cc
