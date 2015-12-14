
//
//


void run_histogrammer(){
	
	// tree file
	
	
	string file="/user/rmeyer/mergedNTuples/DY_v1.root" // mc sample with particle flow candidates
	
	//string file="../../cms_sw/NTupel/DYJetsToLL_M-50_nTuple.root";
	//string file="../../NTupel/DY_FullSim.root";
	//string file="../../cms_sw/NTupel/myTestTree_unpacked_6_pt095charged.root";
	
	//string file="../../cms_sw/myTrees/mcDY_1446802738.root";
		
	const char* selector="Histogrammer.cc+";
	
	
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
	
	cout << "... processing .." << endl;
	ch.Process(selector);
	
	
	cout << "Runtime ~" << 1.*( time(NULL) - start_time)/1 << " sec." << endl;
	
	
	
	
}
