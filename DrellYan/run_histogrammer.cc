
//
//


void run_histogrammer(){
	
	// tree file
	
	
	string file="../../NTupel/DYJetsToLL_M-50_nTuple.root";
	//string file="../../NTupel/DY_FullSim.root";
	
	ifstream f(file.c_str());
	if(!f.good()) {
		cout << "No file found, just compile." << endl;
		TSelector::GetSelector("HistogramProducer.cc+");
		return;
	}

	cout << file.c_str() << endl;
	
	gSystem->Load("pluginTreeWriterTreeWriterAuto.so");
	
		
	// create chain
	TChain ch("TreeWriter/eventTree");
	
	// 
	ch.AddFile( file.c_str() );
		
	double start_time = time(NULL); // measure running time
	
	cout << "... processing .." << endl;
	ch.Process("Histogrammer.cc+" );
	
	
	cout << "Runtime ~" << 1.*( time(NULL) - start_time)/1 << " sec." << endl;
	
	
	
}
