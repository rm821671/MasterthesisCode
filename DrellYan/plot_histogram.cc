
/*
++ root macro to create plots from histograms
++
++
++
 */


// breit wigner function
double func_bw(double *x, double *par){
	// par[0]: sigma0
	// par[1]: s
	// par[2]: Gamma
	return 		par[0]*1. / 
			(
				(par[1]-x[0])*(par[1]-x[0]) + 
				(x[0]*par[2])*(x[0]*par[2])  
			);
}


// fit breit wigner on histogram
int fit_tool(TH1F *h, TCanvas *ct){
	
	

	// read in histograms
	// TF1("f", fit function, range min, range max, number of variables)
	TF1 *f = new TF1("f",func_bw,80.,100.,3);
	
	h->Scale(1/h->GetEntries());

	f->SetParameter(1,6000);
	f->SetParameter(2,2.495);
	f->SetLineColor(kBlack);
	f->SetLineWidth(2);
	
	
		
	
	//cout << h->GetName() << endl;

	//TCanvas *ct = new TCanvas("ct","",700,600);
	ct->cd(3);	

	h->Fit(f);//,"MBR");
	h->Draw("hist");
	f->Draw("same");
	
	return 0;
	
	
} // fit_tool()






/*******************************************************************************************
 * main function
 *
+*/
void plot_histogram(){

	// some style options
	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(111111);
	gStyle->SetPalette(1);
	//TGaxis::SetMaxDigits(3);


	TGaxis::SetMaxDigits(3); // digits per axis, meaning xxxEy e.g. 1.00e03
	
	string infile="my_selector_results_delR01.root";
	//string infile="my_selector_results_DY_FullSim.root";
	


	TFile *f = new TFile(infile.c_str());

	
	/*
	TList *l = f->GetList();

	//string name;
	cout << "here: " << endl;
	cout << l->GetSize() << endl;
	cout << l->At(0) << endl;	

	for(int i=0;i < l->GetSize();i++){
		TObject *obj = (TObject*)l->At(i);
		cout << obj->ClassName() << endl;
		
		//
		if(strncmp(obj->ClassName(),"TH1F",4)==0){
			TH1F *h = (TH1*)obj;
			name = h->GetName();
			cout << name << endl;
		}
		else if(obj->ClassName() == "TH2"){
			TH2 *h = (TH2*)obj;
			name = h->GetName();
		}
		
	} // */

	// ATTENTION:	the saved histograms may have confusing names!
	//		check out the commented initialization in Histogrammer.cc 
/*		
	TH1F *h1 = (TH1F*)f->Get("m_gg_test");
	TH1F *h2 = (TH1F*)f->Get("m_eg_test");
	TH1F *h3 = (TH1F*)f->Get("met_all_test");
	TH1F *h4 = (TH1F*)f->Get("met_eg_test");
	TH1F *h5 = (TH1F*)f->Get("m_gg_all_test");
	TH1F *h6 = (TH1F*)f->Get("m_eg_all_test");	
	TH1F *h7 = (TH1F*)f->Get("met_eg_all_test");
	TH1F *h8 = (TH1F*)f->Get("pt_eg_test");
	TH1F *h9 = (TH1F*)f->Get("Nvtx_test");

	// 
	TH1F *h10 = (TH1F*)f->Get("m_eg_ptcut25_test");	// all probes, "eg" = eg + ee
	TH1F *h11 = (TH1F*)f->Get("m_gg_ptcut25_test");	// photon candidate probes, "gg" = eg
	
	// 
	TH1F *h12 = (TH1F*)f->Get("Htjets_eg_ptcut25_test");	// all probes, "eg" = eg + ee
	TH1F *h13 = (TH1F*)f->Get("Htjets_gg_ptcut25_test");	// photon candidate probes, "gg" = eg
	
	TH1F *h14 = (TH1F*)f->Get("m_ptcut25_0_den_test");	// see histogrammer.cc
	TH1F *h15 = (TH1F*)f->Get("m_ptcut25_0_num_test");	// see histogrammer.cc
	
	TH1F *h16 = (TH1F*)f->Get("m_ptcut25_1_den_test");	// see histogrammer.cc
	TH1F *h17 = (TH1F*)f->Get("m_ptcut25_1_num_test");	// see histogrammer.cc
	
	TH1F *h18 = (TH1F*)f->Get("number_photons_ptcut25_test");
	TH1F *h19 = (TH1F*)f->Get("number_photons_no_ptcut25_test");	
*/	
	TH1F* h20 = (TH1F*)f->Get("f_Pt_total_test");	
	TH1F* h21 = (TH1F*)f->Get("f_Pt_passed_test");	

	TH1F* h22 = (TH1F*)f->Get("f_Nvtx_total_test");	
	TH1F* h23 = (TH1F*)f->Get("f_Nvtx_passed_test");

	TH1F* h24 = (TH1F*)f->Get("f_JetSize_total_test");	
	TH1F* h25 = (TH1F*)f->Get("f_JetSize_passed_test");
	
	TH1F* h26 = (TH1F*)f->Get("f_Eta_total_test");	
	TH1F* h27 = (TH1F*)f->Get("f_Eta_passed_test");

	TH1F* h28 = (TH1F*)f->Get("f_met_total_test");	
	TH1F* h29 = (TH1F*)f->Get("f_met_passed_test");
	
	TH1F* h30 = (TH1F*)f->Get("f_met_total_test");	
	TH1F* h31 = (TH1F*)f->Get("f_met_passed_test");
	
	TH1F* h32 = (TH1F*)f->Get("delR_test");
	TH2F* h33 = (TH2F*)f->Get("delR_plane_test");
	
	

	
	

	
	//h["f_Pt_total"].Fill(p.p.Pt());
	
	//****************************************************************************************************************
	//+
	//+	efficiencies
	//+
	//*****************************************************************************************************************
	TEfficiency* eff1 =0;
	// TEfficiency(const TH1& passed,const TH1& total) 
	if(TEfficiency::CheckConsistency(*h21, *h20)){
		eff1 = new TEfficiency(*h21, *h20);
		
	}
	
	TEfficiency* eff2 =0;
	// TEfficiency(const TH1& passed,const TH1& total) 
	if(TEfficiency::CheckConsistency(*h23, *h22)){
		eff2 = new TEfficiency(*h23, *h22);
		
	}

	TEfficiency* eff3 =0;
	// TEfficiency(const TH1& passed,const TH1& total) 
	if(TEfficiency::CheckConsistency(*h25, *h24)){
		eff3 = new TEfficiency(*h25, *h24);
		
	}
	
	TEfficiency* eff4 =0;
	// TEfficiency(const TH1& passed,const TH1& total) 
	if(TEfficiency::CheckConsistency(*h27, *h26)){
		eff4 = new TEfficiency(*h27, *h26);
		
	}

	TEfficiency* eff5 =0;
	// TEfficiency(const TH1& passed,const TH1& total) 
	if(TEfficiency::CheckConsistency(*h29, *h28)){
		eff5 = new TEfficiency(*h29, *h28);
		
	}
	
	TEfficiency* eff6 =0;
	// TEfficiency(const TH1& passed,const TH1& total) 
	if(TEfficiency::CheckConsistency(*h31, *h30)){
		eff6 = new TEfficiency(*h31, *h30);
		
	}

	// 
	// 
	//TF1 


	// TCanvas *x = new TCanvas("x","title",width,height);
	// TCanvas* c3 = new TCanvas("c3","",1100,500);	

	// ************************************************************************************
	// c1
	/*
	TCanvas* c1 = new TCanvas("c1","",1700,1100);

	// TLegend(x1,y1,x2,y2)
	leg11 = new TLegend(0.75-0.1,0.8-0.1,0.9-0.1,0.9-0.1);
	//leg->SetHeader("The Legend Title");
	
	leg11->AddEntry(h1,"e#gamma","l");
	leg11->AddEntry(h2,"#gamma#gamma","l");
	//leg->AddEntry("gr","Graph with error bars","lep");
	
	h1->SetLineColor(kRed);
	h1->SetLineWidth(2);
	h2->SetLineWidth(2);
	
	h3->SetLineColor(kRed);
	h3->SetLineWidth(2);
	h4->SetLineWidth(2);

	h5->SetLineColor(kRed);
	h5->SetLineWidth(2);
	h6->SetLineWidth(2);

	h7->SetLineWidth(2);
	h8->SetLineWidth(2);
	h9->SetLineWidth(2);

	
	c1->Divide(3,2);
	
	c1->cd(1);	
	h1->DrawClone("hist");
	h2->DrawClone("same hist");
	leg11->DrawClone("same");
	
	c1->cd(2);
	h3->Draw("hist");	
	h4->DrawClone("same hist");

	c1->cd(3);
	h5->DrawClone("hist");
	h6->DrawClone("same hist");
	
	c1->cd(4);
	h3->DrawClone("hist");
	h7->DrawClone("same hist");
	
	c1->cd(5);
	h8->DrawClone("hist");	


	c1->cd(6);
	h9->DrawClone("hist");

	c1->SaveAs("c1.root");
	//c1->Close();
	
	// */
	
	
	// ************************************************************************************
	// c2
	
	TCanvas* c2 = new TCanvas("c2","",1700,1100);
/*	
	h10->SetLineColor(kRed);
	h11->SetLineColor(kBlue);
	h10->SetLineWidth(2);
	h11->SetLineWidth(2);
	
	h14->SetLineWidth(2);
	h15->SetLineWidth(2);
	h14->SetLineColor(kRed);
	h15->SetLineColor(kBlue);
		
	h12->SetLineColor(kRed);
	h13->SetLineColor(kBlue);
	h12->SetLineWidth(2);
	h13->SetLineWidth(2);

	h16->SetLineColor(kRed);
	h17->SetLineColor(kBlue);
	h16->SetLineWidth(2);
	h17->SetLineWidth(2);

	h18->SetLineColor(kRed); // cut
	h19->SetLineColor(kBlue); // nocut
	h18->SetLineWidth(2);
	h19->SetLineWidth(2);
*/	
	h20->SetLineColor(kRed); 
	h21->SetLineColor(kBlue); 
	h20->SetLineWidth(2);
	h21->SetLineWidth(2);

	h22->SetLineColor(kRed); 
	h23->SetLineColor(kBlue); 
	h22->SetLineWidth(2);
	h23->SetLineWidth(2);

	h24->SetLineColor(kRed); 
	h25->SetLineColor(kBlue); 
	h24->SetLineWidth(2);
	h25->SetLineWidth(2);
	
	//TLegend(x1,y1,x2,y2)
	leg21 = new TLegend(0.75-0.1,0.8-0.1,0.9-0.1,0.9-0.1);
	//leg->SetHeader("The Legend Title");
	
	leg21->AddEntry(h20,"N_{g} + N_{e}","l");
	leg21->AddEntry(h21,"N_{g}","l");
	
	// ++++++++++++++++++++++++++++++
	// +	fit efficiencies
	
	// Pt
	//TF1* feff1 = new TF1("feff1","[0]+[1]*TMath::Power((x/[2]+1),-[3])",15,25);
	TF1* feff1 = new TF1("feff1","[0]+[1]/x+[2]/(x*x)",30,200);
	
	// Nvtx
	TF1* feff2 = new TF1("feff2","[0]+[1]*x",0,30);

	// Njets
	TF1* feff3 = new TF1("feff3","[0]+[1]*x",0,25);
	
	//  TText TText(Double_t x, Double_t y, const char* text)
	TText* t1 = new TText(100.,0.05,"p0+p1/x+p2/x^2");
	t1->SetTextAlign(22);
	t1->SetTextFont(43);
	t1->SetTextSize(20);
	t1->SetTextColor(kBlack);
	
	TText* t2 = new TText(15.,0.03,"p0+p1*x");
	t2->SetTextAlign(22);
	t2->SetTextFont(43);
	t2->SetTextSize(20);
	t2->SetTextColor(kBlack);
	
	TText* t3 = new TText(15.,0.03,"p0+p1*x");
	t3->SetTextAlign(22);
	t3->SetTextFont(43);
	t3->SetTextSize(20);
	t3->SetTextColor(kBlack);
	
	
	TGraphAsymmErrors* geff1 = eff1->CreateGraph();
	TGraphAsymmErrors* geff2 = eff2->CreateGraph();
	TGraphAsymmErrors* geff3 = eff3->CreateGraph();
	TGraphAsymmErrors* geff4 = eff4->CreateGraph();
	TGraphAsymmErrors* geff5 = eff5->CreateGraph();
	TGraphAsymmErrors* geff6 = eff6->CreateGraph();
		
	geff1->GetYaxis()->SetTitle("fakerate");
	geff2->GetYaxis()->SetTitle("fakerate");
	geff3->GetYaxis()->SetTitle("fakerate");
	geff4->GetYaxis()->SetTitle("fakerate");
	geff5->GetYaxis()->SetTitle("fakerate");
	geff6->GetYaxis()->SetTitle("fakerate");

	
	//geff1->SetMinimum(0.);
	geff1->SetMaximum(0.1); // y axis range maximum
	
	geff2->SetMinimum(0.015);
	geff2->SetMaximum(0.035); // y axis range maximum
	
	geff3->SetMinimum(0.015);
	geff3->SetMaximum(0.035); // y axis range maximum
	
	
	//feff1->SetParameter(2,14);
	//feff1->SetParameter(3,4.9);
	
	// 1 percent threshold
	//TF1* f1pc = new TF1("f1pc","0.01"
	
	feff1->SetLineColor(kBlue);
	feff2->SetLineColor(kBlue);
	feff3->SetLineColor(kBlue);
	
	geff1->Fit(feff1,"R");
	geff2->Fit(feff2,"R");
	geff3->Fit(feff3,"R");
	
	
	geff1->GetListOfFunctions()->AddFirst(feff1);

		
	c2->Divide(3,2);
	
	
	c2->cd(1);
	h20->DrawClone("hist");
	h21->DrawClone("same hist");
	leg21->DrawClone("same");
	
	c2->cd(4);
	geff1->DrawClone("AP*");
	
	feff1->GetXaxis()->SetRangeUser(30,200);
	feff1->DrawClone("same");
	
	t1->Draw("same");
		
	c2->cd(2);
	h22->DrawClone("hist");
	h23->DrawClone("same hist");
	
	c2->cd(5);
	geff2->DrawClone("AP*");
	feff2->GetXaxis()->SetRangeUser(0,30);
	feff2->DrawClone("same");
	t2->Draw("same");
	
	//feff2->DrawClone("same");
	
	c2->cd(3);
	h24->DrawClone("hist");
	h25->DrawClone("same hist");
	
	c2->cd(6);
	geff3->DrawClone("AP*");
	feff3->GetXaxis()->SetRangeUser(0,25);
	feff3->DrawClone("same");
	t3->Draw("same");
	
	c2->SaveAs("c2.root");
	//c2->Close();
	
	// */

	// ************************************************************************************
	// c3
	TCanvas* c3 = new TCanvas("c3","",1700,1100);
	
	h26->SetLineColor(kRed); 
	h27->SetLineColor(kBlue); 
	h26->SetLineWidth(2);
	h27->SetLineWidth(2);

	h28->SetLineColor(kRed); 
	h29->SetLineColor(kBlue); 
	h28->SetLineWidth(2);
	h29->SetLineWidth(2);


	//TLegend(x1,y1,x2,y2)
	leg31 = new TLegend(0.75-0.1,0.8-0.1,0.9-0.1,0.9-0.1);
	//leg->SetHeader("The Legend Title");
	
	leg31->AddEntry(h26,"N_{g} + N_{e}","l");
	leg31->AddEntry(h27,"N_{g}","l");
		
	c3->Divide(3,2);
	
	
	c3->cd(1);
	h26->DrawClone("hist");
	h27->DrawClone("same hist");
	leg31->DrawClone("same");
	
	c3->cd(4);
	geff4->DrawClone("AP*");		

	c3->cd(2);
	h28->DrawClone("hist");
	h29->DrawClone("same hist");
	
	c3->cd(5);
	geff5->DrawClone("AP*");
	
	c3->cd(3);
	
	
	c3->SaveAs("c3.root");

	//
	// 
	//fit_tool(h11,c2);
	// */


	// ************************************************************************************
	// c4
	TCanvas* c4 = new TCanvas("c4","",600,600);

	//c4->Divide(1,2);
	geff1->SetMaximum(0.04); // y axis range maximum
	c4->cd();
	geff1->DrawClone("AP*");	
	feff1->GetXaxis()->SetRangeUser(30,200);
	feff1->DrawClone("same");
	t1->Draw("same");
	
	// ************************************************************************************
	// c5
	TCanvas* c5 = new TCanvas("c5","",600,600);
		
	c5->cd();
	geff2->DrawClone("AP*");
	feff2->GetXaxis()->SetRangeUser(0,30);
	feff2->DrawClone("same");
	t2->Draw("same");

	// ************************************************************************************
	// c6
	TCanvas* c6 = new TCanvas("c6","",1700,1100);
	c6->Divide(3,2);
	
	c6->cd(1);	
	h32->DrawClone("hist");
	
	c6->cd(2);
	h33->DrawClone();



	// */	
	 
}	// plot_histogram()

















