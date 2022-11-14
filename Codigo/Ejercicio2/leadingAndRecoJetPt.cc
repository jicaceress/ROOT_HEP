#include <iostream>
#include <string>
#include <stdio.h>

void leadingAndRecoJetPt(){
	TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
	TTree *tree = (TTree*) file->Get("JetRecoTree");
	
	float evtw = -1;
	vector<float> *reco_R4_pt = 0;
	vector<float> *truth_R4_pt = 0;

	tree->SetBranchAddress("EventWeight", &evtw);
	tree->SetBranchAddress("RecoJets_R4_pt", &reco_R4_pt);
	tree->SetBranchAddress("TruthJets_R4_pt", &truth_R4_pt);

	TCanvas *canvas_1 = new TCanvas("Canvas1","",1000,1000);
	TCanvas *canvas_2 = new TCanvas("Canvas2","",1000,1000);
    TCanvas *canvas_3 = new TCanvas("Canvas3","",1000,1000);
    TCanvas *canvas_4 = new TCanvas("Canvas4","",1000,1000);
    TCanvas *canvas_5 = new TCanvas("Canvas5","",1000,1000);
    TCanvas *canvas_6 = new TCanvas("Canvas6","",1000,1000);

	TH1F *hist_leadreco_weighted_pt = new TH1F("Lead Reco-jet weighted","Leading weighted jet pT; pT(GeV);Events",50,10,200);
	TH1F *hist_leadreco_pt = new TH1F("Lead Reco-jet","Leading jet pT;pT(GeV);Events",50,10,200);
	TH1F *hist_reco_weighted_pt = new TH1F("Reco-jet weighted","Weighted jet pT; pT(GeV);Events",50,10,200);
    TH1F *hist_reco_pt = new TH1F("Reco-jet","Jet pT; pT(GeV);Events",50,10,200);
	TH1F *hist_leadtruth_weighted_pt = new TH1F("Lead weighted Truth-jet","Leading weighted jet pT; pT(GeV);Events",50,10,200);
	TH1F *hist_truth_weighted_pt = new TH1F("Truth-jet weighted","Weighted jet pT;pT(GeV);Events",50,10,200);
    TH1F *hist_leadtruth_pt = new TH1F("Lead Truth-jet","Leading truth jet pT; pT(GeV);Events",50,10,200);
	TH1F *hist_truth_pt = new TH1F("Truth-jet","Truth jet pT;pT(GeV);Events",50,10,200);

	int nentries, nbytes, i;
	nentries = (Int_t)tree->GetEntries();

	for (i = 0; i < nentries; i++)
	{
    		nbytes = tree->GetEntry(i);
    		if(reco_R4_pt->size()>0){
                hist_leadreco_pt->Fill(reco_R4_pt->at(0)/1000.);
    			hist_leadreco_weighted_pt->Fill(reco_R4_pt->at(0)/1000.,evtw);
    			for(int j=0; j<reco_R4_pt->size(); j++){
                    hist_reco_pt->Fill(reco_R4_pt->at(j)/1000.);
        			hist_reco_weighted_pt->Fill(reco_R4_pt->at(j)/1000.,evtw);
    			}
    		}
            if(truth_R4_pt->size()>0){
                hist_leadtruth_pt->Fill(truth_R4_pt->at(0)/1000.);
    			hist_leadtruth_weighted_pt->Fill(truth_R4_pt->at(0)/1000.,evtw);
    			for(int j=0; j<truth_R4_pt->size(); j++){
                    hist_truth_pt->Fill(truth_R4_pt->at(j)/1000.);
        			hist_truth_weighted_pt->Fill(truth_R4_pt->at(j)/1000.,evtw);
    			}
    		}
	}
	canvas_1->cd();
	hist_reco_weighted_pt->SetMarkerStyle(20);
	hist_reco_weighted_pt->SetMarkerColor(kRed);
	hist_reco_weighted_pt->Draw("");
	hist_leadreco_weighted_pt->SetMarkerStyle(21);
	hist_leadreco_weighted_pt->Draw("same");
    auto legend1 = new TLegend(0.3,0.7,0.48,0.85);
    legend1->AddEntry(hist_reco_weighted_pt,"Jet pT","p");
    legend1->AddEntry(hist_leadreco_weighted_pt,"Leading Jet pT","p");
    legend1->Draw();
	
	canvas_1->SetLogy();
	canvas_1->Print("recoPt.ps");

    canvas_2->cd();
	hist_truth_weighted_pt->SetMarkerStyle(20);
	hist_truth_weighted_pt->SetMarkerColor(kRed);
	hist_truth_weighted_pt->Draw("");
	hist_leadtruth_weighted_pt->SetMarkerStyle(21);
	hist_leadtruth_weighted_pt->Draw("same");
    auto legend2 = new TLegend(0.3,0.7,0.48,0.85);
    legend2->AddEntry(hist_truth_weighted_pt,"Truth Jet pT","p");
    legend2->AddEntry(hist_leadtruth_weighted_pt,"Leading Truth Jet pT","p");
    legend2->Draw();	

	canvas_2->SetLogy();
	canvas_2->Print("truthPt.ps");

    canvas_3->cd();
    hist_leadreco_pt->SetLineColor(2);
	hist_leadreco_pt->Draw("");
	hist_leadtruth_pt->SetLineColor(4);
	hist_leadtruth_pt->Draw("same");
    auto legend3 = new TLegend(0.3,0.7,0.48,0.85);
    legend3->AddEntry(hist_leadreco_pt,"Leading Jet pT","l");
    legend3->AddEntry(hist_leadtruth_pt,"Leading Truth Jet pT","l");
    legend3->Draw();    
	
	canvas_3->SetLogy();
	canvas_3->Print("leadRecoAndTruthPt.ps");

    canvas_4->cd();
	hist_leadreco_weighted_pt->SetMarkerStyle(20);
	hist_leadreco_weighted_pt->SetMarkerColor(kRed);
	hist_leadreco_weighted_pt->Draw("");
	hist_leadtruth_weighted_pt->SetMarkerStyle(21);
	hist_leadtruth_weighted_pt->Draw("same");
    auto legend4 = new TLegend(0.3,0.7,0.48,0.85);
    legend4->AddEntry(hist_leadreco_weighted_pt,"Leading Weighted Jet pT","p");
    legend4->AddEntry(hist_leadtruth_weighted_pt,"Leading Truth Weighted Jet pT","p");
    legend4->Draw();	

	canvas_4->SetLogy();
	canvas_4->Print("leadRecoAndTruthWeightedPt.ps");

    canvas_5->cd();
    hist_reco_pt->SetLineColor(2);
	hist_reco_pt->Draw("");
    hist_truth_pt->SetLineColor(4);
	hist_truth_pt->Draw("same");
    auto legend5 = new TLegend(0.3,0.7,0.48,0.85);
    legend5->AddEntry(hist_reco_pt,"Jet pT","l");
    legend5->AddEntry(hist_truth_pt,"Truth Jet pT","l");
    legend5->Draw();
	
	canvas_5->SetLogy();
	canvas_5->Print("RecoAndTruthPt.ps");

    canvas_6->cd();
	hist_reco_weighted_pt->SetMarkerStyle(20);
	hist_reco_weighted_pt->SetMarkerColor(kRed);
	hist_reco_weighted_pt->Draw("");
	hist_truth_weighted_pt->SetMarkerStyle(21);
    hist_truth_weighted_pt->SetMarkerColor(1);
	hist_truth_weighted_pt->Draw("same");
    auto legend6 = new TLegend(0.3,0.7,0.48,0.85);
    legend6->AddEntry(hist_reco_weighted_pt,"Weighted Jet pT","p");
    legend6->AddEntry(hist_truth_weighted_pt,"Truth Weighted Jet pT","p");
    legend6->Draw();
	
	canvas_6->SetLogy();
	canvas_6->Print("RecoAndTruthWeightedPt.ps");

}
