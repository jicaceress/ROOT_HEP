#include <iostream>
#include <string>
#include <stdio.h>

void leadingRecoTruth(){
	TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
	TTree *tree = (TTree*) file->Get("JetRecoTree");
	
	float evtw = -1;
	vector<float> *reco_R10_pt=0;
	vector<float> *truth_R10_pt=0;
    vector<float> *reco_R10_trimmed_pt=0;
	vector<float> *truth_R10_trimmed_pt=0;

	tree->SetBranchAddress("EventWeight", &evtw);
	tree->SetBranchAddress("RecoJets_R10_pt", &reco_R10_pt);
	tree->SetBranchAddress("TruthJets_R10_pt", &truth_R10_pt);
	tree->SetBranchAddress("RecoJets_R10_Trimmed_pt", &reco_R10_trimmed_pt);
	tree->SetBranchAddress("TruthJets_R10_Trimmed_pt", &truth_R10_trimmed_pt);

    TH1F *hist_reco_R10_pt = new TH1F("Reco Jet","Leading Reco Jet pT;pT(GeV);Events",50,0,600);
    TH1F *hist_truth_R10_pt = new TH1F("Truth Jet","Leading Truth Jet pT;pT(GeV);Events",50,0,600);
    TH1F *hist_reco_R10_trimmed_pt = new TH1F("Reco trimmed Jet","Leading Reco Trimmed Jet pT;pT(GeV);Events",50,0,600);
    TH1F *hist_truth_R10_trimmed_pt = new TH1F("Truth trimmed Jet","Leading Truth Trimmed Jet pT;pT(GeV);Events",50,0,600);
    

	TCanvas *canvas_1 = new TCanvas("Canvas1","",1000,1000);
	TCanvas *canvas_2 = new TCanvas("Canvas2","",1000,1000);
    TCanvas *canvas_3 = new TCanvas("Canvas3","",1000,1000);

    int nentries, nbytes, i;
	nentries = (Int_t)tree->GetEntries();

	for (i = 0; i < nentries; i++){
        nbytes = tree->GetEntry(i);
        if(reco_R10_pt->size()!=0){
            hist_reco_R10_pt->Fill(reco_R10_pt->at(0)/1000.,evtw);
        }
        if(truth_R10_pt->size()!=0){
            hist_truth_R10_pt->Fill(truth_R10_pt->at(0)/1000.,evtw);
        }
        if(reco_R10_trimmed_pt->size()!=0){
            hist_reco_R10_trimmed_pt->Fill(reco_R10_trimmed_pt->at(0)/1000.,evtw);
        }
        if(truth_R10_trimmed_pt->size()!=0){
            hist_truth_R10_trimmed_pt->Fill(truth_R10_trimmed_pt->at(0)/1000.,evtw);
        }
    }

    canvas_1->cd();
    hist_reco_R10_pt->SetLineStyle(1);
    hist_reco_R10_pt->SetLineColor(1);
    hist_reco_R10_pt->Draw("");
    hist_reco_R10_trimmed_pt->SetLineStyle(10);
    hist_reco_R10_trimmed_pt->SetLineColor(2);
    hist_reco_R10_trimmed_pt->Draw("same");
    auto legend1 = new TLegend(0.3,0.7,0.48,0.85);
    legend1->AddEntry(hist_reco_R10_pt,"Leading Reco Jet pT - R10");
    legend1->AddEntry(hist_reco_R10_trimmed_pt,"Leading Reco Trimmed Jet pT");
    legend1->Draw();
    canvas_1->SetLogy();
    canvas_1->Print("recoR10JetpT.ps");

    canvas_2->cd();
    hist_truth_R10_pt->SetMarkerStyle(20);
    hist_truth_R10_pt->SetMarkerColor(kRed);
    hist_truth_R10_pt->Draw("");
    hist_truth_R10_trimmed_pt->SetMarkerStyle(21);
    hist_truth_R10_trimmed_pt->SetMarkerColor(1);
    hist_truth_R10_trimmed_pt->Draw("same");
    auto legend2 = new TLegend(0.3,0.7,0.48,0.85);
    legend2->AddEntry(hist_truth_R10_pt,"Leading Truth Jet pT - R10");
    legend2->AddEntry(hist_truth_R10_trimmed_pt,"Leading Truth Trimmed Jet pT");
    legend2->Draw();
    canvas_2->SetLogy();
    canvas_2->Print("truthR10JetpT.ps");

    canvas_3->cd();
    hist_reco_R10_pt->SetMarkerStyle(20);
    hist_reco_R10_pt->SetMarkerColor(1);
    hist_reco_R10_pt->Draw("");
    hist_reco_R10_trimmed_pt->SetMarkerStyle(21);
    hist_reco_R10_trimmed_pt->SetMarkerColor(2);
    hist_reco_R10_trimmed_pt->Draw("same");
    hist_truth_R10_pt->SetMarkerStyle(22);
    hist_truth_R10_pt->SetMarkerColor(4);
    hist_truth_R10_pt->Draw("same");
    hist_truth_R10_trimmed_pt->SetMarkerStyle(23);
    hist_truth_R10_trimmed_pt->SetMarkerColor(6);
    hist_truth_R10_trimmed_pt->Draw("same");
    auto legend3 = new TLegend(0.3,0.7,0.48,0.85);
    legend3->AddEntry(hist_truth_R10_pt,"Leading Truth Jet pT - R10");
    legend3->AddEntry(hist_truth_R10_trimmed_pt,"Leading Truth Trimmed Jet pT");
    legend3->AddEntry(hist_reco_R10_pt,"Leading Reco Jet pT - R10");
    legend3->AddEntry(hist_reco_R10_trimmed_pt,"Leading Reco Trimmed Jet pT");
    legend3->Draw();
    canvas_3->SetLogy();
    canvas_3->Print("allR10JetpT.ps");
}
