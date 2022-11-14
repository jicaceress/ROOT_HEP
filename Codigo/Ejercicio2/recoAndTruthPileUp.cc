#include <iostream>
#include <string>
#include <stdio.h>

void recoAndTruthPileUp(){
    TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
	TTree *tree = (TTree*) file->Get("JetRecoTree");
    float evtw = -1;
    UInt_t npv = -1; 
    float mu = -1;
    vector<float> *reco_R4_pt = 0;
	vector<float> *truth_R4_pt = 0;
    
    tree->SetBranchAddress("EventWeight", &evtw);
    tree->SetBranchAddress("NPV", &npv);
    tree->SetBranchAddress("mu_average", &mu);
    tree->SetBranchAddress("RecoJets_R4_pt", &reco_R4_pt);
	tree->SetBranchAddress("TruthJets_R4_pt", &truth_R4_pt);

    TCanvas *canvas_1 = new TCanvas("Canvas1","",1000,1000);
    TCanvas *canvas_2 = new TCanvas("Canvas2","",1000,1000);
    TCanvas *canvas_3 = new TCanvas("Canvas3","",1000,1000);
    TCanvas *canvas_4 = new TCanvas("Canvas4","",1000,1000);

    TH2F *hist_jetpt_npv = new TH2F("Reco-jet pT vs. NPV",";NPV; jet pT",50,1,50, 20, 0, 200);
    TH2F *hist_truth_jetpt_npv = new TH2F("Truth-jet pT vs. NPV",";NPV; jet pT",50,1,50, 20, 0, 200);
    TProfile *prof_jetpt_npv = new TProfile("Profile Reco-jet pT vs. NPV",";NPV; jet pT",50,1,50, 0, 200);
    TProfile *prof_truth_jetpt_npv = new TProfile("Profile Truth-jet pT vs. NPV",";NPV; jet pT",50,1,50, 0, 200);

    int nentries, nbytes, i;
    nentries = (Int_t)tree->GetEntries();

    for (i = 0; i < nentries; i++)
    {
        nbytes = tree->GetEntry(i);

        if(reco_R4_pt->size()!=0 && reco_R4_pt->at(0)>20000.){
            for(int j=0; j<reco_R4_pt->size(); j++){
                hist_jetpt_npv->Fill(reco_R4_pt->at(j)/1000.,npv,evtw);
                prof_jetpt_npv->Fill(reco_R4_pt->at(j)/1000.,npv,evtw);
            }
        }
        if(truth_R4_pt->size()!=0 && truth_R4_pt->at(0)>20000.){
            for(int j=0; j<truth_R4_pt->size(); j++){
                hist_truth_jetpt_npv->Fill(truth_R4_pt->at(j)/1000.,npv,evtw);
                prof_truth_jetpt_npv->Fill(truth_R4_pt->at(j)/1000.,npv,evtw);
            }
        }
    }

    canvas_1->cd();
    hist_jetpt_npv->Draw("colz");
    canvas_1->Print("colorMapRecoPileUp.ps");

    canvas_2->cd();
    hist_truth_jetpt_npv->Draw("colz");
    canvas_2->Print("colorMapTruthPileUp.ps");

    canvas_3->cd();
    prof_jetpt_npv->Draw("");
    canvas_3->Print("recoPileUp.ps");

    canvas_4->cd();
    prof_truth_jetpt_npv->Draw("");
    canvas_4->Print("truthPileUp.ps");
}
