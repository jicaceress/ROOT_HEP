#include <iostream>
#include <string>
#include <stdio.h>

void recoAndTruthJVF(){
    TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
	TTree *tree = (TTree*) file->Get("JetRecoTree");
    float evtw = -1;
    UInt_t npv = -1; 
    float mu = -1;
    vector<float> *reco_R4_pt = 0;
	vector<float> *truth_R4_pt = 0;
    vector<float> *reco_R4_jvf = 0;
    vector<float> *reco_tracks_R4_pt = 0;
    
    tree->SetBranchAddress("EventWeight", &evtw);
    tree->SetBranchAddress("NPV", &npv);
    tree->SetBranchAddress("mu_average", &mu);
    tree->SetBranchAddress("RecoJets_R4_pt", &reco_R4_pt);
	tree->SetBranchAddress("TruthJets_R4_pt", &truth_R4_pt);
    tree->SetBranchAddress("TrackJets_R4_pt", &reco_tracks_R4_pt);
    tree->SetBranchAddress("RecoJets_R4_jvf", &reco_R4_jvf);

    TH1F *hist_leadreco_jvf = new TH1F("Lead Reco-jet JVF","Leading jet JVF; JVF;Events",20,-1,1);
    TH1F *hist_leadreco100_jvf = new TH1F("Lead Reco-jet JVF","Leading jet JVF; JVF;Events",20,-1,1);
    TH1F *hist_leadreco_pt = new TH1F("Lead Reco-jet","Leading jet pT;pT(GeV);Events",50,10,200);
    TH1F *hist_leadtruth_pt = new TH1F("Lead Truth-jet","Leading truth jet pT; pT(GeV);Events",50,10,200);
    TH1F *hist_reco_pt = new TH1F("Reco-jet","Jet pT; pT(GeV);Events",50,10,200);
    TH1F *hist_reco_tracks_pt = new TH1F("Lead reco track","Leading truth track pT;pT(GeV);Events",50,10,200);

    TH2F *hist_reco_trackpt_npv = new TH2F("Track pT vs. NPV",";NPV; Track pT",50,1,50, 20, 0, 200);
    TProfile *prof_trackpt_npv = new TProfile("Profile Reco-track pT vs. NPV",";NPV; track pT",50,1,50, 0, 200);

    TCanvas *canvas_1 = new TCanvas("Canvas1","",1000,1000);
    TCanvas *canvas_2 = new TCanvas("Canvas2","",1000,1000);
    TCanvas *canvas_3 = new TCanvas("Canvas3","",1000,1000);
    TCanvas *canvas_4 = new TCanvas("Canvas4","",1000,1000);
    TCanvas *canvas_5 = new TCanvas("Canvas5","",1000,1000);
    

    int nentries, nbytes, i;
    nentries = (Int_t)tree->GetEntries();

    for (i = 0; i < nentries; i++)
    {
        nbytes = tree->GetEntry(i);

        if(reco_R4_pt->size()!=0 && reco_R4_pt->at(0)>20000.){
            hist_leadreco_jvf->Fill(reco_R4_jvf->at(0), evtw);
            if(reco_R4_pt->at(0)>100000.){hist_leadreco100_jvf->Fill(reco_R4_jvf->at(0), evtw);}
        }
        if(reco_tracks_R4_pt->size()!=0 && reco_tracks_R4_pt->at(0)>20000.){
            for(int j=0; j<reco_tracks_R4_pt->size(); j++){
                hist_reco_trackpt_npv->Fill(reco_tracks_R4_pt->at(j)/1000.,npv,evtw);
                prof_trackpt_npv->Fill(reco_tracks_R4_pt->at(j)/1000.,npv,evtw);
            }
        }
        if(reco_R4_pt->size()>0){
			hist_leadreco_pt->Fill(reco_R4_pt->at(0)/1000.,evtw);
            if(fabs(reco_R4_jvf->at(0))>0.5){
    			hist_reco_pt->Fill(reco_R4_pt->at(0)/1000.,evtw);
            }
		}
        if(truth_R4_pt->size()>0){
			hist_leadtruth_pt->Fill(truth_R4_pt->at(0)/1000.,evtw);
		}
        if(reco_tracks_R4_pt->size()>0){
            hist_reco_tracks_pt->Fill(reco_tracks_R4_pt->at(0)/1000.,evtw);
        }
    }
    
    canvas_1->cd();
    hist_leadreco_jvf->SetMarkerStyle(20);
    hist_leadreco_jvf->SetMarkerColor(kRed);
    hist_leadreco_jvf->DrawNormalized("");
    hist_leadreco100_jvf->SetMarkerStyle(22);
    hist_leadreco100_jvf->SetMarkerColor(kBlue);
    hist_leadreco100_jvf->DrawNormalized("same");
    canvas_1->Print("recoJVF.ps");

    canvas_2->cd();
    hist_leadtruth_pt->SetMarkerStyle(20);
    hist_leadtruth_pt->SetMarkerColor(kRed);
    hist_leadtruth_pt->Draw("");
    hist_leadreco_pt->SetMarkerStyle(21);
    hist_leadreco_pt->SetMarkerColor(kBlue);
    hist_leadreco_pt->Draw("same");
    hist_reco_pt->SetMarkerStyle(22);
    hist_reco_pt->SetMarkerColor(1);
    hist_reco_pt->Draw("same");
    auto legend2 = new TLegend(0.3,0.7,0.48,0.85);
    legend2->AddEntry(hist_leadtruth_pt,"Leading Truth Jet pT");
    legend2->AddEntry(hist_leadreco_pt,"Leading Reco Jet pT");
    legend2->AddEntry(hist_reco_pt,"Reco Jet pT");
    legend2->Draw();
    canvas_2->SetLogy();
    canvas_2->Print("comparisonJVF.ps");

    canvas_3->cd();
    hist_leadtruth_pt->SetMarkerStyle(20);
    hist_leadtruth_pt->SetMarkerColor(kRed);
    hist_leadtruth_pt->Draw("");
    hist_leadreco_pt->SetMarkerStyle(21);
    hist_leadreco_pt->SetMarkerColor(kBlue);
    hist_leadreco_pt->Draw("same");
    hist_reco_tracks_pt->SetMarkerStyle(22);
    hist_reco_tracks_pt->SetMarkerColor(6);
    hist_reco_tracks_pt->Draw("same");
    hist_reco_pt->SetMarkerStyle(23);
    hist_reco_pt->SetMarkerColor(1);
    hist_reco_pt->Draw("same");
    auto legend3 = new TLegend(0.3,0.7,0.48,0.85);
    legend3->AddEntry(hist_leadtruth_pt,"Leading Truth Jet pT");
    legend3->AddEntry(hist_leadreco_pt,"Leading Reco Jet pT");
    legend3->AddEntry(hist_reco_pt,"Reco Jet pT");
    legend3->AddEntry(hist_reco_tracks_pt,"Leading Reco Tracks pT");
    legend3->Draw();
    
    canvas_3->SetLogy();
    canvas_3->Print("comparisonTracks.ps");

    canvas_4->cd();
    hist_reco_trackpt_npv->Draw("colz");
    canvas_4->Print("colorMapTrackPileUp.ps");

    canvas_5->cd();
    prof_trackpt_npv->Draw("");
    canvas_5->Print("trackPileUp.ps");
}
