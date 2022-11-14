#include <iostream>
#include <string>
#include <stdio.h>

void deltaRJets(){
    TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
	TTree *tree = (TTree*) file->Get("JetRecoTree");
    float evtw = -1;
    vector<float> *reco_R4_pt = 0;
	vector<float> *truth_R4_pt = 0;
    vector<float> *reco_R4_jvf = 0;
    vector<float> *reco_tracks_R4_pt = 0;
    vector<float> *reco_tracks_R4_eta = 0;
    vector<float> *reco_tracks_R4_phi = 0;
    vector<float> *reco_tracks_R4_m = 0;
    vector<float> *reco_R4_eta = 0;
    vector<float> *reco_R4_phi = 0;
    vector<float> *reco_R4_m = 0;
    vector<float> *truth_R4_eta = 0;
    vector<float> *truth_R4_phi = 0;
    vector<float> *truth_R4_m = 0;

    tree->SetBranchAddress("EventWeight", &evtw);
    tree->SetBranchAddress("RecoJets_R4_pt", &reco_R4_pt);
	tree->SetBranchAddress("TruthJets_R4_pt", &truth_R4_pt);
    tree->SetBranchAddress("TrackJets_R4_pt", &reco_tracks_R4_pt);
    tree->SetBranchAddress("RecoJets_R4_jvf", &reco_R4_jvf);
    tree->SetBranchAddress("TrackJets_R4_eta", &reco_tracks_R4_eta);
    tree->SetBranchAddress("TrackJets_R4_phi", &reco_tracks_R4_phi);
    tree->SetBranchAddress("TrackJets_R4_m", &reco_tracks_R4_m);
    tree->SetBranchAddress("RecoJets_R4_eta", &reco_R4_eta);
    tree->SetBranchAddress("RecoJets_R4_phi", &reco_R4_phi);
    tree->SetBranchAddress("RecoJets_R4_m", &reco_R4_m);
    tree->SetBranchAddress("TruthJets_R4_eta", &truth_R4_eta);
    tree->SetBranchAddress("TruthJets_R4_phi", &truth_R4_phi);
    tree->SetBranchAddress("TruthJets_R4_m", &truth_R4_m);

    TH1F *hist_DR_reco_truth = new TH1F("Delta R reco","Delta R; #Delta R; Events",20,0,2);
    TH1F *hist_DR_recojvf_truth = new TH1F("Delta R reco jvf","Delta R; #Delta R; Events",20,0,2);
    TH1F *hist_DR_track_truth = new TH1F("Delta R tracks","Delta R; #Delta R; Events",20,0,2);
    TH1F *hist_DR_truth20_jet = new TH1F("Ratio Jet","Ratio Jet;jetpT/TruthpT;Events",20,0,3.5);
    TH1F *hist_DR_truth20_track = new TH1F("Ratio Track","Ratio Track;trackpT/TruthpT;Events",20,0,3.5);
    TH1F *hist_DR_truth100_jet = new TH1F("Ratio Jet","Ratio Jet;jetpT/TruthpT;Events",20,0,3.5);
    TH1F *hist_DR_truth100_track = new TH1F("Ratio Track","Ratio Track;trackpT/TruthpT;Events",20,0,3.5);
    TH1F *hist_DR_truth500_jet = new TH1F("Ratio Jet","Ratio Jet;jetpT/TruthpT;Events",20,0,3.5);
    TH1F *hist_DR_truth500_track = new TH1F("Ratio Track","Ratio Track;trackpT/TruthpT;Events",20,0,3.5);

    TCanvas *canvas = new TCanvas("Canvas","",1000,1000);
    TCanvas *canvas_1 = new TCanvas("Canvas1","",1000,1000);
    TCanvas *canvas_2 = new TCanvas("Canvas2","",1000,1000);

    int nentries, nbytes, i;
    nentries = (Int_t)tree->GetEntries();

    for (i = 0; i < nentries; i++)
    {
        nbytes = tree->GetEntry(i);

        if(truth_R4_pt->size()!=0 && truth_R4_pt->at(0)>20000.){
            TLorentzVector truthJet;
            truthJet.SetPtEtaPhiM(truth_R4_pt->at(0),truth_R4_eta->at(0),truth_R4_phi->at(0),truth_R4_m->at(0));        
            if(reco_R4_pt->size()!=0 && fabs(reco_R4_jvf->at(0))>0.5){
                TLorentzVector recoJetJvf;
                recoJetJvf.SetPtEtaPhiM(reco_R4_pt->at(0),reco_R4_eta->at(0),reco_R4_phi->at(0),reco_R4_m->at(0));        
         
                //Plot the Delta R
                hist_DR_recojvf_truth->Fill(truthJet.DeltaR(recoJetJvf),evtw);
            }
            if(reco_R4_pt->size()>0){
                TLorentzVector recoJet;
                recoJet.SetPtEtaPhiM(reco_R4_pt->at(0),reco_R4_eta->at(0),reco_R4_phi->at(0),reco_R4_m->at(0));
                hist_DR_reco_truth->Fill(truthJet.DeltaR(recoJet),evtw);
            }
            if(reco_tracks_R4_pt->size()>0){
                TLorentzVector trackJet;
                trackJet.SetPtEtaPhiM(reco_tracks_R4_pt->at(0),reco_tracks_R4_eta->at(0),reco_tracks_R4_phi->at(0),reco_tracks_R4_m->at(0));
                hist_DR_track_truth->Fill(truthJet.DeltaR(trackJet),evtw);
            }
        }
        if(truth_R4_pt->size()!=0 && truth_R4_pt->at(0)>20000. && truth_R4_pt->at(0)<=100000.){
            TLorentzVector truthJet20;
            truthJet20.SetPtEtaPhiM(truth_R4_pt->at(0),truth_R4_eta->at(0),truth_R4_phi->at(0),truth_R4_m->at(0)); 
            if(reco_R4_pt->size()>0){
                TLorentzVector recoJet2;                
                recoJet2.SetPtEtaPhiM(reco_R4_pt->at(0),reco_R4_eta->at(0),reco_R4_phi->at(0),reco_R4_m->at(0));                
                if(truthJet20.DeltaR(recoJet2)<0.3){
                    hist_DR_truth20_jet->Fill(reco_R4_pt->at(0)/truth_R4_pt->at(0),evtw);
                }
            }
            if(reco_tracks_R4_pt->size()>0){
                TLorentzVector trackJet2;
                trackJet2.SetPtEtaPhiM(reco_tracks_R4_pt->at(0),reco_tracks_R4_eta->at(0),reco_tracks_R4_phi->at(0),reco_tracks_R4_m->at(0));
                if(truthJet20.DeltaR(trackJet2)<0.3){
                    hist_DR_truth20_track->Fill(reco_tracks_R4_pt->at(0)/truth_R4_pt->at(0),evtw);
                }
            }
        }else if(truth_R4_pt->size()!=0 && truth_R4_pt->at(0)>100000. && truth_R4_pt->at(0)<=500000.){
            TLorentzVector truthJet100;
            truthJet100.SetPtEtaPhiM(truth_R4_pt->at(0),truth_R4_eta->at(0),truth_R4_phi->at(0),truth_R4_m->at(0)); 
            if(reco_R4_pt->size()>0){
                TLorentzVector recoJet3;                
                recoJet3.SetPtEtaPhiM(reco_R4_pt->at(0),reco_R4_eta->at(0),reco_R4_phi->at(0),reco_R4_m->at(0));                
                if(truthJet100.DeltaR(recoJet3)<0.3){
                    hist_DR_truth100_jet->Fill(reco_R4_pt->at(0)/truth_R4_pt->at(0),evtw);
                }
            }
            if(reco_tracks_R4_pt->size()>0){
                TLorentzVector trackJet3;
                trackJet3.SetPtEtaPhiM(reco_tracks_R4_pt->at(0),reco_tracks_R4_eta->at(0),reco_tracks_R4_phi->at(0),reco_tracks_R4_m->at(0));
                if(truthJet100.DeltaR(trackJet3)<0.3){
                    hist_DR_truth100_track->Fill(reco_tracks_R4_pt->at(0)/truth_R4_pt->at(0),evtw);
                }
            }
        }else if(truth_R4_pt->size()!=0 && truth_R4_pt->at(0)>500000.){
            TLorentzVector truthJet500;
            truthJet500.SetPtEtaPhiM(truth_R4_pt->at(0),truth_R4_eta->at(0),truth_R4_phi->at(0),truth_R4_m->at(0)); 
            if(reco_R4_pt->size()>0){
                TLorentzVector recoJet4;                
                recoJet4.SetPtEtaPhiM(reco_R4_pt->at(0),reco_R4_eta->at(0),reco_R4_phi->at(0),reco_R4_m->at(0));                
                if(truthJet500.DeltaR(recoJet4)<0.3){
                    hist_DR_truth500_jet->Fill(reco_R4_pt->at(0)/truth_R4_pt->at(0),evtw);
                }
            }
            if(reco_tracks_R4_pt->size()>0){
                TLorentzVector trackJet4;
                trackJet4.SetPtEtaPhiM(reco_tracks_R4_pt->at(0),reco_tracks_R4_eta->at(0),reco_tracks_R4_phi->at(0),reco_tracks_R4_m->at(0));
                if(truthJet500.DeltaR(trackJet4)<0.3){
                    hist_DR_truth500_track->Fill(reco_tracks_R4_pt->at(0)/truth_R4_pt->at(0),evtw);
                }
            }
        }
    }
    canvas->cd();
    hist_DR_recojvf_truth->SetMarkerColor(1);
    hist_DR_recojvf_truth->SetLineColor(1);
    hist_DR_recojvf_truth->Scale(1/hist_DR_recojvf_truth->Integral());
    hist_DR_recojvf_truth->DrawNormalized("");
    hist_DR_reco_truth->Scale(1/hist_DR_reco_truth->Integral());
    hist_DR_reco_truth->SetMarkerColor(2);
    hist_DR_reco_truth->SetLineColor(2);
    hist_DR_reco_truth->DrawNormalized("same");
    hist_DR_track_truth->Scale(1/hist_DR_track_truth->Integral());
    hist_DR_track_truth->SetMarkerColor(4);
    hist_DR_track_truth->SetLineColor(4);
    hist_DR_track_truth->DrawNormalized("same");
    auto legend = new TLegend(0.3,0.7,0.48,0.85);
    legend->AddEntry(hist_DR_recojvf_truth,"Delta R Reco with JVF");
    legend->AddEntry(hist_DR_reco_truth,"Delta R Reco without JVF");
    legend->AddEntry(hist_DR_track_truth,"Delta R Tracks");
    legend->Draw();
    canvas->SetLogy();
    canvas->Print("deltaRJets.ps");

    canvas_1->cd();
    hist_DR_truth20_jet->SetMarkerColor(1);
    hist_DR_truth20_jet->SetLineColor(1);
    hist_DR_truth20_jet->Scale(1/hist_DR_truth20_jet->Integral());
    hist_DR_truth20_jet->DrawNormalized("");
    hist_DR_truth100_jet->SetMarkerColor(2);
    hist_DR_truth100_jet->SetLineColor(2);
    hist_DR_truth100_jet->Scale(1/hist_DR_truth100_jet->Integral());
    hist_DR_truth100_jet->DrawNormalized("same");
    hist_DR_truth500_jet->SetMarkerColor(4);
    hist_DR_truth500_jet->SetLineColor(4);
    hist_DR_truth500_jet->Scale(1/hist_DR_truth500_jet->Integral());
    hist_DR_truth500_jet->DrawNormalized("same");
    auto legend1 = new TLegend(0.55,0.73,0.7,0.85);
    legend1->AddEntry(hist_DR_truth20_jet,"Leading truth pT>20");
    legend1->AddEntry(hist_DR_truth100_jet,"Leading truth pT>100");
    legend1->AddEntry(hist_DR_truth500_jet,"Leading truth pT>500");
    legend1->Draw();
    canvas_1->SetLogy();
    canvas_1->Print("jetRatio.ps");

    canvas_2->cd();
    hist_DR_truth20_track->SetMarkerColor(1);
    hist_DR_truth20_track->SetLineColor(1);
    hist_DR_truth20_track->Scale(1/hist_DR_truth20_track->Integral());
    hist_DR_truth20_track->DrawNormalized("");
    hist_DR_truth100_track->SetMarkerColor(2);
    hist_DR_truth100_track->SetLineColor(2);
    hist_DR_truth100_track->Scale(1/hist_DR_truth100_track->Integral());
    hist_DR_truth100_track->DrawNormalized("same");
    hist_DR_truth500_track->SetMarkerColor(4);
    hist_DR_truth500_track->SetLineColor(4);
    hist_DR_truth500_track->Scale(1/hist_DR_truth500_track->Integral());
    hist_DR_truth500_track->DrawNormalized("same");
    auto legend2 = new TLegend(0.52,0.7,0.7,0.85);
    legend2->AddEntry(hist_DR_truth20_track,"Leading truth pT>20");
    legend2->AddEntry(hist_DR_truth100_track,"Leading truth pT>100");
    legend2->AddEntry(hist_DR_truth500_track,"Leading truth pT>500");
    legend2->Draw();
    canvas_2->SetLogy();
    canvas_2->Print("jetTrackRatio.ps");
}
