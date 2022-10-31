#include <iostream>
#include <string>
#include <stdio.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"
#include "TH1F.h"

void allVariablesTracksAndClusters(){
    TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
    TTree *tree = (TTree*) file->Get("JetRecoTree");
    vector<float> *tracks_pt=0;
    vector<int> *tracks_vtx=0;
    vector<float> *tracks_eta=0;
    vector<float> *tracks_phi=0;
    vector<float> *tracks_mass=0;
    vector<float> *clusters_pt=0;
    vector<float> *clusters_eta=0;
    vector<float> *clusters_phi=0;
    vector<float> *clusters_mass=0;
    double pi = M_PI;
    
    tree->SetBranchAddress("Tracks_pt", &tracks_pt);
    tree->SetBranchAddress("Tracks_vtx", &tracks_vtx);
    tree->SetBranchAddress("Tracks_eta", &tracks_eta);
    tree->SetBranchAddress("Tracks_phi", &tracks_phi);
    tree->SetBranchAddress("Tracks_m", &tracks_mass);
    tree->SetBranchAddress("Clusters_pt", &clusters_pt);
    tree->SetBranchAddress("Clusters_eta", &clusters_eta);
    tree->SetBranchAddress("Clusters_phi", &clusters_phi);
    tree->SetBranchAddress("Clusters_m", &clusters_mass);
    TCanvas *canvas = new TCanvas("Canvas","allVariables",800,600);

    TH1F *hist_lead_track_pt = new TH1F("Lead_Track_pT","Track pT; pT (MeV) ; Events ",400,450,4000);
    TH1F *hist_lead_track_phi = new TH1F("Lead_Track_phi","Example plot: Track phi; phi (°) ; Events ",100,-pi,pi);
    TH1F *hist_lead_track_eta = new TH1F("Lead_Track_eta","Example plot: Track eta; eta(°) ; Events ",100,-5.,5.);
    TH1F *hist_lead_track_mass = new TH1F("Lead_Track_mass","Example plot: Track mass; mass (MeV/c^{2}) ; Events ",50,130,150);
    TH1F *hist_lead_track_vtx = new TH1F("Lead_Track_vtx","Example plot: Track vtx",100,1,100);
    TH1F *hist_clusters_mass = new TH1F("Clusters mass","Example plot: clusters mass; mass (MeV/c^{2}); Events ",1000,-1.,1.);    
    TH1F *hist_clusters_eta = new TH1F("Clusters eta","Example plot: clusters eta; pT (MeV) ; Events ",300,-5.,5.);
    TH1F *hist_clusters_phi = new TH1F("Clusters phi","Example plot: clusters phi; eta (°) ; Events ",100,-pi,pi);
    TH1F *hist_clusters_pt = new TH1F("Clusters pT","Example plot: clusters pT; pT (MeV) ; Events ",1000,0.,10000.);
    
    int nentries, nbytes, i;
    nentries = (Int_t)tree->GetEntries();

    for (i = 0; i < nentries; i++)
    {
        nbytes = tree->GetEntry(i);
        for(int tr=0; tr<tracks_pt->size(); tr++)
        {
            hist_lead_track_pt->Fill(tracks_pt->at(tr));
        }
        for(int tr=0; tr<tracks_phi->size(); tr++){
            hist_lead_track_phi->Fill(tracks_phi->at(tr));
        }
        for(int tr=0; tr<tracks_eta->size(); tr++){
            hist_lead_track_eta->Fill(tracks_eta->at(tr));
        }
        for(int tr=0; tr<tracks_mass->size(); tr++){
            hist_lead_track_mass->Fill(tracks_mass->at(tr));
        }
        for(int tr=0; tr<clusters_pt->size(); tr++){
            hist_clusters_pt->Fill(clusters_pt->at(tr));
        }
        for(int tr=0; tr<clusters_phi->size(); tr++){
            hist_clusters_phi->Fill(clusters_phi->at(tr));
        }
        for(int tr=0; tr<clusters_eta->size(); tr++){
            hist_clusters_eta->Fill(clusters_eta->at(tr));
        }
        for(int tr=0; tr<clusters_eta->size(); tr++){
           hist_clusters_mass->Fill(clusters_mass->at(tr)); 
        }
        for(int tr=0; tr<tracks_vtx->size(); tr++)
        {
            hist_lead_track_vtx->Fill(tracks_vtx->at(tr));
        }
        
    }
    /*
    hist_lead_track_pt->SetFillColor(kBlue);
    hist_lead_track_pt->Draw();
    canvas->Draw();
    
    hist_lead_track_eta->SetFillColor(1);
    hist_lead_track_eta->Draw();
    canvas->Draw();

    hist_lead_track_phi->SetFillColor(2);
    hist_lead_track_phi->Draw();
    canvas->Draw();

    hist_lead_track_mass->SetFillColor(3);
    hist_lead_track_mass->Draw();
    canvas->Draw();

    hist_clusters_pt->SetFillColor(4);
    hist_clusters_pt->Draw();
    canvas->Draw();

    hist_clusters_eta->SetFillColor(5);
    hist_clusters_eta->Draw();
    canvas->Draw();

    hist_clusters_phi->SetFillColor(6);
    hist_clusters_phi->Draw();
    canvas->Draw();

    hist_clusters_mass->SetFillColor(7);
    hist_clusters_mass->Draw();
    canvas->Draw();*/
   
    hist_lead_track_vtx->SetFillColor(kBlue);
    hist_lead_track_vtx->Draw();
    canvas->Draw(); 

}
