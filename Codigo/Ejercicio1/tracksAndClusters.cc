#include <iostream>
#include <string>
#include <stdio.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"
#include "TH1F.h"

void tracksAndClusters(){
    TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
    TTree *tree = (TTree*) file->Get("JetRecoTree");
    Float_t muaverage = -1;
    UInt_t npv = -1;
    vector<float> *tracks_pt=0;
    vector<float> *clusters_pt=0;

    tree->SetBranchAddress("Tracks_pt", &tracks_pt);
    tree->SetBranchAddress("Clusters_pt", &clusters_pt);
    tree->SetBranchAddress("mu_average",&muaverage);
    tree->SetBranchAddress("NPV",&npv);
    TCanvas *canvas = new TCanvas("Canvas","tracksAndClusters",800,600);
    TH2F *hist_npv_vs_tracks = new TH2F("Tracks per PV","NPV Vs. Tracks; npv ;NTracks",50,1,50,270,50,2700);
    TH2F *hist_npv_vs_clusters = new TH2F("Clusters per PV","NPV Vs. Clusters; npv ;NClusters",50,1,50,160,1,1600);
    TH2F *hist_mu_vs_tracks = new TH2F("Tracks per Avg. Int.","Mu Vs. Tracks; mu average ;NTracks",85,1,85,270,50,2700);
    TH2F *hist_mu_vs_clusters = new TH2F("Clusters per Avg. Inter","Mu Vs. Clusters; mu average ;NClusters",60,1,85,160,1,1600);
    TH1F *hist_nTrks = new TH1F("NTracks","Number of tracks; NTracks ; Events ",270,0,2700);
    TH1F *hist_nCltrs = new TH1F("NClusters","Number of clusters; NClusters ; Events ",160,0,1600);
    
    int nentries, nbytes, i;
    nentries = (Int_t)tree->GetEntries();

    for (i = 0; i < nentries; i++)
    {
        nbytes = tree->GetEntry(i);
        hist_nTrks->Fill(tracks_pt->size());
        hist_nCltrs->Fill(clusters_pt->size());
        hist_npv_vs_tracks->Fill(npv,tracks_pt->size());
        hist_npv_vs_clusters->Fill(npv,clusters_pt->size());
        hist_mu_vs_tracks->Fill(muaverage,tracks_pt->size());
        hist_mu_vs_clusters->Fill(muaverage,clusters_pt->size());
    }
/*
    hist_nTrks->SetFillColor(kBlue);
    hist_nTrks->Draw();
    canvas->Draw();

    hist_nCltrs->SetFillColor(kYellow);
    hist_nCltrs->Draw();
    canvas->Draw();

    hist_npv_vs_tracks->SetFillColor(kBlue);
    hist_npv_vs_tracks->Draw("colors");
    canvas->Draw();

    hist_npv_vs_clusters->SetFillColor(kRed);
    hist_npv_vs_clusters->Draw("colors");
    canvas->Draw();

    hist_mu_vs_clusters->SetFillColor(kRed);
    hist_mu_vs_clusters->Draw("colors");
    canvas->Draw();*/

    hist_mu_vs_tracks->SetFillColor(kRed);
    hist_mu_vs_tracks->Draw("colors");
    canvas->Draw();
    
}
