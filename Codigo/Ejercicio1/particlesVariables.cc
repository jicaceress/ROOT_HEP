#include <iostream>
#include <string>
#include <stdio.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"


void particlesVariables(){
    TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
    TTree *tree = (TTree*) file->Get("JetRecoTree");

    vector<float> *particles_pt = 0;    
    vector<float> *particles_eta = 0;
    vector<float> *particles_phi = 0;
    vector<float> *particles_mass = 0;
    vector<int> *particles_pdgid = 0;
    double pi = M_PI;

    tree->SetBranchAddress("Particles_pt", &particles_pt);
    tree->SetBranchAddress("Particles_eta", &particles_eta);
    tree->SetBranchAddress("Particles_phi", &particles_phi);
    tree->SetBranchAddress("Particles_m", &particles_mass);
    tree->SetBranchAddress("Particles_pdgID", &particles_pdgid);

    TCanvas *canvas = new TCanvas("Canvas","allVariables",800,600);
    
    TH1F *hist_particles_pt = new TH1F("Particles_pT","Particles pT; pT (MeV) ; Events ",600,0,6000);
    TH1F *hist_particles_phi = new TH1F("Paricles_phi","Particles phi; phi (°) ; Events ",100,-pi,pi);
    TH1F *hist_particles_eta = new TH1F("Particles_eta","Particles eta; eta(°) ; Events ",100,-5.,5.);
    TH1F *hist_particles_mass = new TH1F("Particles_mass","Particles mass; mass (MeV/c^{2}) ; Events ",50,130,150);
    TH1F *hist_particles_pgdid = new TH1F("Particles_pgdID","Particles pdg ID",100,1,100);

    int nentries, nbytes, i;
    nentries = (Int_t)tree->GetEntries();

    for (i = 0; i < nentries; i++)
    {
        nbytes = tree->GetEntry(i);
        for(int tr=0; tr<particles_pt->size(); tr++)
        {
            hist_particles_pt->Fill(particles_pt->at(tr));
        }
        for(int tr=0; tr<particles_phi->size(); tr++){
            hist_particles_phi->Fill(particles_phi->at(tr));
        }
        for(int tr=0; tr<particles_eta->size(); tr++){
            hist_particles_eta->Fill(particles_eta->at(tr));
        }
        for(int tr=0; tr<particles_mass->size(); tr++){
            hist_particles_mass->Fill(particles_mass->at(tr));
        }   
        for(int tr=0; tr<particles_pgdid->size(); tr++){
            hist_particles_pdgid->Fill(particles_pgdid->at(tr));
        }      
    }
/*
    hist_particles_pt->SetFillColor(kBlue);
    hist_particles_pt->Draw();
    canvas->Draw();
    hist_particles_phi->SetFillColor(kBlue);
    hist_particles_phi->Draw();
    canvas->Draw();
    hist_particles_eta->SetFillColor(kBlue);
    hist_particles_eta->Draw();
    canvas->Draw();
    hist_particles_mass->SetFillColor(kBlue);
    hist_particles_mass->Draw();
    canvas->Draw();*/
    hist_particles_pdgid->SetFillColor(kBlue);
    hist_particles_pdgid->Draw();
    canvas->Draw();
}
