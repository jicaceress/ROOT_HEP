#include <iostream>
#include <string>
#include <stdio.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"

void npvVsMu(){
TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
TTree *tree = (TTree*) file->Get("JetRecoTree");
Float_t muaverage = -1;
UInt_t npv = -1;
tree->SetBranchAddress("mu_average",&muaverage);
tree->SetBranchAddress("NPV",&npv);

TCanvas *canvas = new TCanvas("Canvas","NPV Vs. muAverage",800,600);
TH2F *hist_npv_vs_mu_average = new TH2F("Avg. Int. per PV ","NPV Vs. Mu Average; npv ;mu average",60,1,60,45,1,90);
int nentries, nbytes, i;
    nentries = (Int_t)tree->GetEntries();

    for (i = 0; i < nentries; i++)
    {
        nbytes = tree->GetEntry(i);
        hist_npv_vs_mu_average->Fill(npv,muaverage);
    }
    hist_npv_vs_mu_average->SetFillColor(25);
    hist_npv_vs_mu_average->Draw("colors");
    canvas->Draw();
}
