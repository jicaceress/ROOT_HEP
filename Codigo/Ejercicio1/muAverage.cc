#include <iostream>
#include <string>
#include <stdio.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"

void muAverage(){
    TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/Tracks_Clusters.root");
    TTree *tree = (TTree*) file->Get("JetRecoTree");
    Float_t muaverage = -1;
    tree->SetBranchAddress("mu_average", &muaverage);
    TCanvas *canvas = new TCanvas("Canvas","muAverage",800,600);
    TH1F *hist_mu_average = new TH1F("mu Average","Number of average interactions; mu average ; Events ",30,1,85);
    int nentries, nbytes, i;
    nentries = (Int_t)tree->GetEntries();

    for (i = 0; i < nentries; i++)
    {
        nbytes = tree->GetEntry(i);
        hist_mu_average->Fill(muaverage);
    }
   
    std::cout << "Done!" << std::endl;
    hist_mu_average->SetFillColor(kRed);
    hist_mu_average->Draw();
    canvas->Draw();
    
}
