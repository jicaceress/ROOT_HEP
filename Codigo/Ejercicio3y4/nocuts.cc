#include <iostream>
#include <string>
#include <stdio.h>
#include"TTree.h"
#include"TFile.h"
#include"TCanvas.h"
#include"TH1F.h"
#include"TH2F.h"
#include"TLorentzVector.h"

void nocuts(){


/*opening the root file*/

TFile *file = TFile::Open("/home/fernando/Desktop/HEP/Part2/Data_8TeV.root");
TTree *tree = (TTree*) file->Get("mini");
TFile *file_MC = TFile::Open("/home/fernando/Desktop/HEP/Part2/ttbar_8TeV.root");  
TTree *tree_MC = (TTree*) file->Get("mini");
/*setting the variables*/
Bool_t e_trig;
Bool_t mu_trig;
Bool_t good_vtx;
UInt_t lep_n;
UInt_t jet_n;
UInt_t evtw;
Float_t MET;
Float_t mTc;
Float_t MET_phi;

Float_t jet_MV1[10];
Float_t lep_pt[10];
Float_t lep_eta[10];
Float_t lep_phi[10];
Float_t lep_E[10];
Int_t lep_type[10];
Float_t lep_ptcone30[10];
Float_t lep_etcone20[10];

Float_t jet_pt[10];
Float_t jet_et[10];

Float_t jet_phi[10];
Float_t jet_eta[10];
Float_t jet_jvf[10];


//data variables for the MC

Bool_t e_trig_MC;
Bool_t mu_trig_MC;
Bool_t good_vtx_MC;
UInt_t lep_n_MC;
UInt_t jet_n_MC;
UInt_t evtw_MC;
Float_t MET_MC;
Float_t mTc_MC;
Float_t MET_phi_MC;

Float_t jet_MV1_MC[10];
Float_t lep_pt_MC[10];
Float_t lep_eta_MC[10];
Float_t lep_phi_MC[10];
Float_t lep_E_MC[10];
Int_t lep_type_MC[10];
Float_t lep_ptcone30_MC[10];
Float_t lep_etcone20_MC[10];

Float_t jet_pt_MC[10];
Float_t jet_et_MC[10];

Float_t jet_phi_MC[10];
Float_t jet_eta_MC[10];
Float_t jet_jvf_MC[10];

Float_t scaleFactor_PILEUP;
Float_t scaleFactor_ELE;
Float_t scaleFactor_MUON;
Float_t scaleFactor_BTAG;
Float_t scaleFactor_TRIGGER;
Float_t scaleFactor_JVFSF;
Float_t scaleFactor_ZVERTEX;
Float_t scaleFactor_MCFactor;


/*inicializing the variables*/

tree->SetBranchAddress("trigE", &e_trig);
tree->SetBranchAddress("trigM", &mu_trig);
tree->SetBranchAddress("hasGoodVertex", &good_vtx);
tree->SetBranchAddress("lep_n", &lep_n);
tree->SetBranchAddress("jet_n", &jet_n);
tree->SetBranchAddress("jet_MV1", &jet_MV1);
tree->SetBranchAddress("met_et", &MET);
tree->SetBranchAddress("met_phi", &MET_phi);
tree->SetBranchAddress("jet_E", &jet_et);
tree->SetBranchAddress("jet_phi", &jet_phi);
tree->SetBranchAddress("jet_eta", &jet_eta);
tree->SetBranchAddress("jet_jvf", &jet_jvf);
tree->SetBranchAddress("lep_pt", &lep_pt);
tree->SetBranchAddress("lep_eta", &lep_eta);
tree->SetBranchAddress("lep_phi", &lep_phi);
tree->SetBranchAddress("lep_E", &lep_E);
tree->SetBranchAddress("lep_type", &lep_type);
tree->SetBranchAddress("lep_ptcone30", &lep_ptcone30);
tree->SetBranchAddress("lep_etcone20", &lep_etcone20);
tree->SetBranchAddress("jet_pt", &jet_pt);


/*setting the branch addresses for the MC data set*/

tree_MC->SetBranchAddress("scaleFactor_PILEUP", &scaleFactor_PILEUP);
tree_MC->SetBranchAddress("scaleFactor_ELE", &scaleFactor_ELE);
tree_MC->SetBranchAddress("scaleFactor_MUON", &scaleFactor_MUON);
tree_MC->SetBranchAddress("scaleFactor_BTAG", &scaleFactor_BTAG);
tree_MC->SetBranchAddress("scaleFactor_TRIGGER", &scaleFactor_TRIGGER);
tree_MC->SetBranchAddress("scaleFactor_JVFSF", &scaleFactor_JVFSF);
tree_MC->SetBranchAddress("scaleFactor_ZVERTEX", &scaleFactor_ZVERTEX);
tree_MC->SetBranchAddress("scaleFactor_ZVERTEX", &scaleFactor_MCFactor);

tree_MC->SetBranchAddress("trigE", &e_trig_MC);
tree_MC->SetBranchAddress("trigM", &mu_trig_MC);
tree_MC->SetBranchAddress("hasGoodVertex", &good_vtx_MC);
tree_MC->SetBranchAddress("lep_n", &lep_n_MC);
tree_MC->SetBranchAddress("jet_n", &jet_n_MC);
tree_MC->SetBranchAddress("jet_MV1", &jet_MV1_MC);
tree_MC->SetBranchAddress("met_et", &MET_MC);
tree_MC->SetBranchAddress("met_phi", &MET_phi_MC);
tree_MC->SetBranchAddress("jet_E", &jet_et_MC);
tree_MC->SetBranchAddress("jet_phi", &jet_phi_MC);
tree_MC->SetBranchAddress("jet_eta", &jet_eta_MC);
tree_MC->SetBranchAddress("jet_jvf", &jet_jvf_MC);
tree_MC->SetBranchAddress("lep_pt", &lep_pt_MC);
tree_MC->SetBranchAddress("lep_eta", &lep_eta_MC);
tree_MC->SetBranchAddress("lep_phi", &lep_phi_MC);
tree_MC->SetBranchAddress("lep_E", &lep_E_MC);
tree_MC->SetBranchAddress("lep_type", &lep_type_MC);
tree_MC->SetBranchAddress("lep_ptcone30", &lep_ptcone30_MC);
tree_MC->SetBranchAddress("lep_etcone20", &lep_etcone20_MC);
tree_MC->SetBranchAddress("jet_pt", &jet_pt_MC);
		


/*Creating the Canvas to Draw the histograms*/

TCanvas *canvas_1 = new TCanvas("Canvas1","",1000,1000);
TCanvas *canvas_2 = new TCanvas("Canvas2","",1000,1000);
TCanvas *canvas_3 = new TCanvas("Canvas3","",1000,1000);
TCanvas *canvas_4 = new TCanvas("Canvas4","",1000,1000);
TCanvas *canvas_5 = new TCanvas("Canvas5","",1000,1000);
TCanvas *canvas_6 = new TCanvas("Canvas6","",1000,1000);
TCanvas *canvas_7 = new TCanvas("Canvas7","",1000,1000);
TCanvas *canvas_8 = new TCanvas("Canvas8","",1000,1000);
TCanvas *canvas_9 = new TCanvas("Canvas9","",1000,1000);
TCanvas *canvas_10 = new TCanvas("Canvas10","",1000,1000);
TCanvas *canvas_11 = new TCanvas("Canvas11","",1000,1000);
TCanvas *canvas_12 = new TCanvas("Canvas12","",1000,1000);

TCanvas *canvas_13 = new TCanvas("Canvas13","",1000,1000);
TCanvas *canvas_14 = new TCanvas("Canvas14","",1000,1000);
TCanvas *canvas_15 = new TCanvas("Canvas15","",1000,1000);
TCanvas *canvas_16 = new TCanvas("Canvas16","",1000,1000);
TCanvas *canvas_17 = new TCanvas("Canvas17","",1000,1000);
TCanvas *canvas_18 = new TCanvas("Canvas18","",1000,1000);
TCanvas *canvas_19 = new TCanvas("Canvas19","",1000,1000);
TCanvas *canvas_20 = new TCanvas("Canvas20","",1000,1000);
TCanvas *canvas_21 = new TCanvas("Canvas21","",1000,1000);
TCanvas *canvas_22 = new TCanvas("Canvas22","",1000,1000);
TCanvas *canvas_23 = new TCanvas("Canvas23","",1000,1000);
TCanvas *canvas_24 = new TCanvas("Canvas24","",1000,1000);

/*creating the histograms to draw*/

TH1F *cutflow = new TH1F("Cutflow","Cutflow; Cut; Events",10,0,10);
TH1F *lepton_t_isolation = new TH1F("Lepton Track Isolation","Leptons Track Isolation; ; Events",500,0., 0.25);
TH1F *lepton_isolation = new TH1F("Lepton Calorimeter Isolation","Leptons Calorimeter Isolation; ; Events",500,0.,0.25);
TH1F *lepton_eta = new TH1F("Lepton eta","Leptons eta; ; Events",500,-5., 5.);
TH1F *lepton_pt = new TH1F("Lepton PT","Leptons PT; ; Events",500,-20000., 150000);
TH1F *hist_nbjets = new TH1F("Number of bjets","n-bjets; Jet multiplicity; Events",10,0,10);
TH1F *hist_jets_pt = new TH1F("Jets Pt","jets Pt; ; Events",500,25000.,150000.);
TH1F *hist_jets_eta = new TH1F("Jets eta","jets eta; ; Events",500,-5.,5.);
TH1F *hist_jets_JVF = new TH1F("Jets JVF","jets JVF; ; Events",500,-1.,1.);
TH1F *hist_jets_MV1 = new TH1F("Jets Mv1","jets MV1; ; Events",500,0.,5.);
TH1F *hist_MET = new TH1F("MET","MET; ; Events",500,25000.,150000.);
TH1F *hist_mTc = new TH1F("Transverse Mass","Transverse Mass; ; Events",500,0.,150000.);


TH1F *cutflow_MC = new TH1F("Cutflow_MC","Cutflow_MC; Cut; Events",10,0,10);
TH1F *lepton_t_isolation_MC = new TH1F("Lepton Track Isolation_MC","Leptons Track Isolation; ; Events",500,0., 0.25);
TH1F *lepton_isolation_MC = new TH1F("Lepton Calorimeter Isolation_MC","Leptons Calorimeter Isolation; ; Events",500,0.,0.25);
TH1F *lepton_eta_MC = new TH1F("Lepton eta_MC","Leptons eta; ; Events",500,-5., 5.);
TH1F *lepton_pt_MC = new TH1F("Lepton PT_MC","Leptons PT; ; Events",500,-20000., 150000);
TH1F *hist_nbjets_MC = new TH1F("Number of bjets_MC","n-bjets; Jet multiplicity; Events",10,0,10);
TH1F *hist_jets_pt_MC = new TH1F("Jets Pt_MC","jets Pt; ; Events",500,25000.,150000.);
TH1F *hist_jets_eta_MC = new TH1F("Jets eta_MC","jets eta; ; Events",500,-5.,5.);
TH1F *hist_jets_JVF_MC = new TH1F("Jets JVF_MC","jets JVF; ; Events",500,-1.,1.);
TH1F *hist_jets_MV1_MC = new TH1F("Jets Mv1_MC","jets MV1; ; Events",500,0.,5.);
TH1F *hist_MET_MC = new TH1F("MET_MC","MET; ; Events",500,25000.,150000.);
TH1F *hist_mTc_MC = new TH1F("Transverse Mass_MC","Transverse Mass; ; Events",500,0.,150000.);

/*loop and fill histograms*/

int nentries, nbytes, i;
nentries = (Int_t)tree->GetEntries();

int cut1 = 0;
int cut2 = 0;
int cut3 = 0;
int cut4 = 0;
int cut5 = 0;
int cut6 = 0;
int cut7 = 0;
int cut8 = 0;


for (i = 0; i < nentries; i++)
{
    tree->GetEntry(i);   
     //First cut: Good vertex
    if(!good_vtx) continue;
    cut1++;
    cutflow->Fill(1);

    //Second cut: Trigger
    if(!e_trig && !mu_trig) continue;
    cut2++;
    cutflow->Fill(2);
        
        
// Preselection of good leptons                                                                                
    int n_mu=0;
    int n_el=0;
    int n_lep=0;
    int good_lep = 0;	
    //Loop over leptons
    for(unsigned int j=0; j<lep_n; j++){
        if( lep_pt[j] < 25000.) continue; 
        if( lep_ptcone30[j]/lep_pt[j] > 0.15 ) continue; 
        if( lep_etcone20[j]/lep_pt[j] > 0.15 ) continue;  
        if( lep_type [j]==13 && TMath::Abs(lep_eta[j]) < 2.5 ){
            n_mu++;
            good_lep=j;
            }
	if( lep_type [j]==11 && (TMath::Abs(lep_eta[j]) < 1.37 || (1.52 < TMath::Abs(lep_eta[j]) && TMath::Abs(lep_eta[j]) < 2.47)) ){
	n_el++;
	good_lep=j;
	//std::cout << n_el << "Number of el" << std::endl;
	//}
	}
	
	lepton_pt->Fill(lep_pt[good_lep]);
	lepton_eta->Fill(lep_eta[good_lep]);
	lepton_t_isolation->Fill(lep_ptcone30[good_lep]/lep_pt[good_lep]);
	lepton_isolation->Fill(lep_etcone20[good_lep]/lep_pt[good_lep]);
        /*
        && 
        || 
        To complete: Add electrons and extract the index for the good lepton
        */
    }
    n_lep = n_el + n_mu;
    if(n_lep!=1) continue;
    cutflow->Fill(3); 
    cut3++;
    
     int n_jets=0;
    int n_bjets=0;
    
     //Number of jets distribution
   //hist_njets->Fill(jet_n,evtw);
    
    //at least four jets
    if(jet_n<4) continue; 
    cutflow->Fill(4); 
    cut4++;
     //at least four GOOD jets
     
    for(unsigned int j=0; j<jet_n; j++){
        // To complete: apply jet cuts to find the good jets
        if(jet_pt[j] < 25000.) continue;
        //Eta cut
        if (TMath::Abs(jet_eta[j]) > 2.5 ) continue;
        // JVF Cleaning 
        if(jet_pt[j] < 50000. && TMath::Abs(jet_eta[j]) < 2.4 && jet_jvf[j] < 0.5) continue;
        // cut on 0.7892 MV1 and count the number of b-jets
        //if (jet_MV1 >= 0.7892) continue;
        if (jet_MV1[j] >= 0.7892){
        n_bjets++;
        hist_nbjets->Fill(n_bjets);
        }
        n_jets++;
        hist_jets_pt->Fill(jet_pt[j]);
        hist_jets_eta->Fill(jet_eta[j]); 
        hist_jets_JVF->Fill(jet_jvf[j]);
        hist_jets_MV1->Fill(jet_MV1[j]); 
        }
    
    //Fifth cut: At least 4 good jets
    if(n_jets<4) continue; 
    cutflow->Fill(5); 
    cut5++;
    
    if(n_bjets <= 1) continue;
    cutflow->Fill(6);
    cut6++;
    
    //MET cut
    if(MET<30000.) continue;
    hist_MET->Fill(MET);
    cutflow->Fill(7);
    cut7++;
    // Fede was here and he was kinda gay tho, no homo. 
    //TLorentzVector  Lepton(0.,0.,0.,0.);
    //TLorentzVector  MeT(0.,0.,0.,0.);
    
    if(mTc == sqrt(2*lep_pt[good_lep]*MET*(1-cos(MET_phi-lep_phi[good_lep]))) < 30000.) continue;
    hist_mTc->Fill(mTc);
    cutflow->Fill(8);	
    cut8++;
    
    
    
    
    }

   
std::cout << "Done!" << std::endl;
std::cout << "All events:" << nentries << std::endl;
std::cout << "Cut1:" << cut1 << std::endl;
std::cout << "Cut2:" << cut2 << std::endl;
std::cout << "Cut3:" << cut3 << std::endl;
std::cout << "Cut4:" << cut4 << std::endl;
std::cout << "Cut5:" << cut5 << std::endl;
std::cout << "Cut6:" << cut6 << std::endl;
std::cout << "Cut7:" << cut7 << std::endl;
std::cout << "Cut8:" << cut8 << std::endl;



//resetting the cuts to make them using the tt_bar MC data sets.
    
cut8=cut7=cut6=cut5=cut4=cut3=cut2=cut1=0;

/*
int nentries_MC, nbytes_MC;
nentries_MC = (Int_t)tree_MC->GetEntries();
    
for (i = 0; i < nentries_MC; i++)
{
    tree_MC->GetEntry(i);   
     //First cut: Good vertex
    if(!good_vtx_MC) continue;
    cut1++;
    cutflow_MC->Fill(1);

    //Second cut: Trigger
    if(!e_trig_MC && !mu_trig_MC) continue;
    cut2++;
    cutflow_MC->Fill(2);
        
        
// Preselection of good leptons                                                                                
    int n_mu_MC=0;
    int n_el_MC=0;
    int n_lep_MC=0;
    int good_lep_MC = 0;	
    //Loop over leptons
    for(unsigned int j=0; j<lep_n_MC; j++){
        if( lep_pt_MC[j] < 25000.) continue; 
        if( lep_ptcone30_MC[j]/lep_pt_MC[j] > 0.15 ) continue; 
        if( lep_etcone20_MC[j]/lep_pt_MC[j] > 0.15 ) continue;  
        if( lep_type_MC [j]==13 && TMath::Abs(lep_eta_MC[j]) < 2.5 ){
            n_mu_MC++;
            good_lep_MC=j;
            }
	if( lep_type [j]==11 && (TMath::Abs(lep_eta_MC[j]) < 1.37 || (1.52 < TMath::Abs(lep_eta_MC[j]) && TMath::Abs(lep_eta_MC[j]) < 2.47)) ){
	n_el_MC++;
	good_lep_MC=j;
	//std::cout << n_el << "Number of el" << std::endl;
	//}
	}
	
	lepton_pt_MC->Fill(lep_pt_MC[good_lep_MC]);
	lepton_eta_MC->Fill(lep_eta_MC[good_lep_MC]);
	lepton_t_isolation_MC->Fill(lep_ptcone30_MC[good_lep_MC]/lep_pt_MC[good_lep_MC]);
	lepton_isolation_MC->Fill(lep_etcone20_MC[good_lep_MC]/lep_pt_MC[good_lep_MC]);
        
    }
    n_lep_MC = n_el_MC + n_mu_MC;
    if(n_lep_MC!=1) continue;
    cutflow_MC->Fill(3); 
    cut3++;
    
     int n_jets_MC=0;
    int n_bjets_MC=0;
    
     //Number of jets distribution
   //hist_njets->Fill(jet_n,evtw);
    
    //at least four jets
    if(jet_n_MC<4) continue; 
    cutflow_MC->Fill(4); 
    cut4++;
     //at least four GOOD jets
     
    for(unsigned int j=0; j<jet_n_MC; j++){
        // To complete: apply jet cuts to find the good jets
        if(jet_pt_MC[j] < 25000.) continue;
        //Eta cut
        if (TMath::Abs(jet_eta_MC[j]) > 2.5 ) continue;
        // JVF Cleaning 
        if(jet_pt_MC[j] < 50000. && TMath::Abs(jet_eta_MC[j]) < 2.4 && jet_jvf_MC[j] < 0.5) continue;
        // cut on 0.7892 MV1 and count the number of b-jets
        //if (jet_MV1 >= 0.7892) continue;
        if (jet_MV1_MC[j] >= 0.7892){
        n_bjets_MC++;
        hist_nbjets_MC->Fill(n_bjets_MC);
        }
        n_jets_MC++;
        hist_jets_pt_MC->Fill(jet_pt_MC[j]);
        hist_jets_eta_MC->Fill(jet_eta_MC[j]); 
        hist_jets_JVF_MC->Fill(jet_jvf_MC[j]);
        hist_jets_MV1_MC->Fill(jet_MV1_MC[j]); 
        }
    
    //Fifth cut: At least 4 good jets
    if(n_jets_MC<4) continue; 
    cutflow_MC->Fill(5); 
    cut5++;
    
    if(n_bjets_MC <= 1) continue;
    cutflow_MC->Fill(6);
    cut6++;
    
    //MET cut
    if(MET_MC<30000.) continue;
    hist_MET_MC->Fill(MET_MC);
    cutflow_MC->Fill(7);
    cut7++;
    // Fede was here and he was kinda gay tho, no homo. 
    //TLorentzVector  Lepton(0.,0.,0.,0.);
    //TLorentzVector  MeT(0.,0.,0.,0.);
    
    if(mTc_MC == sqrt(2*lep_pt_MC[good_lep_MC]*MET_MC*(1-cos(MET_phi_MC-lep_phi_MC[good_lep_MC]))) < 30000.) continue;
    hist_mTc_MC->Fill(mTc_MC);
    cutflow_MC->Fill(8);	
    cut8++;
    
    
    
    
    }
    */
std::cout << "Done!" << std::endl;
std::cout << "All events:" << nentries << std::endl;
std::cout << "Cut1:" << cut1 << std::endl;
std::cout << "Cut2:" << cut2 << std::endl;
std::cout << "Cut3:" << cut3 << std::endl;
std::cout << "Cut4:" << cut4 << std::endl;
std::cout << "Cut5:" << cut5 << std::endl;
std::cout << "Cut6:" << cut6 << std::endl;
std::cout << "Cut7:" << cut7 << std::endl;
std::cout << "Cut8:" << cut8 << std::endl;
    
   // we recycle the cuts for the MC data
    

canvas_1->cd();
 cutflow->Draw("colz");
 canvas_1->Print("Cutflow.ps");
 
 canvas_2->cd();
 lepton_t_isolation->Draw("colz");
 canvas_2->Print("lep_Track_Isolation.ps");

 canvas_3->cd();
 lepton_isolation->Draw("colz");
 canvas_3->Print("lepton_Isolation.ps");

 canvas_4->cd();
 lepton_eta->Draw("colz");
 canvas_4->Print("lepton_eta.ps");

canvas_5->cd();
 lepton_pt->Draw("colz");
 canvas_5->Print("lepton_pt.ps");


canvas_6->cd();
 hist_jets_pt->Draw("colz");
 canvas_5->Print("hist_jets_pt.ps");

canvas_7->cd();
 hist_jets_eta->Draw("colz");
 canvas_7->Print("hist_jets_eta.ps");

canvas_8->cd();
 hist_jets_JVF->Draw("colz");
 canvas_8->Print("hist_jets_JVF.ps");


canvas_9->cd();
 hist_jets_MV1->Draw("colz");
 canvas_9->Print("hist_jets_MV1.ps");

canvas_10->cd();
 hist_MET->Draw("colz");
 canvas_8->Print("hist_MET.ps");

canvas_11->cd();
 hist_mTc->Draw("colz");
 canvas_8->Print("hist_transverse_mass.ps");

/*
canvas_12->cd();
 cutflow_MC->Draw("colz");
 canvas_12->Print("Cutflow_MC.ps");
 
 canvas_13->cd();
 lepton_t_isolation_MC->Draw("colz");
 canvas_13->Print("lep_Track_Isolation_MC.ps");

 canvas_14->cd();
 lepton_isolation_MC->Draw("colz");
 canvas_14->Print("lepton_Isolation_MC.ps");

 canvas_15->cd();
 lepton_eta_MC->Draw("colz");
 canvas_15->Print("lepton_eta_MC.ps");

canvas_16->cd();
 lepton_pt_MC->Draw("colz");
 canvas_16->Print("lepton_pt_MC.ps");


canvas_17->cd();
 hist_jets_pt_MC->Draw("colz");
 canvas_17->Print("hist_jets_pt_MC.ps");

canvas_18->cd();
 hist_jets_eta_MC->Draw("colz");
 canvas_18->Print("hist_jets_eta_MC.ps");

canvas_19->cd();
 hist_jets_JVF_MC->Draw("colz");
 canvas_19->Print("hist_jets_JVF_MC.ps");


canvas_20->cd();
 hist_jets_MV1_MC->Draw("colz");
 canvas_20->Print("hist_jets_MV1_MC.ps");

canvas_21->cd();
 hist_MET_MC->Draw("colz");
 canvas_21->Print("hist_MET_MC.ps");

canvas_22->cd();
 hist_mTc_MC->Draw("colz");
 canvas_22->Print("hist_transverse_mass_MC.ps");

*/


/*
hist_njets->Draw();
canvas->Draw();
*/    
    }
    

