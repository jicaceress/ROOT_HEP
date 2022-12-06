#include <iostream>
#include <string>
#include <stdio.h>


void TTbarSelection(){
    TFile *file_data = TFile::Open("/home/ivan/Documentos/Root_HEP/Data_8TeV.root");
    TFile *file = TFile::Open("/home/ivan/Documentos/Root_HEP/ttbar_8TeV.root");

    TTree *tree_data = (TTree*) file_data->Get("mini");
    TTree *tree = (TTree*) file->Get("mini");

    Bool_t e_trig;
    Bool_t mu_trig;
    Bool_t good_vtx;
    UInt_t lep_n;
    UInt_t jet_n;
    Float_t MET;
    Float_t MET_phi;

    Float_t lep_pt[10];
    Float_t lep_eta[10];
    Float_t lep_phi[10];
    Float_t lep_E[10];
    Int_t lep_type[10];
    Float_t lep_ptcone30[10];
    Float_t lep_etcone20[10];

    Float_t jet_pt[10];
    Float_t jet_eta[10];
    Float_t jet_jvf[10];
    Float_t jet_MV1[10];

	Float_t MCweight;
    Float_t scaleFactor_PILEUP;
	Float_t scaleFactor_ELE;
	Float_t scaleFactor_MUON;
	Float_t scaleFactor_BTAG;
	Float_t scaleFactor_TRIGGER;
	Float_t scaleFactor_JVFSF;
	Float_t scaleFactor_ZVERTEX;

	Bool_t e_trig_data;
	Bool_t mu_trig_data;
	Bool_t good_vtx_data;
	UInt_t lep_n_data;
	UInt_t jet_n_data;
	Float_t MET_data;
	Float_t MET_phi_data;

	Float_t lep_pt_data[10];
	Float_t lep_eta_data[10];
	Float_t lep_phi_data[10];
	Float_t lep_E_data[10];
	Int_t lep_type_data[10];
	Float_t lep_ptcone30_data[10];
	Float_t lep_etcone20_data[10];

	Float_t jet_pt_data[10];
	Float_t jet_eta_data[10];
	Float_t jet_jvf_data[10];
	Float_t jet_MV1_data[10];


    tree->SetBranchAddress("trigE", &e_trig);
    tree->SetBranchAddress("trigM", &mu_trig);
    tree->SetBranchAddress("hasGoodVertex", &good_vtx);
    tree->SetBranchAddress("lep_n", &lep_n);
    tree->SetBranchAddress("jet_n", &jet_n);
    tree->SetBranchAddress("met_et", &MET);
    tree->SetBranchAddress("met_phi", &MET_phi);

    tree->SetBranchAddress("lep_pt", &lep_pt);
    tree->SetBranchAddress("lep_eta", &lep_eta);
    tree->SetBranchAddress("lep_phi", &lep_phi);
    tree->SetBranchAddress("lep_E", &lep_E);
    tree->SetBranchAddress("lep_type", &lep_type);
    tree->SetBranchAddress("lep_ptcone30", &lep_ptcone30);
    tree->SetBranchAddress("lep_etcone20", &lep_etcone20);

    tree->SetBranchAddress("jet_pt", &jet_pt);
    tree->SetBranchAddress("jet_eta", &jet_eta);
    tree->SetBranchAddress("jet_jvf", &jet_jvf);
    tree->SetBranchAddress("jet_MV1", &jet_MV1);

	tree_data->SetBranchAddress("trigE", &e_trig_data);
	tree_data->SetBranchAddress("trigM", &mu_trig_data);
	tree_data->SetBranchAddress("hasGoodVertex", &good_vtx_data);
	tree_data->SetBranchAddress("lep_n", &lep_n_data);
	tree_data->SetBranchAddress("jet_n", &jet_n_data);
	tree_data->SetBranchAddress("met_et", &MET_data);
	tree_data->SetBranchAddress("met_phi", &MET_phi_data);

	tree_data->SetBranchAddress("lep_pt", &lep_pt_data);
	tree_data->SetBranchAddress("lep_eta", &lep_eta_data);
	tree_data->SetBranchAddress("lep_phi", &lep_phi_data);
	tree_data->SetBranchAddress("lep_E", &lep_E_data);
	tree_data->SetBranchAddress("lep_type", &lep_type_data);
	tree_data->SetBranchAddress("lep_ptcone30", &lep_ptcone30_data);
	tree_data->SetBranchAddress("lep_etcone20", &lep_etcone20_data);

	tree_data->SetBranchAddress("jet_pt", &jet_pt_data);
	tree_data->SetBranchAddress("jet_eta", &jet_eta_data);
	tree_data->SetBranchAddress("jet_jvf", &jet_jvf_data);
	tree_data->SetBranchAddress("jet_MV1", &jet_MV1_data);

	tree->SetBranchAddress("scaleFactor_PILEUP", &scaleFactor_PILEUP);
	tree->SetBranchAddress("scaleFactor_ELE", &scaleFactor_ELE);
	tree->SetBranchAddress("scaleFactor_MUON", &scaleFactor_MUON);
	tree->SetBranchAddress("scaleFactor_BTAG", &scaleFactor_BTAG);
	tree->SetBranchAddress("scaleFactor_TRIGGER", &scaleFactor_TRIGGER);
	tree->SetBranchAddress("scaleFactor_JVFSF", &scaleFactor_JVFSF);
	tree->SetBranchAddress("scaleFactor_ZVERTEX", &scaleFactor_ZVERTEX);
    tree->SetBranchAddress("mcWeight",&MCweight);


    TH1F *cutflow = new TH1F("Cutflow","Cutflow; Cut; Events",10,0,10);
    TH1F *hist_njets = new TH1F("Number of jets","n-jets; Jet multiplicity; Events",10,0,10);
    TH1F *hist_leppt =new TH1F("pT of Leptons","Leptons pT; pT[GeV] ; Events",100,0,100);
    TH1F *hist_trackisolation = new TH1F("Track isolation","Track Isolation; jetcone30pT/jetpT ; Event",100,0,0.5);
    TH1F *hist_calisolation = new TH1F("Calorimeter Isolation","Calorimeter Isolation ; jetcone20pT/jetpT ; Event",100,0,0.5);
    TH1F *hist_lepeta = new TH1F("Eta of Leptons","Leptons Eta ; Eta ; Event",100,-3,3);
    TH1F *hist_jetpt = new TH1F("pT of Jets","Jets pT ; pT [GeV] ; Event", 100, 0, 100);
    TH1F *hist_jeteta = new TH1F("Eta of Jets","Jets Eta ; eta ;Event",100,-3,3);
    TH1F *hist_jetJVF = new TH1F("JVF","JVF ; Jets JVF ;Event",100,0,1);
    TH1F *hist_jetMV1 = new TH1F("MV1","Jets MV1",100,0,1);
    TH1F *hist_nbjets = new TH1F("Number of b jets","Number of b-jets ; b-jets ; Events",6,0,6);
    TH1F *hist_MET = new TH1F("MET","Missing Transverse Energy ; MET[GeV] ; Events",200,0,200);
    TH1F *hist_mTW =new TH1F("mTW","Transverse missing energy of W ; MTW[GeV] ; Events",100,30,200);


    int nentries, nentriesMC, i;

    int cut1 = 0;
    int cut2 = 0;
    int cut3 = 0;
    int cut4 = 0;
    int cut5 = 0;
    int cut6 = 0;
    int cut7 = 0;
    int cut8 = 0;
    nentries = (Int_t)tree_data->GetEntries();
    for (i = 0; i < nentries; i++)
    {
        tree_data->GetEntry(i);
        //First cut: Good vertex
        if(!good_vtx_data) continue;
        cut1++;
        cutflow->Fill(1);

        //Second cut: Trigger
        if(!e_trig_data && !mu_trig_data) continue;
        cut2++;
        cutflow->Fill(2);

        // Preselection of good leptons
        int n_mu=0;
        int n_el=0;
        int nlep=0;
        int good_lep;

        //Loop over leptons
        for(unsigned int l=0; l<lep_n_data; l++){
            if( lep_pt_data[l] < 25000.) continue;
            if( lep_ptcone30_data[l]/lep_pt_data[l] > 0.15 ) continue;
            if( lep_etcone20_data[l]/lep_pt_data[l] > 0.15 ) continue;
            if( lep_type_data [l]==13 && TMath::Abs(lep_eta_data[l]) < 2.5){
                n_mu++; good_lep=l;
                }
            if(lep_type_data[l]==11){
                if(TMath::Abs(lep_eta_data[l])<1.37 || (1.52<TMath::Abs(lep_eta_data[l])<2.47)){
                    n_el++;good_lep=l;
                }
            }
        }
        nlep=n_mu+n_el;
        //Select events with only 1 good lepton and fill the cutflow histogram
        //Example:
        //Third cut (one good lepton):
        if(nlep!=1) continue;
        cutflow->Fill(3);
        cut3++;


        int n_jets=0;
        int n_bjets=0;

        //Fourth cut: At least 4 jets
        if(jet_n_data<4) continue;
        cutflow->Fill(4);
        cut4++;


        for(unsigned int j=0; j<jet_n_data; j++){
            // To complete: apply jet cuts to find the good jets
            if(jet_pt_data[j] < 25000.) continue;
            if(TMath::Abs(jet_eta_data[j])>2.5) continue;
            if(jet_pt_data[j]<50000. && TMath::Abs(jet_eta_data[j])<2.4 && jet_jvf_data[j]<0.5) continue;
            if(jet_MV1_data[j]>=0.7892){
                n_bjets++;
                }
            n_jets++;
            }
        //Fifth cut: At least 4 good jets
        if(n_jets<4) continue;
        cutflow->Fill(5);
        cut5++;

        //Sixth cut: at least two b-jet
        if(n_bjets<=1) continue;
        cutflow->Fill(6);
        cut6++;


        //Seventh cut: MET > 30 GeV
        if(MET_data<30000.) continue;
        cutflow->Fill(7);
        cut7++;

        // TLorentzVector definitions
        // TLorentzVector Lepton  = TLorentzVector();
        // TLorentzVector  MeT  = TLorentzVector();

        // Lepton.SetPtEtaPhiE(lep_pt_data[good_lep],)
        //To complete: Lorentz vector for the lepton and MET. Use SetPtEtaPhiE().

        //Calculation of the mTW using TLorentz vectors
        float mTW_data = sqrt(2*lep_pt_data[good_lep]*MET_data*(1-cos(MET_phi_data-lep_phi_data[good_lep])));

        //Eight cut: mTW > 30 GeV
        if(mTW_data< 30000.) continue;
        cutflow->Fill(8);
        cut8++;

        hist_leppt->Fill(lep_pt_data[good_lep]/1000.,1.);
        hist_lepeta->Fill(lep_eta_data[good_lep],1.);
        hist_trackisolation->Fill(lep_ptcone30_data[good_lep]/lep_pt_data[good_lep],1.);
        hist_calisolation->Fill(lep_etcone20_data[good_lep]/lep_pt_data[good_lep],1.);

        hist_jetpt->Fill(jet_pt_data[0]/1000.,1.);
        hist_jeteta->Fill(jet_eta_data[0],1.);
        hist_jetJVF->Fill(jet_jvf_data[0],1.);
        hist_jetMV1->Fill(jet_MV1_data[0],1.);

        hist_njets->Fill(n_jets,1.);
        hist_nbjets->Fill(n_bjets,1.);

        hist_MET->Fill(MET_data/1000.,1.);
        hist_mTW->Fill(mTW_data/1000.,1.);
    }

    std::cout << "\nTTbar selection for data\n" << std::endl;
    std::cout << "All events:" << nentries << std::endl;
    std::cout << "Cut1:" << cut1 << std::endl;
    std::cout << "Cut2:" << cut2 << std::endl;
    std::cout << "Cut3:" << cut3 << std::endl;
    std::cout << "Cut4:" << cut4 << std::endl;
    std::cout << "Cut5:" << cut5 << std::endl;
    std::cout << "Cut6:" << cut6 << std::endl;
    std::cout << "Cut7:" << cut7 << std::endl;
    std::cout << "Cut8:" << cut8 << std::endl;

    TCanvas *canvas_cutflow =new TCanvas("Canvas cutflow","a");
    gStyle->SetOptStat(0);
    canvas_cutflow->SetLogy();
    cutflow->Draw("BOX");
    canvas_cutflow->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/cutflow.png");
    canvas_cutflow->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_njets = new TCanvas("Canvas njets","b");
    gStyle->SetOptStat(0);
    canvas_njets->SetLogy();
    hist_njets->Draw();
    canvas_njets->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/nJets.png");
    canvas_njets->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_leppt =new TCanvas("Canvas lep pT","c");
    gStyle->SetOptStat(0);
    hist_leppt->Draw();
    canvas_leppt->Draw();
    canvas_leppt->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/LeptonpT.png");
    canvas_leppt->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_trackisolation = new TCanvas("Canvas TrackIsolation","d");
	gStyle->SetOptStat(0);
    hist_trackisolation->Draw();
    canvas_trackisolation->Draw();
    canvas_trackisolation->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/TrackIsolation.png");
    canvas_trackisolation->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_calisolation = new TCanvas("Canvas CalIsolation","e");
	gStyle->SetOptStat(0);
    canvas_calisolation->SetLogy();
    hist_calisolation->Draw();
    canvas_calisolation->Draw();
    canvas_calisolation->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/CalIsolation.png");
    canvas_calisolation->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_lepeta = new TCanvas("Canvas Lepeta","f");
	gStyle->SetOptStat(0);
    hist_lepeta->Draw();
    canvas_lepeta->Draw();
    canvas_lepeta->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/LeptonEta.png");
    canvas_lepeta->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_jetpt = new TCanvas("Canvas jetpT","g");
	gStyle->SetOptStat(0);
    hist_jetpt->Draw();
    canvas_jetpt->Draw();
    canvas_jetpt->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/JetpT.png");
    canvas_jetpt->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_jeteta = new TCanvas("Canvas jeteta","h");
	gStyle->SetOptStat(0);
    hist_jeteta->Draw();
    canvas_jeteta->Draw();
    canvas_jeteta->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/JetEta.png");
    canvas_jeteta->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_jetJVF = new TCanvas("Canvas JVF","i");
	gStyle->SetOptStat(0);
    hist_jetJVF->Draw();
    canvas_jetJVF->Draw();
    canvas_jetJVF->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/JetJVF.png");
    canvas_jetJVF->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_jetMV1 = new TCanvas("Canvas MV1","j");
	gStyle->SetOptStat(0);
    hist_jetMV1->Draw();
    canvas_jetMV1->Draw();
    canvas_jetMV1->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/JetMV1.png");
    canvas_jetMV1->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_nbjets = new TCanvas("Canvas nbjets","k");
    canvas_nbjets->SetLogy();
	gStyle->SetOptStat(0);
    hist_nbjets->Draw();
    canvas_nbjets->Draw();
    canvas_nbjets->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/nbjets.png");
    canvas_nbjets->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_MET = new TCanvas("Canvas MET","l");
	gStyle->SetOptStat(0);
    hist_MET->Draw();
    canvas_MET->Draw();
    canvas_MET->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/MET.png");
    canvas_MET->Close(); 
    gSystem->ProcessEvents(); 


    TCanvas *canvas_mTW = new TCanvas("Canvas mTW","m");
    hist_mTW->Draw();
    canvas_mTW->Draw();
    canvas_mTW->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Data_selection_after_cuts/mTW.png");
    canvas_mTW->Close(); 
    gSystem->ProcessEvents(); 


    //Reset the cuts and the histograms
    cut1=cut2=cut3=cut4=cut5=cut6=cut7=cut8=0.;
    TH1F *cutflowMC = new TH1F("Cutflow","Cutflow; Cut; Events",10,0,10);
    TH1F *hist_njetsMC = new TH1F("Number of jets","n-jets; Jet multiplicity; Events",10,0,10);
    TH1F *hist_lepptMC =new TH1F("pT of Leptons","Leptons pT; pT[GeV] ; Events",100,0,100);
    TH1F *hist_trackisolationMC = new TH1F("Track isolation","Track Isolation; jetcone30pT/jetpT ; Event",100,0,0.5);
    TH1F *hist_calisolationMC = new TH1F("Calorimeter Isolation","Calorimeter Isolation ; jetcone20pT/jetpT ; Event",100,0,0.5);
    TH1F *hist_lepetaMC = new TH1F("Eta of Leptons","Leptons Eta ; Eta ; Event",100,-3,3);
    TH1F *hist_jetptMC = new TH1F("pT of Jets","Jets pT ; pT [GeV] ; Event", 100, 0, 100);
    TH1F *hist_jetetaMC = new TH1F("Eta of Jets","Jets Eta ; eta ;Event",100,-3,3);
    TH1F *hist_jetJVFMC = new TH1F("JVF","JVF ; Jets JVF ;Event",100,0,1);
    TH1F *hist_jetMV1MC = new TH1F("MV1","Jets MV1",100,0,1);
    TH1F *hist_nbjetsMC= new TH1F("Number of b jets","Number of b-jets ; b-jets ; Events",6,0,6);
    TH1F *hist_METMC = new TH1F("MET","Missing Transverse Energy ; MET[GeV] ; Events",200,0,200);
    TH1F *hist_mTWMC =new TH1F("mTW","Transverse missing energy of W ; MTW[GeV] ; Events",100,30,200);





    nentriesMC = (Int_t)tree->GetEntries();
    for (i = 0; i < nentriesMC; i++)
    {
        tree->GetEntry(i);   
        Float_t SF = scaleFactor_PILEUP * scaleFactor_ELE * scaleFactor_MUON * scaleFactor_BTAG * scaleFactor_TRIGGER * scaleFactor_JVFSF * scaleFactor_ZVERTEX;
        Float_t ScaleWeight = 1000./(49761200.21*0.072212854/137.29749);
        Float_t evtw = SF * MCweight*ScaleWeight;

        //First cut: Good vertex
        if(!good_vtx) continue;
        cut1++;
        cutflowMC->Fill(1);

        //Second cut: Trigger
        if(!e_trig && !mu_trig) continue;
        cut2++;
        cutflowMC->Fill(2);

        // Preselection of good leptons
        int n_mu=0;
        int n_el=0;
        int nlep=0;
        int good_lep;

        //Loop over leptons
        for(unsigned int l=0; l<lep_n; l++){
            if( lep_pt[l] < 25000.) continue; 
            if( lep_ptcone30[l]/lep_pt[l] > 0.15 ) continue; 
            if( lep_etcone20[l]/lep_pt[l] > 0.15 ) continue;  
            if( lep_type [l]==13 && TMath::Abs(lep_eta[l]) < 2.5){
                n_mu++; good_lep=l;
                }
            if(lep_type[l]==11){
                if(TMath::Abs(lep_eta[l])<1.37 || (1.52<TMath::Abs(lep_eta[l])<2.47)){
                    n_el++;good_lep=l;
                }
            }
        }
        nlep=n_mu+n_el;
        //Select events with only 1 good lepton and fill the cutflow histogram 
        //Example:
        //Third cut (one good lepton):
        if(nlep!=1) continue;
        cutflowMC->Fill(3); 
        cut3++;
        
        
        int n_jets=0;
        int n_bjets=0;
        
        //Fourth cut: At least 4 jets
        if(jet_n<4) continue;
        cutflowMC->Fill(4); 
        cut4++;

    
        for(unsigned int j=0; j<jet_n; j++){
            // To complete: apply jet cuts to find the good jets
            if(jet_pt[j] < 25000.) continue;
            if(TMath::Abs(jet_eta[j])>2.5) continue;
            if(jet_pt[j]<50000. && TMath::Abs(jet_eta[j])<2.4 && jet_jvf[j]<0.5) continue;                                                                                    
            if(jet_MV1[j]>=0.7892){ 
                n_bjets++;
                }
            n_jets++;
            }
        //Fifth cut: At least 4 good jets
        if(n_jets<4) continue; 
        cutflowMC->Fill(5); 
        cut5++;
        

        //Sixth cut: at least two b-jet
        if(n_bjets<=1) continue;
        cutflowMC->Fill(6);
        cut6++;



        //Seventh cut: MET > 30 GeV
        if(MET<30000.) continue;
        cutflowMC->Fill(7);
        cut7++;
        
        // TLorentzVector definitions                                                               
        // TLorentzVector Lepton  = TLorentzVector();
        // TLorentzVector  MeT  = TLorentzVector();

        // Lepton.SetPtEtaPhiE(lep_pt[good_lep],)
        //To complete: Lorentz vector for the lepton and MET. Use SetPtEtaPhiE().

        //Calculation of the mTW using TLorentz vectors             
        float mTW = sqrt(2*lep_pt[good_lep]*MET*(1-cos(MET_phi-lep_phi[good_lep])));
        //Eight cut: mTW > 30 GeV
        if(mTW < 30000.) continue;
        cutflowMC->Fill(8);
        cut8++;

        hist_lepptMC->Fill(lep_pt[good_lep]/1000.,evtw);
        hist_lepetaMC->Fill(lep_eta[good_lep],evtw);
        hist_trackisolationMC->Fill(lep_ptcone30[good_lep]/lep_pt[good_lep],evtw);
        hist_calisolationMC->Fill(lep_etcone20[good_lep]/lep_pt[good_lep],evtw);

        hist_njetsMC->Fill(n_jets,evtw);
        hist_nbjetsMC->Fill(n_bjets,evtw);

        hist_jetptMC->Fill(jet_pt[0]/1000.,evtw);
        hist_jetetaMC->Fill(jet_eta[0],evtw);
        hist_jetJVFMC->Fill(jet_jvf[0],evtw);
        hist_jetMV1MC->Fill(jet_MV1[0],evtw);

        hist_METMC->Fill(MET/1000.);
        hist_mTWMC->Fill(mTW/1000.);
    }

    std::cout << "\nTTbar selection for Montecarlo Simulation:\n" << std::endl;
    std::cout << "All events:" << nentries << std::endl;
    std::cout << "Cut1:" << cut1 << std::endl;
    std::cout << "Cut2:" << cut2 << std::endl;
    std::cout << "Cut3:" << cut3 << std::endl;
    std::cout << "Cut4:" << cut4 << std::endl;
    std::cout << "Cut5:" << cut5 << std::endl;
    std::cout << "Cut6:" << cut6 << std::endl;
    std::cout << "Cut7:" << cut7 << std::endl;
    std::cout << "Cut8:" << cut8 << std::endl;

    TCanvas *canvas_cutflowMC =new TCanvas("Canvas cutflow","a");
    gStyle->SetOptStat(0);
    canvas_cutflowMC->SetLogy();
    cutflowMC->Draw(",h");
    canvas_cutflowMC->Draw();
    canvas_cutflowMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/cutflow.png");
    canvas_cutflowMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_njetsMC = new TCanvas("Canvas njets","b");
    gStyle->SetOptStat(0);
    canvas_njetsMC->SetLogy();
    hist_njetsMC->Draw(",h");
    canvas_njetsMC->Draw();
    canvas_njetsMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/nJets.png");
    canvas_njetsMC->Close(); 
    gSystem->ProcessEvents(); 
    
    TCanvas *canvas_lepptMC =new TCanvas("Canvas lep pT","c");
    gStyle->SetOptStat(0);
    hist_lepptMC->Draw(",h");
    canvas_lepptMC->Draw();
    canvas_lepptMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/LeptonpT.png");
    canvas_lepptMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_trackisolationMC = new TCanvas("Canvas TrackIsolation","d");
	gStyle->SetOptStat(0);
    hist_trackisolationMC->Draw(",h");
    canvas_trackisolationMC->Draw();
    canvas_trackisolationMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/TrackIsolation.png");
    canvas_trackisolationMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_calisolationMC = new TCanvas("Canvas CalIsolation","e");
	gStyle->SetOptStat(0);
    canvas_calisolationMC->SetLogy();
    hist_calisolationMC->Draw(",h");
    canvas_calisolationMC->Draw();
    canvas_calisolationMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/CalIsolation.png");
    canvas_calisolationMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_lepetaMC = new TCanvas("Canvas Lepeta","f");
	gStyle->SetOptStat(0);
    hist_lepetaMC->Draw(",h");
    canvas_lepetaMC->Draw();
    canvas_lepetaMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/LeptonEta.png");
    canvas_lepetaMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_jetptMC = new TCanvas("Canvas jetpT","g");
	gStyle->SetOptStat(0);
    hist_jetptMC->Draw(",h");
    canvas_jetptMC->Draw();
    canvas_jetptMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/JetpT.png");
    canvas_jetptMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_jetetaMC = new TCanvas("Canvas jeteta","h");
	gStyle->SetOptStat(0);
    hist_jetetaMC->Draw(",h");
    canvas_jetetaMC->Draw();
    canvas_jetetaMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/JetEta.png");
    canvas_jetetaMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_jetJVFMC = new TCanvas("Canvas JVF","i");
	gStyle->SetOptStat(0);
    hist_jetJVFMC->Draw(",h");
    canvas_jetJVFMC->Draw();
    canvas_jetJVFMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/JetJVF.png");
    canvas_jetJVFMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_jetMV1MC = new TCanvas("Canvas MV1","j");
	gStyle->SetOptStat(0);
    hist_jetMV1MC->Draw(",h");
    canvas_jetMV1MC->Draw();
    canvas_jetMV1MC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/JetMV1.png");
    canvas_jetMV1MC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_nbjetsMC= new TCanvas("Canvas nbjets","k");
    canvas_nbjetsMC->SetLogy();
	gStyle->SetOptStat(0);
    hist_nbjetsMC->Draw(",h");
    canvas_nbjetsMC->Draw();
    canvas_nbjetsMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/nbjets.png");
    canvas_nbjetsMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_METMC= new TCanvas("Canvas MET","l");
	gStyle->SetOptStat(0);
    hist_METMC->Draw(",h");
    canvas_METMC->Draw();
    canvas_METMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/MET.png");
    canvas_METMC->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_mTWMC= new TCanvas("Canvas mTW","m");
    hist_mTWMC->Draw(",h");
    canvas_mTWMC->Draw();
    canvas_mTWMC->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Montecarlo_selection_after_cuts_without_weight/mTW.png");
    canvas_mTWMC->Close(); 
    gSystem->ProcessEvents(); 

    //Comparision histograms

    TCanvas *canvas_cutflowComparision =new TCanvas("Canvas cutflow","a");
    gStyle->SetOptStat(0);
    canvas_cutflowComparision->SetLogy();
    cutflow->SetMarkerStyle(21);
    cutflow->SetMarkerColor(kBlack);
    cutflow->Draw("P");
    cutflowMC->SetMarkerStyle(34);
    cutflowMC->SetMarkerColor(kRed);
    cutflowMC->Draw("P,same");
	auto legend_cutflow_comparison = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_cutflow_comparison->SetBorderSize(0);
	legend_cutflow_comparison->AddEntry(cutflow, "Data", "lep");
	legend_cutflow_comparison->AddEntry(cutflowMC, "MC", "lep");
	legend_cutflow_comparison->Draw();
    canvas_cutflowComparision->Draw();
    canvas_cutflowComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/cutflow.png");
    canvas_cutflowComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_njetsComparision = new TCanvas("Canvas njets","b");
    gStyle->SetOptStat(0);
    canvas_njetsComparision->SetLogy();
    hist_njets->SetMarkerStyle(21);
    hist_njets->SetMarkerColor(kBlack);
    hist_njets->Draw("P");
    hist_njetsMC->SetMarkerStyle(34);
    hist_njetsMC->SetMarkerColor(kRed);
    hist_njetsMC->Draw("P,same");
	auto legend_njets = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_njets->SetBorderSize(0);
	legend_njets->AddEntry(hist_njets, "Data", "lep");
	legend_njets->AddEntry(hist_njetsMC, "MC", "lep");
	legend_njets->Draw();
    canvas_njetsComparision->Draw();
    canvas_njetsComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/nJets.png");
    canvas_njetsComparision->Close(); 
    gSystem->ProcessEvents(); 
    
    TCanvas *canvas_lepptComparision =new TCanvas("Canvas lep pT","c");
    gStyle->SetOptStat(0);
    hist_leppt->SetMarkerStyle(21);
    hist_leppt->SetMarkerColor(kBlack);
    hist_leppt->Draw("P");
    hist_lepptMC->SetMarkerStyle(34);
    hist_lepptMC->SetMarkerColor(kRed);
    hist_lepptMC->Draw("P,same");
    auto legend_leppt = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_leppt->SetBorderSize(0);
	legend_leppt->AddEntry(hist_leppt, "Data", "lep");
	legend_leppt->AddEntry(hist_lepptMC, "MC", "lep");
	legend_leppt->Draw();
    canvas_lepptComparision->Draw();
    canvas_lepptComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/LeptonpT.png");
    canvas_lepptComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_trackComparision = new TCanvas("Canvas TrackIsolation","d");
	gStyle->SetOptStat(0);
    hist_trackisolation->SetMarkerStyle(21);
    hist_trackisolation->SetMarkerColor(kBlack);
    hist_trackisolation->Draw("P");
    hist_trackisolationMC->SetMarkerStyle(34);
    hist_trackisolationMC->SetMarkerColor(kRed);
    hist_trackisolationMC->Draw("P,same");
    auto legend_trackiso = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_trackiso->SetBorderSize(0);
	legend_trackiso->AddEntry(hist_trackisolation, "Data", "lep");
	legend_trackiso->AddEntry(hist_trackisolationMC, "MC", "lep");
    legend_trackiso->Draw();
    canvas_trackComparision->Draw();
    canvas_trackComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/TrackIsolation.png");
    canvas_trackComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_calComparision = new TCanvas("Canvas CalIsolation","e");
	gStyle->SetOptStat(0);
    canvas_calComparision->SetLogy();
    hist_calisolation->SetMarkerStyle(21);
    hist_calisolation->SetMarkerColor(kBlack);
    hist_calisolation->Draw("P");
    hist_calisolationMC->SetMarkerStyle(34);
    hist_calisolationMC->SetMarkerColor(kRed);
    hist_calisolationMC->Draw("P,same");
    auto legend_caliso = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_caliso->SetBorderSize(0);
	legend_caliso->AddEntry(hist_calisolation, "Data", "lep");
	legend_caliso->AddEntry(hist_calisolationMC, "MC", "lep");
    legend_caliso->Draw();
    canvas_calComparision->Draw();
    canvas_calComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/CalIsolation.png");
    canvas_calComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_lepetaComparision = new TCanvas("Canvas Lepeta","f");
	gStyle->SetOptStat(0);
    hist_lepeta->SetMarkerStyle(21);
    hist_lepeta->SetMarkerColor(kBlack);
    hist_lepeta->Draw("P");
    hist_lepetaMC->SetMarkerStyle(34);
    hist_lepetaMC->SetMarkerColor(kRed);
    hist_lepetaMC->Draw("P,same");
    auto legend_lepeta = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_lepeta->SetBorderSize(0);
	legend_lepeta->AddEntry(hist_lepeta, "Data", "lep");
	legend_lepeta->AddEntry(hist_lepetaMC, "MC", "lep");
    legend_lepeta->Draw();
    canvas_lepetaComparision->Draw();
    canvas_lepetaComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/LeptonEta.png");
    canvas_lepetaComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_jetptComparision = new TCanvas("Canvas jetpT","g");
	gStyle->SetOptStat(0);
    hist_jetpt->SetMarkerStyle(21);
    hist_jetpt->SetMarkerColor(kBlack);
    hist_jetpt->Draw("P");
    hist_jetptMC->SetMarkerStyle(34);
    hist_jetptMC->SetMarkerColor(kRed);
    hist_jetptMC->Draw("P,same");
    auto legend_jetpt = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_jetpt->SetBorderSize(0);
	legend_jetpt->AddEntry(hist_jetpt, "Data", "lep");
	legend_jetpt->AddEntry(hist_jetptMC, "MC", "lep");
    legend_jetpt->Draw();
    canvas_jetptComparision->Draw();
    canvas_jetptComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/JetpT.png");
    canvas_jetptComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_jetetaComparision = new TCanvas("Canvas jeteta","h");
	gStyle->SetOptStat(0);
    hist_jeteta->SetMarkerStyle(21);
    hist_jeteta->SetMarkerColor(kBlack);
    hist_jeteta->Draw("P");
    hist_jetetaMC->SetMarkerStyle(34);
    hist_jetetaMC->SetMarkerColor(kRed);
    hist_jetetaMC->Draw("P,same");
    auto legend_jeteta = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_jeteta->SetBorderSize(0);
	legend_jeteta->AddEntry(hist_jeteta, "Data", "lep");
	legend_jeteta->AddEntry(hist_jetetaMC, "MC", "lep");
    legend_jeteta->Draw();
    canvas_jetetaComparision->Draw();
    canvas_jetetaComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/JetEta.png");
    canvas_jetetaComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_jetJVFComparision = new TCanvas("Canvas JVF","i");
	gStyle->SetOptStat(0);
    hist_jetJVF->SetMarkerStyle(21);
    hist_jetJVF->SetMarkerColor(kBlack);
    hist_jetJVF->Draw("P");
    hist_jetJVFMC->SetMarkerStyle(34);
    hist_jetJVFMC->SetMarkerColor(kRed);
    hist_jetJVFMC->Draw("P,same");
    auto legend_jetJVF= new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_jetJVF->SetBorderSize(0);
	legend_jetJVF->AddEntry(hist_jetJVF, "Data", "lep");
	legend_jetJVF->AddEntry(hist_jetJVFMC, "MC", "lep");
    legend_jetJVF->Draw();
    canvas_jetJVFComparision->Draw();
    canvas_jetJVFComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/JetJVF.png");
    canvas_jetJVFComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_jetMV1Comparision = new TCanvas("Canvas MV1","j");
	gStyle->SetOptStat(0);
    hist_jetMV1->SetMarkerStyle(21);
    hist_jetMV1->SetMarkerColor(kBlack);
    hist_jetMV1->Draw("P");
    hist_jetMV1MC->SetMarkerStyle(34);
    hist_jetMV1MC->SetMarkerColor(kRed);
    hist_jetMV1MC->Draw("P,same");
    auto legend_jetMV1 = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_jetMV1->SetBorderSize(0);
	legend_jetMV1->AddEntry(hist_jetMV1, "Data", "lep");
	legend_jetMV1->AddEntry(hist_jetMV1MC, "MC", "lep");
    legend_jetMV1->Draw();
    canvas_jetMV1Comparision->Draw();
    canvas_jetMV1Comparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/JetMV1.png");
    canvas_jetMV1Comparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_nbjetsComparision= new TCanvas("Canvas nbjets","k");
    canvas_nbjetsComparision->SetLogy();
	gStyle->SetOptStat(0);
    hist_nbjets->SetMarkerStyle(21);
    hist_nbjets->SetMarkerColor(kBlack);
    hist_nbjets->Draw("P");
    hist_nbjetsMC->SetMarkerStyle(34);
    hist_nbjetsMC->SetMarkerColor(kRed);
    hist_nbjetsMC->Draw("P,same");
    auto legend_nbjets = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_nbjets->SetBorderSize(0);
	legend_nbjets->AddEntry(hist_nbjets, "Data", "lep");
	legend_nbjets->AddEntry(hist_nbjetsMC, "MC", "lep");
    legend_nbjets->Draw();
    canvas_nbjetsComparision->Draw();
    canvas_nbjetsComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/nbjets.png");
    canvas_nbjetsComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_METComparision= new TCanvas("Canvas MET","l");
	gStyle->SetOptStat(0);
    hist_MET->SetMarkerStyle(21);
    hist_MET->SetMarkerColor(kBlack);
    hist_MET->Draw("P");
    hist_METMC->SetMarkerStyle(34);
    hist_METMC->SetMarkerColor(kRed);
    hist_METMC->Draw("P,same");
    auto legend_MET = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_MET->SetBorderSize(0);
	legend_MET->AddEntry(hist_MET, "Data", "lep");
	legend_MET->AddEntry(hist_METMC, "MC", "lep");
    legend_MET->Draw();
    canvas_METComparision->Draw();
    canvas_METComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/MET.png");
    canvas_METComparision->Close(); 
    gSystem->ProcessEvents(); 

    TCanvas *canvas_mTWComparision= new TCanvas("Canvas mTW","m");
    hist_mTW->SetMarkerStyle(21);
    hist_mTW->SetMarkerColor(kBlack);
    hist_mTW->Draw("P");
    hist_mTWMC->SetMarkerStyle(34);
    hist_mTWMC->SetMarkerColor(kRed);
    hist_mTWMC->Draw("P,same");
    auto legend_mTW = new TLegend(0.7, 0.6, 0.89, 0.5);
	legend_mTW->SetBorderSize(0);
	legend_mTW->AddEntry(hist_mTW, "Data", "lep");
	legend_mTW->AddEntry(hist_mTWMC, "MC", "lep");
    legend_mTW->Draw();
    canvas_mTWComparision->Draw();
    canvas_mTWComparision->SaveAs("/home/ivan/Documentos/Root_HEP/ROOT_HEP/Imagenes/Ejercicio3y4/Comparision/mTW.png");
    canvas_mTWComparision->Close(); 
    gSystem->ProcessEvents(); 
}
