#define TopAnalysis_cxx
// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.

#include "WZAnalysis.h"
#include "WZhistograms.h"
#include <iostream>
#include <cstring>
#include <string>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLorentzVector.h>

string name;

void WZAnalysis::Begin(TTree *)
{
	// The Begin() function is called at the start of the query.
	// When running with PROOF Begin() is only called on the client.
	// The tree argument is deprecated (on PROOF 0 is passed).
}

void WZAnalysis::SlaveBegin(TTree *)
{
	// The SlaveBegin() function is called after the Begin() function.
	// When running with PROOF SlaveBegin() is called on each slave server.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();
	printf("Starting analysis with process option: %s \n", option.Data());

	name = option;

	define_histograms();

	FillOutputList();
}

Bool_t WZAnalysis::Process(Long64_t entry)
{
	// The Process() function is called for each entry in the tree (or possibly
	// keyed object in the case of PROOF) to be processed. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	// When processing keyed objects with PROOF, the object is already loaded
	// and is available via the fObject pointer.
	//
	// This function should contain the \"body\" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	//
	// Use fStatus to set the return value of TTree::Process().
	//
	// The return value is currently not used.

	fChain->GetTree()->GetEntry(entry);
	//  int cut1_mc = 0;

	if (fChain->GetTree()->GetEntries() > 0)
	{
		// Do analysis

		// SF
		Float_t scaleFactor = scaleFactor_ELE * scaleFactor_MUON * scaleFactor_TRIGGER;
		// EventW
		Float_t eventWeight = mcWeight * scaleFactor_PILEUP * scaleFactor_ZVERTEX;
		// weight = SF * EventW
		Double_t weight = scaleFactor * eventWeight;

		// Make difference between data and MC
		if (weight == 0. && channelNumber != 110090 && channelNumber != 110091)
			weight = 1.;

		// Missing Et of the event in GeV
		Float_t missingEt = met_et / 1000.;

		// Preselection cut for electron/muon trigger, Good Run List, and good vertex
		if (trigE || trigM)
		{
			if (passGRL)
			{
				if (hasGoodVertex)
				{
					// Find the good leptons
					int goodlep_index[3]; // list for good leptons, we need exactly three good ones.
					int goodlep_n = 0;	  // Number for the good lepton index
					int lep_index = 0;	  // number for the lepton index

					for (unsigned int i = 0; i < lep_n; i++)
					{
						if (lep_pt->at(i) > 25000. && (lep_ptcone30[i] / lep_pt->at(i)) < 0.15 && (lep_etcone20[i] / lep_pt->at(i)) < 0.15) // isolated track
						{
							// electron selection in fiducial region excluding candidates in the transition region between the barrel and endcap electromagnetic calorimeters
							if (lep_type->at(i) == 11 && TMath::Abs(lep_eta->at(i)) < 2.47 && (TMath::Abs(lep_eta->at(i)) < 1.37 || TMath::Abs(lep_eta->at(i)) > 1.52))
							{
								goodlep_n++;
								goodlep_index[lep_index] = i;
								lep_index++;
							}
							if (lep_type->at(i) == 13 && TMath::Abs(lep_eta->at(i)) < 2.5)
							{ // muon selection
								goodlep_n++;
								goodlep_index[lep_index] = i;
								lep_index++;
							}
						}
					}

					// Exactly three good leptons
					if (goodlep_n == 3)
					{
						// We take exactly three good leptons with the good lepton index and we set the transverse mass condition.
						int goodlep1_index = goodlep_index[0];
						int goodlep2_index = goodlep_index[1];
						int goodlep3_index = goodlep_index[2];

						// TLorentzVector definitions
						TLorentzVector Lepton_1 = TLorentzVector(0., 0., 0., 0.);
						TLorentzVector Lepton_2 = TLorentzVector(0., 0., 0., 0.);
						TLorentzVector Lepton_3 = TLorentzVector(0., 0., 0., 0.);

						TLorentzVector MeT = TLorentzVector(0., 0., 0., 0.);
						TLorentzVector Lepton1_MeT = TLorentzVector(0., 0., 0., 0.);

						Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep1_index), lep_eta->at(goodlep1_index), lep_phi->at(goodlep1_index), lep_E->at(goodlep1_index));
						Lepton_2.SetPtEtaPhiE(lep_pt->at(goodlep2_index), lep_eta->at(goodlep2_index), lep_phi->at(goodlep2_index), lep_E->at(goodlep2_index));
						Lepton_3.SetPtEtaPhiE(lep_pt->at(goodlep3_index), lep_eta->at(goodlep3_index), lep_phi->at(goodlep3_index), lep_E->at(goodlep3_index));

						MeT.SetPtEtaPhiE(met_et, 0, met_phi, met_et);

						// Now we defirne a TLV with a sum of the three leptons we have found with the cuts. This is to find the invariant mass of the two
						// and then procure a cut for it, with the invariant mass of Z.
						TLorentzVector Lepton_1plus2 = TLorentzVector(0., 0., 0., 0.);
						TLorentzVector Lepton_1plus3 = TLorentzVector(0., 0., 0., 0.);
						TLorentzVector Lepton_2plus3 = TLorentzVector(0., 0., 0., 0.);

						Lepton_1plus2 = Lepton_1 + Lepton_2;
						Lepton_1plus3 = Lepton_1 + Lepton_3;
						Lepton_2plus3 = Lepton_2 + Lepton_3;

						// WE now find the invariant mass

						float InvMass_lep12 = Lepton_1plus2.Mag() / 100.;
						float InvMass_lep13 = Lepton_1plus3.Mag() / 100.;
						float InvMass_lep23 = Lepton_2plus3.Mag() / 100.;

						// After this we would apply the cuts

						// we define a variable such that the difference of it with the invariant mass of the Z lepton is t cuton the differen
						// this would be one of the final cuts, but it will be needed

						float delta_lep12 = 0;
						float delta_lep13 = 0;
						float delta_lep23 = 0;

						// now, since we need a SINGLE electron/muon on the endstate, we need to check the charge of the leptons
						// since we are searchig for two with the opposite charge and one charge one, so the end state has
						// one unique single ele/muon
						//
						if (lep_charge->at(goodlep1_index) * lep_charge->at(goodlep2_index) < 0) // check for the charge
						{
							if (lep_type->at(goodlep1_index) == lep_type->at(goodlep2_index)) // Make sure leptons are the same
							{
								delta_lep12 = TMath::Abs(InvMass_lep12 - 91.18); // do the cut
							}
						}
						// and we repeat
						if (lep_charge->at(goodlep1_index) * lep_charge->at(goodlep3_index) < 0) // check for the charge
						{
							if (lep_type->at(goodlep1_index) == lep_type->at(goodlep3_index)) // Make sure leptons are the same
							{
								delta_lep13 = TMath::Abs(InvMass_lep13 - 91.18); // do the cut
							}
						}

						if (lep_charge->at(goodlep2_index) * lep_charge->at(goodlep3_index) < 0) // check for the charge
						{
							if (lep_type->at(goodlep2_index) == lep_type->at(goodlep3_index)) // Make sure leptons are the same
							{
								delta_lep23 = TMath::Abs(InvMass_lep13 - 91.18); // do the cut
							}
						}

						int caso = 0;
						float cut_Ztran = 0; // this is to store the value of the transverse mass difference for any isntance

						if (delta_lep12 < 0 && delta_lep23 < 0 && delta_lep13 > 0)
						{ // case inv mass 13, excludes 23 and 12
							cut_Ztran = delta_lep13;
							caso = 2;
						}

						if (delta_lep12 < 0 && delta_lep13 < 0 && delta_lep23 > 0)
						{ // repeat
							cut_Ztran = delta_lep23;
							caso = 1;
						}

						if (delta_lep23 < 0 && delta_lep13 < 0 && delta_lep12 > 0)
						{ // and repeat
							cut_Ztran = delta_lep12;
							caso = 3;
						}

						if ((delta_lep12 > 0 && delta_lep23 > 0) && delta_lep13 == 0 && (delta_lep12 < delta_lep23))
						{

							// case where lep 12 and 23 exist but lep 13 does not, we choose the one where 12 is smaller
							cut_Ztran = delta_lep12;
							caso = 3;
						}

						if ((delta_lep12 > 0 && delta_lep23 > 0) && delta_lep13 == 0 && (delta_lep12 > delta_lep23))
						{
							cut_Ztran = delta_lep23;
							caso = 1;
						}

						if ((delta_lep12 > 0 && delta_lep13 > 0) && delta_lep23 == 0 && (delta_lep12 < delta_lep13))
						{
							cut_Ztran = delta_lep12;
							caso = 3;
						}

						if ((delta_lep12 > 0 && delta_lep13 > 0) && delta_lep23 == 0 && (delta_lep12 > delta_lep13))
						{
							cut_Ztran = delta_lep13;
							caso = 2;
						}

						if ((delta_lep13 > 0 && delta_lep23 > 0) && delta_lep12 == 0 && (delta_lep13 < delta_lep23))
						{
							cut_Ztran = delta_lep13;
							caso = 2;
						}

						if ((delta_lep13 > 0 && delta_lep23 > 0) && delta_lep12 == 0 && (delta_lep13 > delta_lep23))
						{
							cut_Ztran = delta_lep23;
							caso = 1;
						}

						// now that the cases are done, what we need to do is build the transverse mass from each: The Lepton from the W boson and the
						// lepton pair combined, for the W it would also include the missing since there is a neutrino involved.

						float mass_leplep = 0;
						float pt_leplep = 0;
						float lep_W = 0;
						float met_W = 0;

						// sort the lepton variables, depending on the case.
						if (caso == 1)
						{
							mass_leplep = InvMass_lep23;
							pt_leplep = (Lepton_2 + Lepton_3).Pt() / 1000.;
							lep_W = Lepton_1.Pt() / 1000.;
						}
						if (caso == 2)
						{
							mass_leplep = InvMass_lep13;
							pt_leplep = (Lepton_1 + Lepton_3).Pt() / 1000.;
							lep_W = Lepton_2.Pt() / 1000.;
						}
						if (caso == 3)
						{
							mass_leplep = InvMass_lep12;
							pt_leplep = (Lepton_1 + Lepton_2).Pt() / 1000.;
							lep_W = Lepton_3.Pt() / 1000.;
						}

						if (caso == 1)
						{
							met_W = sqrt(2 * Lepton_1.Pt() * MeT.Et() * (1 - cos(Lepton_1.DeltaPhi(MeT))));
						}
						if (caso == 2)
						{
							met_W = sqrt(2 * Lepton_2.Pt() * MeT.Et() * (1 - cos(Lepton_2.DeltaPhi(MeT))));
						}
						if (caso == 3)
						{
							met_W = sqrt(2 * Lepton_3.Pt() * MeT.Et() * (1 - cos(Lepton_3.DeltaPhi(MeT))));
						}

								// now we place the cut on cut_Ztran, which is redundant but that's what we're going for. It's also supposed to be Wtran
								// but we misnamed the variable.
						if (cut_Ztran < 10.)
						{

							// now, we have made the cut for the ONE good lepton with 25GeV, we need to do the same for the other two leptons, that also have to be good. 		// we also make the cut for missing wt
							if (met_W > 30000. && (lep_pt->at(goodlep1_index) / 1000. > 25 || lep_pt->at(goodlep2_index) / 1000. > 25 || lep_pt->at(goodlep3_index) / 1000. > 25))
							{

								// at this point, all of the cuts are done, and we fill the histograms.
								// this part is confusing
								double names_of_global_variable[] = {mass_leplep, pt_leplep, met_et / 1000., met_W / 1000.};

								double names_of_leadlep_variable[] = {Lepton_1.Pt() / 1000., Lepton_1.Eta(), Lepton_1.E() / 1000., Lepton_1.Phi(), (double)lep_charge->at(goodlep1_index), (double)lep_type->at(goodlep1_index)};

								double names_of_secondlep_variable[] = {Lepton_2.Pt() / 1000., Lepton_2.Eta(), Lepton_2.E() / 1000., Lepton_2.Phi(), (double)lep_charge->at(goodlep2_index), (double)lep_type->at(goodlep2_index)};

								double names_of_thirdlep_variable[] = {Lepton_3.Pt() / 1000., Lepton_3.Eta(), Lepton_3.E() / 1000., Lepton_3.Phi(), (double)lep_charge->at(goodlep3_index), (double)lep_type->at(goodlep3_index)};

								// starting to fill the histograms
								// first the names
								TString histonames_of_global_variable[] = {"hist_mass_leplep", "hist_pt_leplep", "hist_etmiss", "hist_met_W"};

								TString histonames_of_lep_variable[] = {"hist_leptons_pt", "hist_leptons_eta", "hist_leptons_E", "hist_leptons_phi", "hist_lepttons_ch", "hist_leptons_ID"};

								// now how many

								int length_global = sizeof(names_of_global_variable) / sizeof(names_of_global_variable[0]);
								int length_leadlep = sizeof(names_of_leadlep_variable) / sizeof(names_of_leadlep_variable[0]);
								int length_secondlep = sizeof(names_of_secondlep_variable) / sizeof(names_of_secondlep_variable[0]);
								int length_thirdlep = sizeof(names_of_thirdlep_variable) / sizeof(names_of_thirdlep_variable[0]);

								// now to fill them

								for (int i = 0; i < length_global; i++)
								{
									FillHistogramsGlobal(names_of_global_variable[i], weight, histonames_of_global_variable[i]);
								}

								for (int i = 0; i < length_leadlep; i++)
								{

									FillHistogramsLeadlept(names_of_leadlep_variable[i], weight, histonames_of_lep_variable[i]);
									FillHistogramsLeadlept(names_of_secondlep_variable[i], weight, histonames_of_lep_variable[i]);
									FillHistogramsLeadlept(names_of_thirdlep_variable[i], weight, histonames_of_lep_variable[i]);
								}
							}
						}
							
						
					}
				}
			}
		}
	}

	// may have to establish two leptons

	// leptons of opposite charge

	//  std::cout<<cut1_mc<<std::endl;
	return kTRUE;
}

void WZAnalysis::SlaveTerminate()
{
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.
}

void WZAnalysis::Terminate()
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.

	name = "output_WZ/" + name + ".root";

	const char *filename = name.c_str();

	TFile physicsoutput_Top(filename, "recreate");
	WriteHistograms();
	physicsoutput_Top.Close();
}

