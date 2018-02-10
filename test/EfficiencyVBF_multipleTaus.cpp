#include <sstream>
#include <fstream>
#include <map>
#include <iostream>
#include <fstream>
#include <utility>
#include <tuple>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "../interface/object.h"
#include "../interface/utils.h"

using namespace std;






int main(int argc, char** argv){


  double Mjj_cut_on = 620; 
  double lead_cut_on = 115;
  double sub_cut_on = 40;


  double Mjj_cut_off = 900; 
  double lead_cut_off = 150;
  double sub_cut_off = 60;
  //TString fileList = "fileLists/SingleMuon.txt";
  //TString fileList = "fileLists/SignalMC/TagAndProbe_forData_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer17MiniAOD-NZSFlatPU28to62_SUS01_92X_upgrade2017_realistic_v10-v1.txt";
  TString fileList = "fileLists/TagAndProbe_forData_W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v1.txt";
 
  //ZeroBias sample L1
  TString directory = "/data_CMS/cms/amendola/TauTagAndProbeNtuples/W2JetsToLNu_Summer17_MiniAOD_multipleTaus_HLTfilt_08_02_18/";
   // TString directory = "/data_CMS/cms/amendola/TauTagAndProbeNtuples/WJetsToLNu_Summer17_MiniAOD_multipleTaus_HLTfilt_08_02_18/";
   



  
  TString fOutNameVBF;
  
  
  
  //  fOutNameVBF = "VBFEff_run2017H.root";
  // fOutNameVBF = directory+"WJetsToNuEff_2017.root";
   fOutNameVBF = directory+"W2JetsEff_2017.root";

  TChain * tInput;

  tInput= new TChain ("Ntuplizer_noTagAndProbe_multipleTaus/TagAndProbe");
  appendFromFileList(tInput, fileList);





  ULong64_t       EventNumber;

  Int_t           lumi;
  std::vector<float>  * tauPt=0;
  std::vector<float>  * tauEta=0;
  std::vector<float>  * tauPhi=0;
  Int_t           JetsNumber;
  std::vector<float>   *jets_px=0;
  std::vector<float>   *jets_py=0;
  std::vector<float>   *jets_pz=0;
  std::vector<float>   *jets_e=0;
  std::vector<float>   *jets_leptonDeltaR= 0;
  std::vector<float>   *jets_PUJetID = 0;
  Int_t L1_tauN;
  std::vector<Float_t>* L1_tauEt=0;
  std::vector<Float_t>* L1_tauEta=0;
  std::vector<Float_t>* L1_tauPhi=0;
  std::vector<int> *L1_tauIso=0;  
   
  Int_t L1_jetN;
  std::vector<Float_t> *L1_jetEt=0;
  std::vector<Float_t> *L1_jetEta=0;
  std::vector<Float_t> *L1_jetPhi=0;


  
  TBranch        *b_EventNumber;   //!

  TBranch        *b_lumi;   //!
  TBranch        *b_tauPt;   //!
  TBranch        *b_tauPhi;   //!
  TBranch        *b_tauEta;   //!

  TBranch        *b_JetsNumber;   //!
  TBranch        *b_jets_px;   //!
  TBranch        *b_jets_py;   //!
  TBranch        *b_jets_pz;   //!
  TBranch        *b_jets_e;   //!
  TBranch        *b_jets_leptonDeltaR;
  TBranch        *b_jets_PUJetID;


  TBranch *b_L1_tauEt;
  TBranch *b_L1_tauEta;
  TBranch *b_L1_tauPhi;
  TBranch *b_L1_tauIso;
  
  TBranch *b_L1_jetEt;
  TBranch *b_L1_jetEta;
  TBranch *b_L1_jetPhi;
  


  tInput->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);

  tInput->SetBranchAddress("lumi", &lumi, &b_lumi);
  tInput->SetBranchAddress("tauPt", &tauPt, &b_tauPt);
  tInput->SetBranchAddress("tauEta", &tauEta, &b_tauEta);
  tInput->SetBranchAddress("tauPhi", &tauPhi, &b_tauPhi);


  
  tInput->SetBranchAddress("l1tPtJet", &L1_jetEt,  &b_L1_jetEt);
  tInput->SetBranchAddress("l1tEtaJet", &L1_jetEta, &b_L1_jetEta);
  tInput->SetBranchAddress("l1tPhiJet", &L1_jetPhi, &b_L1_jetPhi);

  tInput->SetBranchAddress("JetsNumber", &JetsNumber, &b_JetsNumber);
  tInput->SetBranchAddress("jets_px", &jets_px, &b_jets_px);
  tInput->SetBranchAddress("jets_py", &jets_py, &b_jets_py);
  tInput->SetBranchAddress("jets_pz", &jets_pz, &b_jets_pz);
  tInput->SetBranchAddress("jets_e", &jets_e, &b_jets_e);
  tInput->SetBranchAddress("jets_leptonDeltaR",&jets_leptonDeltaR,&b_jets_leptonDeltaR);
  tInput->SetBranchAddress("jets_PUJetID", &jets_PUJetID,&b_jets_PUJetID);


  
  TFile* fOutVBF = new TFile (Form("%s",fOutNameVBF.Data()), "recreate");
  
 
  

  TH1D* Off_passVBF_pt = new TH1D ("Off_passVBF_pt", "", 100,50,300);
  TH1D* On_passVBF_pt = new TH1D ("On_passVBF_pt", "", 100,50,300);

  TH1D* Off_passVBF_mjj = new TH1D ("Off_passVBF_mjj", "", 100,500,1500);
  TH1D* On_passVBF_mjj = new TH1D ("On_passVBF_mjj", "", 100,500,1500);

  TH1D* Mjj_Off = new TH1D ("Mjj_Off", "", 100,500,1500);
  TH1D* Mjj_On = new TH1D ("Mjj_On", "", 100,500,1500);
  
  
  
  // analyze data    
  
  std::vector<object> jetL1;

  std::vector<TLorentzVector> jetOff;
  std::vector<TLorentzVector> tauOff;

  std::vector<TLorentzVector> jetOffmatch;

 

  std::vector< tuple<double, int,int> > mjj_on;
  
  std::vector< tuple<double, int,int> > mjj_off;
  std::vector< tuple<double, int,int> > mjj_off_match;

  



  for (Long64_t iEv = 0 ;true ; ++iEv){
  
    int got = 0;

    got = tInput->GetEntry(iEv);


    
    if (got == 0) break;
   


    if (iEv%1000 == 0) cout << iEv << endl;

    
    jetL1.clear();

    jetOff.clear();
    jetOffmatch.clear();
    tauOff.clear();

    mjj_on.clear();
    mjj_off.clear();
    mjj_off_match.clear();


    for (long int iL1 = 0; iL1 < L1_jetEt->size(); iL1++){ //loop on jets emulated
      // selections
      double jetPt  = L1_jetEt->at(iL1);
      if(jetPt>sub_cut_on) jetL1.push_back(object(L1_jetEt->at(iL1),L1_jetEta->at(iL1),L1_jetPhi->at(iL1),-999)) ;
    }
    std::sort(jetL1.begin(),jetL1.end());





    for (long int itau = 0; itau < tauPt->size(); itau++){ //loop on jets offline
      TLorentzVector tau;

      tau.SetPtEtaPhiM(tauPt->at(itau), tauPhi->at(itau), tauEta->at(itau),0.);
      tauOff.push_back(tau);
    }
    
    std::sort(tauOff.begin(),tauOff.end(),SortByPt);
    
    for (long int ijet = 0; ijet < jets_px->size(); ijet++){ //loop on jets offline
      TLorentzVector jet;
      if(jets_PUJetID->at(ijet)<0.) continue;
      if(fabs(jets_leptonDeltaR->at(ijet))<0.5) continue;
      jet.SetPxPyPzE(jets_px->at(ijet),jets_py->at(ijet),jets_pz->at(ijet),jets_e->at(ijet));
      double jetPt  = jet.Pt();
      if(jetPt>sub_cut_off){
       	if(tauOff.size()>0){
	  if (tauOff[0].Pt()>sub_cut_off){
	    if(jet.DeltaR(tauOff[0])>0.2)   jetOff.push_back(jet);
	    
	  }else{
	    jetOff.push_back(jet);
	  }
	}else{
	  jetOff.push_back(jet);
	}
	for (long int iL1 = 0; iL1 < L1_jetEt->size(); iL1++){ //loop on jets emulated
	  // selections
	  double jetPt  = L1_jetEt->at(iL1);
	  if(jetPt>20) {
	    TLorentzVector L1jet;
	    L1jet.SetPtEtaPhiM(L1_jetEt->at(iL1),L1_jetEta->at(iL1),L1_jetPhi->at(iL1),0);
	    if (L1jet.DeltaR(jet)<0.5){
	      if(tauOff.size()>0){
		if (tauOff[0].Pt()>sub_cut_off){
		  if(jet.DeltaR(tauOff[0])>0.2)   jetOffmatch.push_back(jet);
		  
		}else{
		  jetOffmatch.push_back(jet);
		}
	      }else{
		jetOffmatch.push_back(jet);
	      }
	    }
	  }
	}
      }
    }    
    
    
    
    std::sort(jetOff.begin(),jetOff.end(),SortByPt);
    
    
    //unpacked  
    if (jetL1.size() >= 2){
      for (int iJet = 0; iJet <jetL1.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jetL1.size(); kJet++){      
	  if (kJet!=iJet) {
	    
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jetL1[iJet].Et(),
			      jetL1[iJet].Eta(),
			      jetL1[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jetL1[kJet].Et(),
			      jetL1[kJet].Eta(),
			      jetL1[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_on.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    
	  }
	  
	}
	
      }
      
      std::sort(mjj_on.begin(),mjj_on.end());
      
    }
    
    
    
    if(jetOff.size()>=2){
      for (int iJet = 0; iJet <jetOff.size(); iJet++){
	    
	for (int kJet = iJet+1; kJet <jetOff.size(); kJet++){
	  if (kJet!=iJet) {
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jetOff[iJet].Pt(),
			      jetOff[iJet].Eta(),
			      jetOff[iJet].Phi(),
			      jetOff[iJet].M()
			      );
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jetOff[kJet].Pt(),
			      jetOff[kJet].Eta(),
			      jetOff[kJet].Phi(),
			      jetOff[kJet].M()
			      );
	    TLorentzVector jetPair = ijet+kjet;
	      mjj_off.push_back(make_tuple(jetPair.M(),iJet,kJet));

	  }

	}
      }
        
      std::sort(mjj_off.begin(),mjj_off.end());
      double Mjj_off = std::get<0>(*(mjj_off.rbegin()));
      double lead_off = jetOff[0].Pt();

      double Mjj_on = 0;
      if (mjj_on.size()>0)Mjj_on= std::get<0>(*(mjj_on.rbegin()));

      double lead_on = 0;
      if (jetL1.size()>1) lead_on = jetL1[0].Et();
      bool passOff_mjj = false;
      bool passOff_lead = false;
      if(Mjj_off>Mjj_cut_off) {
	passOff_mjj = true;
	Off_passVBF_pt->Fill(lead_off);
      }
      if(lead_off>lead_cut_off) {
	passOff_lead = true;
	Off_passVBF_mjj->Fill(Mjj_off);

      }
      //      if (passOff_mjj == true && Mjj_on>Mjj_cut_on && lead_on>lead_cut_on) On_passVBF_pt->Fill(lead_off);
      // if (passOff_lead == true && Mjj_on>Mjj_cut_on && lead_on>lead_cut_on) On_passVBF_mjj->Fill(Mjj_off);
      


      
      
    }

    if(jetOffmatch.size()>=2){
      for (int iJet = 0; iJet <jetOffmatch.size(); iJet++){
	
	for (int kJet = iJet+1; kJet <jetOffmatch.size(); kJet++){
	  if (kJet!=iJet) {
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jetOffmatch[iJet].Pt(),
			      jetOffmatch[iJet].Eta(),
			      jetOffmatch[iJet].Phi(),
			      jetOffmatch[iJet].M()
			      );
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jetOffmatch[kJet].Pt(),
			      jetOffmatch[kJet].Eta(),
			      jetOffmatch[kJet].Phi(),
			      jetOffmatch[kJet].M()
			      );
	    TLorentzVector jetPair = ijet+kjet;
	      mjj_off_match.push_back(make_tuple(jetPair.M(),iJet,kJet));

	  }

	}
      }
        
      std::sort(mjj_off_match.begin(),mjj_off_match.end());
      double Mjj_off = std::get<0>(*(mjj_off_match.rbegin()));
      double lead_off = jetOffmatch[0].Pt();

      double Mjj_on = 0;
      if (mjj_on.size()>0)Mjj_on= std::get<0>(*(mjj_on.rbegin()));

      double lead_on = 0;
      if (jetL1.size()>1) lead_on = jetL1[0].Et();
      bool passOff_mjj = false;
      bool passOff_lead = false;
      if(Mjj_off>Mjj_cut_off) {
	passOff_mjj = true;
	//Off_passVBF_pt->Fill(lead_off);
      }
      if(lead_off>lead_cut_off) {
	passOff_lead = true;
	//Off_passVBF_mjj->Fill(Mjj_off);
	Mjj_Off->Fill(Mjj_off);
      }
      if (passOff_mjj == true && Mjj_on>Mjj_cut_on && lead_on>lead_cut_on) On_passVBF_pt->Fill(lead_off);
      if (passOff_lead == true && lead_on>lead_cut_on) {
	if(Mjj_on>Mjj_cut_on) On_passVBF_mjj->Fill(Mjj_off);
	Mjj_On->Fill(Mjj_on);
      }
    }
  }


  fOutVBF->cd();

  fOutVBF->Write();
  cout<<fOutNameVBF<<endl;
}
