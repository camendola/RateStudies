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
#include "interface/object.h"
#include "interface/utils.h"

using namespace std;


int main(int argc, char** argv){


  


  //ZeroBias sample L1
  //TString directory = "/data_CMS/cms/amendola/TauTagAndProbeNtuples/VBFHTauTau_Summer17_MiniAOD_multipleTaus_04_02_18/";
  TString directory = "/data_CMS/cms/amendola/TauTagAndProbeNtuples/GluGluHTauTau_Summer17_MiniAOD_multipleTaus_04_02_18/";
  //TString directory = "/data_CMS/cms/amendola/TauTagAndProbeNtuples/DYToLL_Summer17_MiniAOD_multipleTaus_04_02_18/";



  
  TString fOutNameVBF;
  
  bool emu = false;
  
  fOutNameVBF = directory+"ggHsig_L1Plots.root";
  //  TFile *file = TFile::Open(Form("%sNtuple_VBF_L1_RECO_WithMay2017_Jets_MultipleTaus_17_05_17_0.root",directory.Data()),"read");
  TFile *file = TFile::Open(Form("%sGluGluHTauTau_Summer17_MiniAOD_multipleTaus_04_02_18.root",directory.Data()),"read");
  
  TTree * tInput = (TTree*) file->Get("Ntuplizer_noTagAndProbe_multipleTaus/TagAndProbe");
  
  
  ULong64_t       EventNumber;

  Int_t           lumi;
  std::vector<float>   *tauPt=0;
  std::vector<float>   *tauEta=0;
  std::vector<float>   *tauPhi=0;
  Int_t           JetsNumber;
  std::vector<float>   *jets_px=0;
  std::vector<float>   *jets_py=0;
  std::vector<float>   *jets_pz=0;
  std::vector<float>   *jets_e=0;
  Int_t L1_tauN;
  std::vector<Float_t>* L1_tauEt=0;
  std::vector<Float_t>* L1_tauEta=0;
  std::vector<Float_t>* L1_tauPhi=0;
  std::vector<int> *L1_tauIso=0;  
   
  Int_t L1_jetN;
  std::vector<Float_t> *L1_jetEt=0;
  std::vector<Float_t> *L1_jetEta=0;
  std::vector<Float_t> *L1_jetPhi=0;
  Int_t L1Emu_tauN;
  std::vector<Float_t>* L1Emu_tauEt=0;
  std::vector<Float_t>* L1Emu_tauEta=0;
  std::vector<Float_t>* L1Emu_tauPhi=0;
  std::vector<int> *L1Emu_tauIso=0;  
   
  Int_t L1Emu_jetN;
  std::vector<Float_t> *L1Emu_jetEt=0;
  std::vector<Float_t> *L1Emu_jetEta=0;
  std::vector<Float_t> *L1Emu_jetPhi=0;

  
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


  TBranch *b_L1_tauEt;
  TBranch *b_L1_tauEta;
  TBranch *b_L1_tauPhi;
  TBranch *b_L1_tauIso;
  
  TBranch *b_L1_jetEt;
  TBranch *b_L1_jetEta;
  TBranch *b_L1_jetPhi;

  TBranch *b_L1Emu_tauEt;
  TBranch *b_L1Emu_tauEta;
  TBranch *b_L1Emu_tauPhi;
  TBranch *b_L1Emu_tauIso;
  	     
  TBranch *b_L1Emu_jetEt;
  TBranch *b_L1Emu_jetEta;
  TBranch *b_L1Emu_jetPhi;

   
  tInput->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);

  tInput->SetBranchAddress("lumi", &lumi, &b_lumi);
  tInput->SetBranchAddress("tauPt", &tauPt, &b_tauPt);
  tInput->SetBranchAddress("tauEta", &tauEta, &b_tauEta);

  tInput->SetBranchAddress("tauPhi", &tauPhi, &b_tauPhi);
  
  tInput->SetBranchAddress("JetsNumber", &JetsNumber, &b_JetsNumber);
  tInput->SetBranchAddress("jets_px", &jets_px, &b_jets_px);
  tInput->SetBranchAddress("jets_py", &jets_py, &b_jets_py);
  tInput->SetBranchAddress("jets_pz", &jets_pz, &b_jets_pz);
  tInput->SetBranchAddress("jets_e", &jets_e, &b_jets_e);
  
  if(emu){  
  tInput ->SetBranchAddress("l1tEmuEta", &L1Emu_tauEta, &b_L1Emu_tauEta);
  tInput ->SetBranchAddress("l1tEmuPhi", &L1Emu_tauPhi, &b_L1Emu_tauPhi);
  tInput ->SetBranchAddress("l1tEmuPt", &L1Emu_tauEt, &b_L1Emu_tauEt);
  tInput ->SetBranchAddress("l1tEmuIso", &L1Emu_tauIso, &b_L1Emu_tauIso);
  tInput ->SetBranchAddress("l1tEmuEtaJet", &L1Emu_jetEta, &b_L1Emu_jetEta);
  tInput ->SetBranchAddress("l1tEmuPhiJet", &L1Emu_jetPhi, &b_L1Emu_jetPhi);
  tInput ->SetBranchAddress("l1tEmuPtJet", &L1Emu_jetEt, &b_L1Emu_jetEt);   
  }
  tInput ->SetBranchAddress("l1tEta", &L1_tauEta, &b_L1_tauEta);
  tInput ->SetBranchAddress("l1tPhi", &L1_tauPhi, &b_L1_tauPhi);
  tInput ->SetBranchAddress("l1tPt", &L1_tauEt, &b_L1_tauEt);
  tInput ->SetBranchAddress("l1tIso", &L1_tauIso, &b_L1_tauIso);
  tInput ->SetBranchAddress("l1tEtaJet", &L1_jetEta, &b_L1_jetEta);
  tInput ->SetBranchAddress("l1tPhiJet", &L1_jetPhi, &b_L1_jetPhi);
  tInput ->SetBranchAddress("l1tPtJet", &L1_jetEt, &b_L1_jetEt);   
  
  
  
  
  
  ///////////////
  //// HISTO ////
  ///////////////
  
  TFile* fOutVBF = new TFile (Form("%s",fOutNameVBF.Data()), "recreate");
  
 
  TH1D* SubJetMjj620_OffMatchUn = new TH1D ("SubJetMjj620_OffMatchUn", "", 100,25,125);
  TH1D* SubJetMjj620_OffMatchEmu = new TH1D ("SubJetMjj620_OffMatchEmu", "", 100,25,125);
  TH1D* SubJetMjj620_Emu = new TH1D ("SubJetMjj620_Emu", "", 100,20,120);
  TH1D* SubJetMjj620_Un = new TH1D ("SubJetMjj620_Un", "", 100,20,120);


  TH1D* DeltaEta = new TH1D ("DeltaEta", "", 50,0,10);
  TH1D* Eta = new TH1D ("Eta", "", 50,-5,5);

  TH1D* Eta_OffMatchEmu = new TH1D ("Eta_OffMatchEmu", "", 50,-5,5);
  TH1D* Eta_OffMatchUn = new TH1D ("Eta_OffMatchUn", "", 50,-5,5);
  TH1D* DeltaEta_OffMatchUn = new TH1D ("DeltaEta_OffMatchUn", "", 50,0,10);
  TH1D* Eta_Mjj620_OffMatchEmu = new TH1D ("Eta_Mjj620_OffMatchEmu", "", 50,-5,5);
  TH1D* Eta_Mjj620_OffMatchUn = new TH1D ("Eta_Mjj620_OffMatchUn", "", 50,-5,5);

  
  TH1D* Et_OffMatchEmu = new TH1D ("Et_OffMatchEmu", "", 100,25,125);
  TH1D* Et_OffMatchUn = new TH1D ("Et_OffMatchUn", "", 100,25,125);

  TH1D* Mjj_Off = new TH1D ("Mjj_Off", "", 100,0,800);
  TH1D* LeadJet_Off = new TH1D ("LeadJet_Off", "", 120,30,150);
  TH1D* SubLeadJet_Off = new TH1D ("SubLeadJet_Off", "", 120,30,150);
  TH1D* Mjj_OffMatchUn = new TH1D ("Mjj_OffMatchUn", "", 100,0,800);

  TH1D* Off_passVBF = new TH1D ("Off_passVBF", "", 100,30,130);
  TH1D* On_passVBF = new TH1D ("On_passVBF", "", 100,30,130);
  TH1D* Int_Off_passVBF = new TH1D ("Int_Off_passVBF", "", 100,30,130);
  TH1D* Int_On_passVBF = new TH1D ("Int_On_passVBF", "", 100,30,130);

  TH2D* Off_passVBF_XY = new TH2D ("Off_passVBF_XY", "", 30,35,65,100,30,130);
  TH2D* On_passVBF_XY = new TH2D ("On_passVBF_XY", "", 30,35,65,100,30,130);
  TH2D* Int_Off_passVBF_XY = new TH2D ("Int_Off_passVBF_XY", "", 30,35,65,100,30,130);
  TH2D* Int_On_passVBF_XY = new TH2D ("Int_On_passVBF_XY", "", 30,35,65,100,30,130);


  TH2D* Off_passVBF_XZ = new TH2D ("Off_passVBF_XZ", "", 800,0,800,100,30,130);
  TH2D* On_passVBF_XZ = new TH2D ("On_passVBF_XZ", "", 800,0,800,100,30,130);
  TH2D* Int_Off_passVBF_XZ = new TH2D ("Int_Off_passVBF_XZ", "", 800,0,800,100,30,130);
  TH2D* Int_On_passVBF_XZ = new TH2D ("Int_On_passVBF_XZ", "", 800,0,800,100,30,130);

  TGraphAsymmErrors* VBF1Defficiency = new  TGraphAsymmErrors();
  VBF1Defficiency->SetName("VBF1Defficiency");



  
  


  
  
  // analyze data    
  
  std::vector<object> jet20un;
  std::vector<object> jet20emu;    
  std::vector<TLorentzVector> jetOff;
  std::vector<TLorentzVector> tauOff;
  std::vector<TLorentzVector> jetOff_matchEmu;
  std::vector<TLorentzVector> jetOff_matchUn;
  std::vector< tuple<double, int,int> > mjj_un;
  std::vector< tuple<double, int,int> > mjj_emu;
  std::vector< tuple<double, int,int> > mjj30_emu;
  std::vector< tuple<double, int,int,double,double> > mjj30_emu_sortbyPt;
  std::vector< tuple<double, int,int> > mjj40_off;
  std::vector< tuple<double, int,int,double,double> > mjj40_off_sortbyPt;
  std::vector< tuple<double, int,int> > mjj_off_matchEmu;
  std::vector< tuple<double, int,int> > mjj_off_matchUn;
  
  int Nentries =     tInput->GetEntries();
  for (Long64_t iEv = 0 ;iEv<Nentries ; ++iEv){
  
   
    tInput->GetEntry(iEv);
   

     
    if (iEv%1000 == 0) cout << iEv << endl;
     
          
    jet20un.clear();
    jet20emu.clear();
    jetOff.clear();
    tauOff.clear();

    jetOff_matchEmu.clear();
    jetOff_matchUn.clear();
    mjj_un.clear();
    mjj_emu.clear();
    mjj30_emu.clear();
    mjj40_off.clear();
    mjj30_emu_sortbyPt.clear();
    mjj40_off_sortbyPt.clear();
    mjj_off_matchEmu.clear();
    mjj_off_matchUn.clear();
    
    for (long int iL1 = 0; iL1 < L1_jetEt->size(); iL1++){ //loop on jets unpacked
      // selections
      double jetPtUn  = L1_jetEt->at(iL1);
      if(jetPtUn>20.) jet20un.push_back(object(L1_jetEt->at(iL1),L1_jetEta->at(iL1),L1_jetPhi->at(iL1),-999)) ;
    }
    if (emu){
      for (long int iL1emu = 0; iL1emu < L1Emu_jetEt->size(); iL1emu++){ //loop on jets emulated
	double jetPtEmu  = L1Emu_jetEt->at(iL1emu);
	if(jetPtEmu>20.) jet20emu.push_back(object(L1Emu_jetEt->at(iL1emu),L1Emu_jetEta->at(iL1emu),L1Emu_jetPhi->at(iL1emu),-999)) ;
      }
    }

    for(int itau = 0; itau<tauPt->size(); itau++ ){
	
      TLorentzVector tau;
      tau.SetPtEtaPhiM(tauPt->at(itau), tauPhi->at(itau), tauEta->at(itau),0.);
      if(tau.Pt()>25.)  tauOff.push_back(tau);
    }
        
    std::sort(tauOff.begin(),tauOff.end(),SortByPt);
    
    for (long int ijet = 0; ijet < jets_px->size(); ijet++){ //loop on jets offline
      TLorentzVector jet;
      jet.SetPxPyPzE(jets_px->at(ijet),jets_py->at(ijet),jets_pz->at(ijet),jets_e->at(ijet));
      double jetPt  = jet.Pt();
      if(jetPt>25.){
	int nTaus = tauOff.size();
	if (nTaus > 2) nTaus = 2;
	for(int itau = 0; itau < nTaus; itau++){
	  if (jet.DeltaR(tauOff[itau])>0.2){
	    jetOff.push_back(jet);
	    object::SortByDeltaR DeltaRJetOff(jetOff[ijet].Phi(),jetOff[ijet].Eta());
	    if(jet20emu.size()>0){
	      std::sort(jet20emu.begin(),jet20emu.end(),DeltaRJetOff);
	      if(jet20emu[0].DeltaRTLV(jetOff[ijet])<0.2) jetOff_matchEmu.push_back(jet) ;
	    }
	    std::sort(jet20un.begin(),jet20un.end(),DeltaRJetOff);
	    if(jet20un[0].DeltaRTLV(jetOff[ijet])<0.2) jetOff_matchUn.push_back(jet) ;
	  }
	}
      }
    }

    std::sort(jet20un.begin(),jet20un.end());

    std::sort(jet20emu.begin(),jet20emu.end());

    std::sort(jetOff.begin(),jetOff.end(),SortByPt);
    std::sort(jetOff_matchEmu.begin(),jetOff_matchEmu.begin(),SortByPt);
    std::sort(jetOff_matchUn.begin(),jetOff_matchUn.begin(),SortByPt);

    //unpacked  
    if (jet20un.size() >= 2){
      for (int iJet = 0; iJet <jet20un.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet20un.size(); kJet++){      
	  if (kJet!=iJet) {
	    
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jet20un[iJet].Et(),
			      jet20un[iJet].Eta(),
			      jet20un[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jet20un[kJet].Et(),
			      jet20un[kJet].Eta(),
			      jet20un[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_un.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    if(jetPair.M()>620) SubJetMjj620_Un->Fill(kjet.Pt());
	  }
	  
	}
	
      }
      
      std::sort(mjj_un.begin(),mjj_un.end());
      
    }
    
    //emulated
    if (jet20emu.size() >= 2){
      for (int iJet = 0; iJet <jet20emu.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet20emu.size(); kJet++){      
	  if (kJet!=iJet) {
	    
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jet20emu[iJet].Et(),
			      jet20emu[iJet].Eta(),
			      jet20emu[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jet20emu[kJet].Et(),
			      jet20emu[kJet].Eta(),
			      jet20emu[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;

	    mjj_emu.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    if (jet20emu[kJet].Et()>35) {
	      mjj30_emu.push_back(make_tuple(jetPair.M(),iJet,kJet));
	      mjj30_emu_sortbyPt.push_back(make_tuple(jetPair.M(),iJet,kJet,jet20emu[iJet].Et(),jet20emu[kJet].Et()));
	    }
	    if(jetPair.M()>620) SubJetMjj620_Emu->Fill(kjet.Pt());
	  }
	}
      }
      std::sort(mjj_emu.begin(),mjj_emu.end());
      std::sort(mjj30_emu.begin(),mjj30_emu.end());
      std::sort(mjj30_emu_sortbyPt.begin(),mjj30_emu_sortbyPt.end(),SortMjjByJetThreshold);
    }
  
    //offline
    
    if(jetOff_matchUn.size()>=2){
      for (int iJet = 0; iJet <jetOff_matchUn.size(); iJet++){
	    
	for (int kJet = iJet+1; kJet <jetOff_matchUn.size(); kJet++){
	  if (kJet!=iJet) {
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jetOff_matchUn[iJet].Pt(),
			      jetOff_matchUn[iJet].Eta(),
			      jetOff_matchUn[iJet].Phi(),
			      jetOff_matchUn[iJet].M()
			      );
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jetOff_matchUn[kJet].Pt(),
			      jetOff_matchUn[kJet].Eta(),
			      jetOff_matchUn[kJet].Phi(),
			      jetOff_matchUn[kJet].M()
			      );
	    TLorentzVector jetPair = ijet+kjet;
	    if (jetOff_matchUn[kJet].Pt()>30){
	    mjj_off_matchUn.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    }
	    if(jetPair.M()>700) {
	      SubJetMjj620_OffMatchUn->Fill(kjet.Pt());
	      Eta_Mjj620_OffMatchUn->Fill(kjet.Eta());    
	    }

	  }
	}
	std::sort(mjj_off_matchUn.begin(),mjj_off_matchUn.end());

	if(mjj_off_matchUn.size()>0){
	  if (std::get<0>(*(mjj_off_matchUn.rbegin()))> 80){
	    Mjj_OffMatchUn->Fill(std::get<0>(*(mjj_off_matchUn.rbegin())));
	    DeltaEta_OffMatchUn->Fill(fabs(jetOff[std::get<1>(*(mjj_off_matchUn.rbegin()))].Eta()-jetOff[std::get<2>(*(mjj_off_matchUn.rbegin()))].Eta()));
	    
	  }
	  
	}
      }
      if(jetOff_matchEmu.size()>=2){
      for (int iJet = 0; iJet <jetOff_matchEmu.size(); iJet++){
	    
	for (int kJet = iJet+1; kJet <jetOff_matchEmu.size(); kJet++){
	  if (kJet!=iJet) {
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jetOff_matchEmu[iJet].Pt(),
			      jetOff_matchEmu[iJet].Eta(),
			      jetOff_matchEmu[iJet].Phi(),
			      jetOff_matchEmu[iJet].M()
			      );
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jetOff_matchEmu[kJet].Pt(),
			      jetOff_matchEmu[kJet].Eta(),
			      jetOff_matchEmu[kJet].Phi(),
			      jetOff_matchEmu[kJet].M()
			      );
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_off_matchEmu.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    if(jetPair.M()>700) {
	      SubJetMjj620_OffMatchEmu->Fill(kjet.Pt());
	      Eta_Mjj620_OffMatchEmu->Fill(kjet.Eta());    
	    }
	  }

	}
      }
      std::sort(mjj_off_matchEmu.begin(),mjj_off_matchEmu.end());
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
	    if (jetOff[kJet].Pt()>30){
	      mjj40_off.push_back(make_tuple(jetPair.M(),iJet,kJet));
	      mjj40_off_sortbyPt.push_back(make_tuple(jetPair.M(),iJet,kJet,jetOff[iJet].Pt(),jetOff[kJet].Pt()));
	      
	    }
	  }
	}
	}
    
      std::sort(mjj40_off.begin(),mjj40_off.end());
      std::sort(mjj40_off_sortbyPt.begin(),mjj40_off_sortbyPt.end(),SortMjjByJetThreshold);
      if(mjj40_off.size()>0){
	if (std::get<0>(*(mjj40_off.rbegin()))> 80){
	Mjj_Off->Fill(std::get<0>(*(mjj40_off.rbegin())));
	
	DeltaEta->Fill(fabs(jetOff[std::get<1>(*(mjj40_off.rbegin()))].Eta()-jetOff[std::get<2>(*(mjj40_off.rbegin()))].Eta()));
	Eta->Fill(jetOff[std::get<1>(*(mjj40_off.rbegin()))].Eta());
	}
      }
      LeadJet_Off->Fill(jetOff[0].Pt());
      SubLeadJet_Off->Fill(jetOff[1].Pt());
    }
    }
    for(int ijet = 0; ijet<jetOff_matchEmu.size();ijet++){
      Eta_OffMatchEmu->Fill(jetOff_matchEmu[ijet].Eta());
      Et_OffMatchEmu->Fill(jetOff_matchEmu[ijet].Et());    
    }
    
    for(int ijet = 0; ijet<jetOff_matchUn.size();ijet++){
      Eta_OffMatchUn->Fill(jetOff_matchUn[ijet].Eta());
      Et_OffMatchUn->Fill(jetOff_matchUn[ijet].Et());    
    }

    //plots1D
    if(mjj40_off.size()>0 && jetOff.size()>0){
      if(std::get<0>(*(mjj40_off.rbegin()))>720 && std::get<4>(*(mjj40_off_sortbyPt.rbegin()))>45) Off_passVBF->Fill(jetOff[0].Pt());
      if(mjj30_emu.size()>0 && jet20emu.size()>0){
	if(jet20emu[0].Et()>jetOff[0].Pt()-10){
	  if( std::get<0>(*(mjj30_emu.rbegin()))>620 && std::get<0>(*(mjj40_off.rbegin()))>720 && std::get<4>(*(mjj40_off_sortbyPt.rbegin()))>45 && std::get<4>(*(mjj30_emu_sortbyPt.rbegin()))>35 )  On_passVBF->Fill(jetOff[0].Pt());
	    
	}
      }
    }
    //plots2D
    if(mjj40_off.size()>0 && jetOff.size()>0){
      if(std::get<0>(*(mjj40_off.rbegin()))>720) Off_passVBF_XY->Fill(std::get<4>(*(mjj40_off_sortbyPt.rbegin())),jetOff[0].Pt());
      if (std::get<4>(*(mjj40_off_sortbyPt.rbegin()))>45)Off_passVBF_XZ->Fill(std::get<0>(*(mjj40_off.rbegin())),jetOff[0].Pt());
      if(mjj30_emu.size()>0 && jet20emu.size()>0){
	if((jet20emu[0].Et()>jetOff[0].Pt()-10) &&(std::get<4>(*(mjj30_emu_sortbyPt.rbegin()))>std::get<4>(*(mjj40_off_sortbyPt.rbegin()))-10)&&(std::get<0>(*(mjj30_emu.rbegin()))>std::get<0>(*(mjj40_off.rbegin()))-100)){
	  if( std::get<0>(*(mjj30_emu.rbegin()))>620 && std::get<0>(*(mjj40_off.rbegin()))>720)  On_passVBF_XY->Fill(std::get<4>(*(mjj40_off_sortbyPt.rbegin())),jetOff[0].Pt());
	  if(std::get<4>(*(mjj40_off_sortbyPt.rbegin()))>45 && std::get<4>(*(mjj30_emu_sortbyPt.rbegin()))>35 ) On_passVBF_XZ->Fill(std::get<0>(*(mjj40_off.rbegin())),jetOff[0].Pt());
	}
      }
    }

      
      
  }

  for (int i = 1; i<= On_passVBF->GetNbinsX(); i++){
    Int_On_passVBF->SetBinContent(i,On_passVBF->Integral(i,On_passVBF->GetNbinsX()+1));
    Int_Off_passVBF->SetBinContent(i,Off_passVBF->Integral(i,Off_passVBF->GetNbinsX()+1));
  }

  for (int i = 1; i<= On_passVBF_XY->GetNbinsX(); i++){
    for (int j = 1; j<= On_passVBF_XY->GetNbinsY(); j++){
      Int_On_passVBF_XY->SetBinContent(i,j,On_passVBF_XY->Integral(i,On_passVBF_XY->GetNbinsX()+1,j,On_passVBF_XY->GetNbinsY()+1));
      Int_Off_passVBF_XY->SetBinContent(i,j,Off_passVBF_XY->Integral(i,Off_passVBF_XY->GetNbinsX()+1,j,Off_passVBF_XY->GetNbinsY()+1));
    }
  }
  for (int i = 1; i<= On_passVBF_XZ->GetNbinsX(); i++){
    for (int j = 1; j<= On_passVBF_XZ->GetNbinsY(); j++){
      Int_On_passVBF_XZ->SetBinContent(i,j,On_passVBF_XZ->Integral(i,On_passVBF_XZ->GetNbinsX()+1,j,On_passVBF_XZ->GetNbinsY()+1));
      Int_Off_passVBF_XZ->SetBinContent(i,j,Off_passVBF_XZ->Integral(i,Off_passVBF_XZ->GetNbinsX()+1,j,Off_passVBF_XZ->GetNbinsY()+1));
    }
  }
  
  VBF1Defficiency->BayesDivide(Int_On_passVBF,Int_Off_passVBF);
  
  TH2D* VBF2Defficiency_XY = (TH2D*)Int_On_passVBF_XY->Clone("VBF2Defficiency_XY");
  VBF2Defficiency_XY->Divide(Int_Off_passVBF_XY);
  TH2D* VBF2Defficiency_XZ = (TH2D*)Int_On_passVBF_XZ->Clone("VBF2Defficiency_XZ");
  VBF2Defficiency_XZ->Divide(Int_Off_passVBF_XZ);

    //TEfficiency* VBF2Defficiency_XZ = new  TEfficiency(Int_On_passVBF_XZ,Int_Off_passVBF_XZ);
    //VBF2Defficiency_XZ->SetName("VBF2Defficiency_XZ");


  fOutVBF->cd();
  VBF1Defficiency->Write();
  //VBF2Defficiency_XY->Write();
  //VBF2Defficiency_XZ->Write();
  fOutVBF->Write();
  cout<<"Output saved in "<<directory<<endl;
  }
