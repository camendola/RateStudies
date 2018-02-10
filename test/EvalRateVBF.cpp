// evaluates rate  - compile with c++ -lm -o EvalRateVBF test/EvalRateVBF.cpp `root-config --glibs --cflags`
// and launch the executable 
#include <sstream>
#include <fstream>
#include <map>
#include <iostream>
#include <fstream>
#include <utility>
#include <regex>
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
#include "../interface/object.h"

using namespace std;

void appendFromFileList (TChain* chain, TString filename)
{
 
  std::ifstream infile(filename.Data());
  std::string line;
  while (std::getline(infile, line))
    {
      line = line.substr(0, line.find("#", 0)); // remove comments introduced by #
      while (line.find(" ") != std::string::npos) line = line.erase(line.find(" "), 1); // remove white spaces
      while (line.find("\n") != std::string::npos) line = line.erase(line.find("\n"), 1); // remove new line characters
      while (line.find("\r") != std::string::npos) line = line.erase(line.find("\r"), 1); // remove carriage return characters
      if (!line.empty()) // skip empty lines
	chain->Add(line.c_str());
    }
  return;
}


int main(int argc, char** argv){

  int Mjj_cut = atoi(argv[1]); //for 1-d plots
  int leadjet_cut = atoi(argv[2]); //for 1-d plots
  int subleadjet_cut = atoi(argv[3]); //for 1-d plots
  TString directory = argv[4];
  TString fileList  = argv[5];
  TString PUperLumiFile = argv[6];
  ifstream PUFile(PUperLumiFile);
  int trigger_ID = atoi(argv[7]);
  bool emulated  =false;
  if(atoi(argv[8])==1) emulated = true; 
  int Run = atoi(argv[9]);
  int Fill = atoi(argv[10]);
  int maxevents = atoi(argv[11]);
  int lumimin = atoi(argv[12]);
  int lumimax = atoi(argv[13]);
  std::stringstream max_ss;
  max_ss << maxevents;
  TString max_str = max_ss.str();

  int njob = atoi(argv[14]);

  bool reweight = false;
  
  


  
  TTree * tInput;


  TChain * cInput;
  TChain * lInput;
  TChain * ugtInput; 



  
  if (emulated){
    cInput = new TChain ("l1UpgradeEmuTree/L1UpgradeTree");
    appendFromFileList(cInput, fileList);
    lInput = new TChain ("l1EventTree/L1EventTree");
    appendFromFileList(lInput, fileList);
  }else{
    cInput= new TChain ("l1UpgradeTree/L1UpgradeTree");
    appendFromFileList(cInput, fileList);
    lInput = new TChain ("l1EventTree/L1EventTree");
    appendFromFileList(lInput, fileList);

  }

      ugtInput = new TChain ("l1uGTTree/L1uGTTree");
      appendFromFileList(ugtInput, fileList);


    
    
  TString fOutNameVBF = "VBF_rate" ;
  if(emulated) fOutNameVBF += "_emu";
  fOutNameVBF += Form("_%d",njob);
  fOutNameVBF += ".root";
  
  
  lInput->SetMakeClass(1);
  cInput->SetMakeClass(1);  
  ugtInput->SetMakeClass(1);
  
  Int_t lumi ;
  UInt_t run ;
  
  UShort_t stage2_jetN ;
  std::vector<float> stage2_jetEt;
  std::vector<float> stage2_jetEta;
  std::vector<float> stage2_jetPhi;
  std::vector<short> stage2_jetBx;
  
  std::vector<float> stage2_tauEt;
  std::vector<float> stage2_tauEta;
  std::vector<float> stage2_tauPhi;
  std::vector<short> stage2_tauIso;
  std::vector<short> stage2_tauBx;
  
  
  std::vector<bool>  trigger_algo;
  
  // set branch and variables
  TBranch *b_lumi ;
  TBranch *b_run ;

  
  TBranch *b_stage2_jetN ;
  TBranch *b_stage2_jetEt;
  TBranch *b_stage2_jetEta;
  TBranch *b_stage2_jetPhi;
  TBranch *b_stage2_jetBx;
  
  TBranch *b_stage2_tauEt;
  TBranch *b_stage2_tauEta;
  TBranch *b_stage2_tauPhi;
  TBranch *b_stage2_tauIso;
  TBranch *b_stage2_tauBx;
  
  TBranch *b_algoDecisionInitial;
  
  
  
  lInput ->SetBranchAddress("lumi", &lumi , &b_lumi );
  lInput ->SetBranchAddress("run", &run , &b_run );

  cInput ->SetBranchAddress("nJets", &stage2_jetN , &b_stage2_jetN);
  cInput ->SetBranchAddress("jetEta", &stage2_jetEta, &b_stage2_jetEta);
  cInput ->SetBranchAddress("jetPhi", &stage2_jetPhi, &b_stage2_jetPhi);
  cInput ->SetBranchAddress("jetEt", &stage2_jetEt, &b_stage2_jetEt);
  cInput ->SetBranchAddress("jetBx", &stage2_jetBx, &b_stage2_jetBx);
  
  cInput ->SetBranchAddress("tauEta", &stage2_tauEta, &b_stage2_tauEta);
  cInput ->SetBranchAddress("tauPhi", &stage2_tauPhi, &b_stage2_tauPhi);
  cInput ->SetBranchAddress("tauEt", &stage2_tauEt, &b_stage2_tauEt);
  cInput ->SetBranchAddress("tauBx", &stage2_tauBx, &b_stage2_tauBx);
  cInput ->SetBranchAddress("tauIso", &stage2_tauIso, &b_stage2_tauIso);
  
  
  ugtInput ->SetBranchAddress("m_algoDecisionInitial", &trigger_algo, &b_algoDecisionInitial);
  
  
  TFile* fOutVBF = new TFile (Form("%s/%s",directory.Data(),fOutNameVBF.Data()), "recreate");
  
  fOutVBF->cd(); 
  //VBF
  
  TH2D* DiJet2D_Pass = new TH2D ("DiJet2D", "events firing, 2-dim (lead, mjj)", 800, 0, 800, 100,30, 130);
  TH2D* Ratio_DiJet2D = new TH2D ("Ratio_DiJet2D", "VBF 2-dim ratio", 800, 0, 800, 100, 30, 130);
  
  
  TH2D* DiJet2D_Sub_Pass = new TH2D ("DiJet2D_Sub", "events firing, 2-dim (lead, sublead)", 30, 30, 60, 100,30, 130);
  TH2D* Ratio_Sub_DiJet2D = new TH2D ("Ratio_Sub_DiJet2D", "VBF 2-dim ratio", 30, 30, 60, 100, 30, 130);
  
  
  TH2D* DiJet2D_Pass_rej = new TH2D ("DiJet2D_rej", "events firing, 2-dim (lead, mjj), no TT28", 800, 0, 800, 100,30, 130);
  TH2D* Ratio_DiJet2D_rej = new TH2D ("Ratio_DiJet2D_rej", "VBF 2-dim ratio", 800, 0, 800, 100, 30, 130);
  
  
  TH2D* DiJet2D_Sub_Pass_rej = new TH2D ("DiJet2D_Sub_rej", "events firing, 2-dim (lead, sublead), no TT28", 30, 30, 60, 100,30, 130);
  TH2D* Ratio_Sub_DiJet2D_rej = new TH2D ("Ratio_Sub_DiJet2D_rej", "VBF 2-dim ratio", 30, 30, 60, 100, 30, 130);
  
  
  TH1D* VBF_lead = new TH1D ("VBF_lead", "events firing, 1-dim (lead jet)", 100,50, 150);
  TH1D* Ratio_VBF_lead = new TH1D ("Ratio_VBF_lead", "VBF 1-dim ratio", 100,50, 150);
  
  
  TH1D* VBFpure_lead = new TH1D ("VBFpure_lead", "events firing, 1-dim (lead jet)", 100,50, 150);
  TH1D* Ratio_VBFpure_lead = new TH1D ("Ratio_VBFpure_lead", "VBF 1-dim ratio", 100,50, 150);
  
  
  TH1D* VBF_mjj = new TH1D ("VBF_mjj", "events firing, 1-dim (mjj)", 800, 0, 800);
  TH1D* Ratio_VBF_mjj = new TH1D ("Ratio_VBF_mjj", "VBF 1-dim ratio", 800, 0, 800);

  TH1D* Npass = new TH1D ("Npass","Npass",1,0,1);
  TH1D* uGT = new TH1D ("uGT", "ugt rate", 1,0,1);
  
  TH1D* Mjj30 = new TH1D ("Mjj30", "", 100, 0, 800);
  TH1D* jetsRes = new TH1D ("jetsRes", "", 60, 0, 30);
  TH1D* Mjj35 = new TH1D ("Mjj35", "", 100, 0, 800);
  TH1D* Mjj40 = new TH1D ("Mjj40", "", 100, 0, 800);
  TH1D* Mjj45 = new TH1D ("Mjj45", "", 100, 0, 800);
  TH1D* Mjj50 = new TH1D ("Mjj50", "", 100, 0, 800);
  TH1D* Mjj55 = new TH1D ("Mjj55", "", 100, 0, 800);

   
 
  std::map<Int_t,Float_t> PU_per_LS;
  std::string str; 
  while (std::getline(PUFile, str))
    {
      TString temp(str);
      
      regex reg(Form("%d,%d,",Fill,Run));
      temp = regex_replace(str, reg, "");
      int pos_coma = temp.First(",");
      TString LS_str(temp,pos_coma);
      // cout<<LS_str<<endl;
      TString Replacing = LS_str ;
      Replacing += ",";
      temp.ReplaceAll(Replacing.Data(),"");
      TString PU_str = temp;
      //      cout<<PU_str<<endl;
      std::istringstream ss_LS(LS_str.Data());
      Int_t LS ;
      ss_LS >> LS;
      std::istringstream ss_PU(PU_str.Data());
      Float_t PU ;
      ss_PU >> PU;     
      PU_per_LS.insert(std::pair<Int_t,Float_t>(LS , PU ));
      
    }
  // analyze data    
  long int nEvents = 0;

  int nEventsPass = 0;


  
  std::vector<object> jet30;
  std::vector<object> tau32;   
  std::vector<object> jet30rej;   
  

  
  std::vector< std::tuple<double,int,int> > mjj_pass; //VBF
  std::vector< std::tuple<double,int,int,double, double> > mjj_pass_sortPt; //VBF
  std::vector< std::tuple<double,int,int> > mjj_rej; //VBF
  std::vector< std::tuple<double,int,int,double, double> > mjj_rej_sortPt; //VBF
  
  bool L1_DoubleIsoTau32er = false;
  
  for (Long64_t iEv =0 ; true ; ++iEv){
    lumi = 0;			   
    run = 0;
    
    stage2_jetN =0;		   
    stage2_jetEt.clear(); 
    stage2_jetEta.clear();
    stage2_jetPhi.clear();
    stage2_jetBx.clear();
    
    
    stage2_tauEt.clear(); 
    stage2_tauEta.clear();
    stage2_tauPhi.clear();
    stage2_tauBx.clear();
    stage2_tauIso.clear();
    
    
    trigger_algo.clear();
    int got = 0;
    
    lInput->GetEntry(iEv);
    got = cInput->GetEntry(iEv);
    

    if (got == 0) break;
    ugtInput->GetEntry(iEv);
   
    
    
    if (iEv%1000 == 0) cout << iEv << " / " << nEvents << endl;
    
    jet30.clear();
    tau32.clear();
    jet30rej.clear();
    mjj_pass.clear();
    mjj_pass_sortPt.clear();
    mjj_rej.clear();
    mjj_rej_sortPt.clear();
    

    
    

    if(lumi<lumimin || lumi>lumimax) continue; 
    
    if(run>Run) continue;
    
    if(PU_per_LS.find(lumi)==PU_per_LS.end()) continue;
    Float_t weight = 0;
    if (reweight){
      weight = PU_per_LS[56]/PU_per_LS[lumi];
    }else{
      weight =  1.;
    }
    
    nEventsPass ++;

    Npass->Fill(0.);
    cout <<"@@@@@@"<< nEventsPass<<endl;
    
    if(trigger_ID!=-999){
      cout <<"before"<<endl;
      if (trigger_algo[trigger_ID]){
	uGT->Fill(0);
      }else{
	uGT->Fill(-1);
      }
     cout <<"after"<<endl;
    }
   
    for (long int iL1 = 0; iL1 < stage2_jetEt.size(); iL1++){ //loop on jets
      // selections
      double jetPt  = stage2_jetEt.at(iL1);
      double jetEta  = fabs(stage2_jetEta.at(iL1));
      int jetBx = fabs(stage2_jetBx.at(iL1));
      cout <<"jets "<<jetPt<<" "<<jetBx<<" "<<jetEta<<" "<<endl;
      if(jetPt>subleadjet_cut && jetBx == 0) {
	jet30.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	if(jetEta>2.7 && jetEta<3.0){
	  if(jetPt>60.) jet30rej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	}else{
	  jet30rej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	}
      }
      
    }

    for (long int iL1 = 0; iL1 < stage2_tauEt.size(); iL1++){ //loop on taus
      // selections
      
      double tauPt  = stage2_tauEt.at(iL1);
      double tauEta  = fabs(stage2_tauEta.at(iL1));
      double tauIso  =stage2_tauIso.at(iL1);
      int tauBx = fabs(stage2_tauBx.at(iL1));
      cout <<"taus"<<tauPt<<" "<<tauBx<<" "<<tauEta<<" "<<tauIso<<endl;
      if(tauPt>32. && tauBx == 0 && tauEta<2.1 && tauIso>0.5) {
	tau32.push_back(object(stage2_tauEt.at(iL1),stage2_tauEta.at(iL1),stage2_tauPhi.at(iL1),stage2_tauIso.at(iL1))) ;
      }
      
    }
    


    std::sort (jet30.begin(),jet30.end());
    std::sort (tau32.begin(),tau32.end());
    std::sort (jet30rej.begin(),jet30rej.end());
    
    L1_DoubleIsoTau32er = false;
    if (tau32.size()>1)     L1_DoubleIsoTau32er =true;
    if (tau32.size()>1)     cout<<tau32[1].Et()<< "########### tau pass"<<endl;
    cout<<endl;
    //VBF
    
    
    if (jet30.size() >= 2){
      for (int iJet = 0; iJet <jet30.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet30.size(); kJet++){      
	  if (kJet!=iJet) {
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jet30[iJet].Et(),
			      jet30[iJet].Eta(),
			      jet30[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jet30[kJet].Et(),
			      jet30[kJet].Eta(),
			      jet30[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    cout<<"mjj "<<jetPair.M()<<" k "<<jet30[kJet].Et()<<" i "<<jet30[iJet].Et()<<endl;
	    mjj_pass.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    if(jetPair.M()>=620)    mjj_pass_sortPt.push_back(make_tuple(jetPair.M(),iJet,kJet,jet30[iJet].Et(),jet30[kJet].Et()));
	    double jetDiff = jet30[iJet].Et()-jet30[kJet].Et();
	    double jetSum  = jet30[iJet].Et()+jet30[kJet].Et();
	    jetsRes->Fill(jetSum/jetDiff);
	    Mjj30->Fill(jetPair.M());
	    if(jet30[kJet].Et()>35) Mjj35->Fill(jetPair.M());
	    if(jet30[kJet].Et()>40) Mjj40->Fill(jetPair.M());
	    if(jet30[kJet].Et()>45) Mjj45->Fill(jetPair.M());
	    if(jet30[kJet].Et()>50) Mjj50->Fill(jetPair.M());
	    if(jet30[kJet].Et()>55) Mjj55->Fill(jetPair.M());
	    
	  }
	  
	}
	
      }
      std::sort(mjj_pass.begin(),mjj_pass.end());
      std::sort(mjj_pass_sortPt.begin(),mjj_pass_sortPt.end(),SortMjjByJetThreshold);
      double Mjj_max = std::get<0>(*(mjj_pass.rbegin()));      
      DiJet2D_Pass -> Fill (Mjj_max,jet30[0].Et(),weight);        
      if(mjj_pass_sortPt.size()>0) {
	double subJet_min = std::get<4>(*(mjj_pass_sortPt.rbegin()));
	if (Mjj_max>Mjj_cut){ cout <<"### mjjpass"<<endl; DiJet2D_Sub_Pass -> Fill (subJet_min,jet30[0].Et(),weight); } 
	if (subJet_min>subleadjet_cut && Mjj_max>Mjj_cut){ cout<<"### sub pass"<<endl;  VBF_lead->Fill(jet30[0].Et(),weight);}

	if (subJet_min>subleadjet_cut && Mjj_max>Mjj_cut && (!L1_DoubleIsoTau32er)) VBFpure_lead->Fill(jet30[0].Et(),weight);
	
	if (subJet_min>subleadjet_cut && jet30[0].Et()>leadjet_cut) VBF_mjj->Fill(Mjj_max,weight); 
      }else{
	DiJet2D_Sub_Pass->Fill(-1,-1);
      }
    } else {
      DiJet2D_Pass -> Fill (-1, -1);
      DiJet2D_Sub_Pass->Fill(-1,-1);
    }


    //VBF rejecting TT28
    
    if (jet30rej.size() >= 2){
      for (int iJet = 0; iJet <jet30rej.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet30rej.size(); kJet++){      
	  if (kJet!=iJet) {
	    
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jet30rej[iJet].Et(),
			      jet30rej[iJet].Eta(),
			      jet30rej[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jet30rej[kJet].Et(),
			      jet30rej[kJet].Eta(),
			      jet30rej[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_rej.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    if(jetPair.M()>=620) mjj_rej_sortPt.push_back(make_tuple(jetPair.M(),iJet,kJet,jet30rej[iJet].Et(),jet30rej[kJet].Et()));
	  }
	}
      }
      std::sort(mjj_rej.begin(),mjj_rej.end());
      std::sort(mjj_rej_sortPt.begin(),mjj_rej_sortPt.end(),SortMjjByJetThreshold);
      double Mjj_max =  std::get<0>(*(mjj_rej.rbegin()));
      DiJet2D_Pass_rej -> Fill ( Mjj_max,jet30rej[0].Et(),weight);  
      if(mjj_rej_sortPt.size()>0) {
	double subJet_min =  std::get<4>(*(mjj_rej_sortPt.rbegin()));
	DiJet2D_Sub_Pass_rej -> Fill (subJet_min,jet30rej[0].Et(),weight);  
      }else{
	DiJet2D_Sub_Pass_rej->Fill(-1,-1);
      }
    } else {
      DiJet2D_Pass_rej -> Fill (-1, -1);
      DiJet2D_Sub_Pass_rej->Fill(-1,-1);
    }  
      
  }

  // compute ratio plots
  //scaling to be made afterwards!
  
    
  cout << endl;
  cout << "Computing ratio..." << endl; 

  //VBF
  if(nEventsPass>0){
  for (int i = 1; i <=DiJet2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiJet2D_Pass->GetNbinsY();j++ ){
      double binDiJet2D = 1.*(DiJet2D_Pass->Integral(i, DiJet2D_Pass->GetNbinsX()+1,j,DiJet2D_Pass->GetNbinsY()+1));
      Ratio_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
            
      //rej
      binDiJet2D = 1.*(DiJet2D_Pass_rej->Integral(i, DiJet2D_Pass_rej->GetNbinsX()+1,j,DiJet2D_Pass_rej->GetNbinsY()+1));
      Ratio_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        
           
    }
  }
  
  for (int j = 1; j<= VBF_lead->GetNbinsX();j++ ){
    double binDiJet = 1.*(VBF_lead->Integral(j, VBF_lead->GetNbinsX()+1));
      Ratio_VBF_lead -> SetBinContent (j, binDiJet);
  }

  for (int j = 1; j<= VBFpure_lead->GetNbinsX();j++ ){
    double binDiJet = 1.*(VBFpure_lead->Integral(j, VBFpure_lead->GetNbinsX()+1));
      Ratio_VBFpure_lead -> SetBinContent (j, binDiJet);
  }
   for (int j = 1; j<= VBF_mjj->GetNbinsX();j++ ){
     double binDiJet = 1.*(VBF_mjj->Integral(j, VBF_mjj->GetNbinsX()+1));
      Ratio_VBF_mjj -> SetBinContent (j, binDiJet);        
  }
  for (int i = 1; i <=DiJet2D_Sub_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiJet2D_Sub_Pass->GetNbinsY();j++ ){
      double binDiJet2D = 1.*(DiJet2D_Sub_Pass->Integral(i, DiJet2D_Sub_Pass->GetNbinsX()+1,j,DiJet2D_Sub_Pass->GetNbinsY()+1));
      Ratio_Sub_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      //rej
      binDiJet2D = 1.*(DiJet2D_Sub_Pass_rej->Integral(i, DiJet2D_Sub_Pass_rej->GetNbinsX()+1,j,DiJet2D_Sub_Pass_rej->GetNbinsY()+1));
      Ratio_Sub_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        
    }
  }

  }
  fOutVBF->cd();
  fOutVBF -> Write();
  cout <<"Output saved in "<<fOutNameVBF<<endl;
}

