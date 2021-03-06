// evaluates rate  - compile with c++ -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`
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
#include "interface/object.h"

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

 
  //  2016
  //   float nbStudiedRun = 224;
  //2017
  //float nbStudiedRun = 96;
  //float scale = 0.001*(nbStudiedRun*11245.6);

    float nbStudiedRun =96.;
  float thisLumiRun = 1750E30;
  float scaleToLumi = 2.E34;
  float scale = 0.001*(nbStudiedRun*11245.6)*scaleToLumi/thisLumiRun;  
  
  float add = 0;



  TString LumiTarget = "scaleToLumi20E33";
  TString year = "2017";
  
  cout << "Scale factor: " << scale << endl;
  //ZeroBias sample L1
  //  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/ZeroBiasRun2017A/";
  //2016
  //    TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias26Apr2017HighPU_ibx0_BunchTrain0-5_2016H9Nov_-1Events/";
  //2017
  //    TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesHighPU_2017_fill5824/";
        TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesHighPU_2017_fill6194/";
  //  TString directory = "/eos/user/c/camendol/local/";
  //TFile *file = TFile::Open(Form("%sL1total.root", directory.Data()),"read");
   //2016
    //  TString fileList = "fileLists/L1NtuplesL1Menu2017ZeroBiasRun2016H_BunchTrainsX.txt";
  //2017
  //TString fileList = "fileLists/L1NtuplesL1Menu2017ZeroBiasFill5824.txt";
	TString fileList = "fileLists/ZeroBias2017_HighPU.list";
  
    
    bool emulated  =false;
  bool reweight = false;

  //     ifstream PUFile("utils/PU_per_LS.txt"); //2016
  //    ifstream PUFile("utils/PU_per_LS_fill5824_2017.txt"); //2017
      ifstream PUFile("utils/PU_per_LS_fill6194_2017.txt"); //2017

  
  TTree * tInput;
  TChain * cInput;
  TChain * lInput; 



  
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

  TString fOutNameVBF;
  TString fOutNameXtrigger;
  TString fOutNameTaus;
  TString fOutNameTausBjet;

  if(emulated){
    fOutNameVBF = directory+"Emu_rateL1_VBF_"+LumiTarget+"_"+year+".root";
    fOutNameXtrigger = directory+"Emu_rateL1_xtrigger_"+LumiTarget+"_"+year+".root";
    fOutNameTaus = directory+"Emu_rateL1_taus_"+LumiTarget+"_"+year+".root";
    fOutNameTausBjet = directory+"Emu_rateL1_ditau_bjet_"+LumiTarget+"_"+year+".root";
  }else{
    fOutNameVBF = directory+"rateL1_VBF_"+LumiTarget+"_"+year+".root";
    fOutNameXtrigger = directory+"rateL1_xtrigger_"+LumiTarget+"_"+year+".root";
    fOutNameTaus = directory+"rateL1_taus_"+LumiTarget+"_"+year+".root";
    fOutNameTausBjet = directory+"rateL1_ditau_bjet_"+LumiTarget+"_"+year+".root";
  }
  
  
  lInput->SetMakeClass(1);
  cInput->SetMakeClass(1);  
  
  Int_t lumi ;
  UInt_t run ;
  UShort_t stage2_tauN ;
  std::vector<float> stage2_tauEt;
  std::vector<float> stage2_tauEta;
  std::vector<float> stage2_tauPhi;
  std::vector<short> stage2_tauIso;  


  UShort_t stage2_muonN ;
  std::vector<float> stage2_muonEt;
  std::vector<float> stage2_muonEta;
  std::vector<float> stage2_muonPhi;
  std::vector<short> stage2_muonIso;  

  UShort_t stage2_jetN ;
  std::vector<float> stage2_jetEt;
  std::vector<float> stage2_jetEta;
  std::vector<float> stage2_jetPhi;

  
  
  /*
    UInt_t lumi ;
    UShort_t stage2_tauN ;
    std::vector<float>* stage2_tauEt=0;
    std::vector<float>* stage2_tauEta=0;
    std::vector<float>* stage2_tauPhi=0;
    std::vector<short> *stage2_tauIso=0;  


    UShort_t stage2_muonN ;
    std::vector<float>* stage2_muonEt=0;
    std::vector<float>* stage2_muonEta=0;
    std::vector<float>* stage2_muonPhi=0;
    std::vector<short> *stage2_muonIso=0;  
  
    UShort_t stage2_jetN ;
    std::vector<float> *stage2_jetEt=0;
    std::vector<float> *stage2_jetEta=0;
    std::vector<float> *stage2_jetPhi=0;
  */
  
  // set branch and variables
  TBranch *b_lumi ;
  TBranch *b_run ;

  TBranch *b_stage2_tauN ;
  TBranch *b_stage2_tauEt;
  TBranch *b_stage2_tauEta;
  TBranch *b_stage2_tauPhi;
  TBranch *b_stage2_tauIso;

  TBranch *b_stage2_muonN ;
  TBranch *b_stage2_muonEt;
  TBranch *b_stage2_muonEta;
  TBranch *b_stage2_muonPhi;
  TBranch *b_stage2_muonIso;
  
  TBranch *b_stage2_jetN ;
  TBranch *b_stage2_jetEt;
  TBranch *b_stage2_jetEta;
  TBranch *b_stage2_jetPhi;
  
  

  

    lInput ->SetBranchAddress("lumi", &lumi , &b_lumi );
    lInput ->SetBranchAddress("run", &run , &b_run );
    cInput ->SetBranchAddress("nTaus", &stage2_tauN , &b_stage2_tauN );
    cInput ->SetBranchAddress("tauEta", &stage2_tauEta, &b_stage2_tauEta);
    cInput ->SetBranchAddress("tauPhi", &stage2_tauPhi, &b_stage2_tauPhi);
    cInput ->SetBranchAddress("tauEt", &stage2_tauEt, &b_stage2_tauEt);
    cInput ->SetBranchAddress("tauIso", &stage2_tauIso, &b_stage2_tauIso);
    cInput ->SetBranchAddress("nMuons", &stage2_muonN , &b_stage2_muonN );
    cInput ->SetBranchAddress("muonEta", &stage2_muonEta, &b_stage2_muonEta);
    cInput ->SetBranchAddress("muonPhi", &stage2_muonPhi, &b_stage2_muonPhi);
    cInput ->SetBranchAddress("muonEt", &stage2_muonEt, &b_stage2_muonEt);
    cInput ->SetBranchAddress("muonIso", &stage2_muonIso, &b_stage2_muonIso);
    cInput ->SetBranchAddress("nJets", &stage2_jetN , &b_stage2_jetN);
    cInput ->SetBranchAddress("jetEta", &stage2_jetEta, &b_stage2_jetEta);
    cInput ->SetBranchAddress("jetPhi", &stage2_jetPhi, &b_stage2_jetPhi);
    cInput ->SetBranchAddress("jetEt", &stage2_jetEt, &b_stage2_jetEt);

  
  ///////////////
  //// HISTO ////
  ///////////////
  
  TFile* fOutVBF = new TFile (Form("%s",fOutNameVBF.Data()), "recreate");
  TFile* fOutTaus = new TFile (Form("%s",fOutNameTaus.Data()), "recreate");
  TFile* fOutXtrigger = new TFile (Form("%s",fOutNameXtrigger.Data()), "recreate");
  TFile* fOutTausBjet = new TFile (Form("%s",fOutNameTausBjet.Data()), "recreate");
  fOutTaus->cd(); 
  //taus 
  TH1D* SubleadTauPt_Pass = new TH1D ("SubleadTauPt_Pass", "SubleadTauPt_Pass", 100, 0, 100);
  TH1D* SubleadTauPt_Boost20_Pass = new TH1D ("SubleadTauPt_Boost20_Pass", "SubleadTauPt_Boost20_Pass", 100, 0, 100);
  TH1D* SubleadTauPt_Boost30_Pass = new TH1D ("SubleadTauPt_Boost30_Pass", "SubleadTauPt_Boost30_Pass", 100, 0, 100);
  TH1D* SubleadTauPt_Boost40_Pass = new TH1D ("SubleadTauPt_Boost40_Pass", "SubleadTauPt_Boost40_Pass", 100, 0, 100);
  TH1D* SubleadTauPt_Boost50_Pass = new TH1D ("SubleadTauPt_Boost50_Pass", "SubleadTauPt_Boost50_Pass", 100, 0, 100);
  //  TH1D* SubleadTauPt_Pass_jet = new TH1D ("SubleadTauPt_Pass_jet", "SubleadTauPt_Pass", 100, 0, 100);
  //  TH1D* DiTauPt_Pass = new TH1D ("DiTauPt_Pass", "DiTauPt_Pass", 100, 0, 100);
  TH1D* Rate_diTauSublead = new TH1D ("Rate_diTauSublead", "rate - double tau", 100, 0, 100);
  TH1D* Rate_diTauSublead_Boost20 = new TH1D ("Rate_diTauSublead_Boost20", "Rate_diTauSublead_Boost20", 100, 0, 100);
  TH1D* Rate_diTauSublead_Boost30 = new TH1D ("Rate_diTauSublead_Boost30", "Rate_diTauSublead_Boost30", 100, 0, 100);
  TH1D* Rate_diTauSublead_Boost40 = new TH1D ("Rate_diTauSublead_Boost40", "Rate_diTauSublead_Boost40", 100, 0, 100);
  TH1D* Rate_diTauSublead_Boost50 = new TH1D ("Rate_diTauSublead_Boost50", "Rate_diTauSublead_Boost50", 100, 0, 100);
  TH1D* Ratio_diTauSublead = new TH1D ("Ratio_diTauSublead", "ratio - double tau", 100, 0, 100);
  TH1D* Ratio_diTauSublead_Boost20 = new TH1D ("Ratio_diTauSublead_Boost20", "Ratio_diTauSublead_Boost20", 100, 0, 100);
  TH1D* Ratio_diTauSublead_Boost30 = new TH1D ("Ratio_diTauSublead_Boost30", "Ratio_diTauSublead_Boost30", 100, 0, 100);
  TH1D* Ratio_diTauSublead_Boost40 = new TH1D ("Ratio_diTauSublead_Boost40", "Ratio_diTauSublead_Boost40", 100, 0, 100);
  TH1D* Ratio_diTauSublead_Boost50 = new TH1D ("Ratio_diTauSublead_Boost50", "Ratio_diTauSublead_Boost50", 100, 0, 100);
  //  TH1D* Rate_diTauSublead_jet = new TH1D ("Rate_diTauSublead_jet", "rate - double tau", 100, 0, 100);
  //  TH1D* Ratio_diTauSublead_jet = new TH1D ("Ratio_diTauSublead_jet", "ratio - double tau", 100, 0, 100);
  TH2D* DiTau2D_Pass = new TH2D ("DiTau2D_Pass", "double tau",  400, 0,400, 100, 25,125); 
  TH2D* Ratio_DiTau2D = new TH2D ("Ratio_DiTau2D", "Ratio - 2 taus", 400, 0,400, 100, 25,125); 
  TH2D* Rate_DiTau2D = new TH2D ("Rate_DiTau2D", "Rate - 2 taus", 400, 0,400, 100, 25,125); 
  /*    TH2D* PureDiTau2D_Pass_wrt30 = new TH2D ("PureDiTau2D_Pass_wrt30", "double tau",  400, 0,400, 100, 25,125); 
	TH2D* PureRatio_DiTau2D_wrt30 = new TH2D ("PureRatio_DiTau2D_wrt30", "Ratio - 2 taus", 400, 0,400, 100, 25,125); 
	TH2D* PureRate_DiTau2D_wrt30 = new TH2D ("PureRate_DiTau2D_wrt30", "Rate - 2 taus", 400, 0,400, 100, 25,125);
	TH2D* PureDiTau2D_Pass_wrt31 = new TH2D ("PureDiTau2D_Pass_wrt31", "double tau",  400, 0,400, 100, 25,125); 
	TH2D* PureRatio_DiTau2D_wrt31 = new TH2D ("PureRatio_DiTau2D_wrt31", "Ratio - 2 taus", 400, 0,400, 100, 25,125); 
	TH2D* PureRate_DiTau2D_wrt31 = new TH2D ("PureRate_DiTau2D_wrt31", "Rate - 2 taus", 400, 0,400, 100, 25,125); 
  */

  TH2D* PureDiTau2D_Pass[8] ;
  TH2D* PureRatio_DiTau2D[8];
  TH2D* PureRate_DiTau2D[8] ;
  for (int xx=0; xx<8; xx++){
    PureDiTau2D_Pass[xx]  = new TH2D ("PureDiTau2D_Pass", "pass pure rate",  400, 0,400, 100, 25,125); 
    PureRatio_DiTau2D[xx] = new TH2D ("PureRatio_DiTau2D", "Ratio - 2 taus pass pure rate", 400, 0,400, 100, 25,125); 
    PureRate_DiTau2D[xx] = new TH2D ("PureRate_DiTau2D", "Rate - 2 taus pass pure rate", 400, 0,400, 100, 25,125);
  }
  TH1D* PtTauTau = new TH1D ("PtTauTau", "", 100, 0, 400);
  TH1D* MTauTau = new TH1D ("MTauTau", "", 100, 0, 400);
  
  fOutVBF->cd(); 
  //VBF

  TH2D* DiJet2D_Pass = new TH2D ("DiJet2D", "double jet const", 800, 0, 800, 100,30, 130);
  TH2D* Ratio_DiJet2D = new TH2D ("Ratio_DiJet2D", "Ratio - 2 jets, E_{T} > 30 GeV ", 800, 0, 800, 100, 30, 130);
  TH2D* Rate_DiJet2D = new TH2D ("Rate_DiJet2D", "Rate - 2 jets, E_{T} > 30 GeV ", 800, 0, 800, 100, 30, 130);

  TH2D* DiJet2D_Sub_Pass = new TH2D ("DiJet2D_Sub", "double jet const", 30, 30, 60, 100,30, 130);
  TH2D* Ratio_Sub_DiJet2D = new TH2D ("Ratio_Sub_DiJet2D", "Ratio - 2 jets, E_{T} > 30 GeV ", 30, 30, 60, 100, 30, 130);
  TH2D* Rate_Sub_DiJet2D = new TH2D ("Rate_Sub_DiJet2D", "Rate - 2 jets, E_{T} > 30 GeV ", 30, 30, 60, 100, 30, 130); 

  TH2D* DiJet2D_Pass_rej = new TH2D ("DiJet2D_rej", "double jet const", 800, 0, 800, 100,30, 130);
  TH2D* Ratio_DiJet2D_rej = new TH2D ("Ratio_DiJet2D_rej", "Ratio - 2 jets, E_{T} > 30 GeV ", 800, 0, 800, 100, 30, 130);
  TH2D* Rate_DiJet2D_rej = new TH2D ("Rate_DiJet2D_rej", "Rate - 2 jets, E_{T} > 30 GeV ", 800, 0, 800, 100, 30, 130);

  TH2D* DiJet2D_Sub_Pass_rej = new TH2D ("DiJet2D_Sub_rej", "double jet const", 30, 30, 60, 100,30, 130);
  TH2D* Ratio_Sub_DiJet2D_rej = new TH2D ("Ratio_Sub_DiJet2D_rej", "Ratio - 2 jets, E_{T} > 30 GeV ", 30, 30, 60, 100, 30, 130);
  TH2D* Rate_Sub_DiJet2D_rej = new TH2D ("Rate_Sub_DiJet2D_rej", "Rate - 2 jets, E_{T} > 30 GeV ", 30, 30, 60, 100, 30, 130); 
  
  TH1D* VBF_lead = new TH1D ("VBF_lead", "double jet const", 100,30, 130);
  TH1D* Ratio_VBF_lead = new TH1D ("Ratio_VBF_lead", "double jet const", 100,30, 130);
  TH1D* Rate_VBF_lead = new TH1D ("Rate_VBF_lead", "double jet const", 100,30, 130);
  
  TH1D* Mjj30 = new TH1D ("Mjj30", "", 100, 0, 800);
  TH1D* jetsRes = new TH1D ("jetsRes", "", 60, 0, 30);
  TH1D* Mjj35 = new TH1D ("Mjj35", "", 100, 0, 800);
  TH1D* Mjj40 = new TH1D ("Mjj40", "", 100, 0, 800);
  TH1D* Mjj45 = new TH1D ("Mjj45", "", 100, 0, 800);
  TH1D* Mjj50 = new TH1D ("Mjj50", "", 100, 0, 800);
  TH1D* Mjj55 = new TH1D ("Mjj55", "", 100, 0, 800);

  /* TH1D* ETA = new TH1D ("ETA", "", 1000, -4, 4);
  TH1D* ETA_strip1 = new TH1D ("ETA_strip400", "", 1000, -4, 4);
  TH1D* ETA_strip2 = new TH1D ("ETA_strip467", "", 1000, -4, 4);
  
  TH2D* DiJet2D_Pass_strips = new TH2D ("DiJet2D_strips", "check", 100, 0, 800, 100,30, 130);
  TH2D* DiJet2D_Pass_strips_nolead = new TH2D ("DiJet2D_strips_nolead", "check", 100, 0, 800, 100,30, 130);
  */


  
  fOutXtrigger->cd(); 
  //Di-Tau+jet  
  TH2D* DiTauJet2D_Pass = new TH2D ("DiTauJet2D", "double tau + jet", 100, 0, 100, 100,0, 100);
  TH2D* Ratio_DiTauJet2D = new TH2D ("Ratio_DiTauJet2D", "Ratio - 2 taus, jet E_{T} > 30 GeV", 100, 0, 100, 100, 0, 100);
  TH2D* Rate_DiTauJet2D = new TH2D ("Rate_DiTauJet2D", "Rate - 2 taus, jet E_{T} > 30 GeV ", 100, 0, 100, 100, 0, 100); 
  TH2D* PureDiTauJet2D_Pass_jet = new TH2D ("PureDiTauJet2D_jet", "double tau + jet", 100, 0, 100, 100,0, 100);
  TH2D* PureRatio_DiTauJet2D_jet = new TH2D ("PureRatio_DiTauJet2D_jet", "Ratio - 2 taus, jet E_{T} > 30 GeV", 100, 0, 100, 100, 0, 100);
  TH2D* PureRate_DiTauJet2D_jet = new TH2D ("PureRate_DiTauJet2D_jet", "Rate - 2 taus, jet E_{T} > 30 GeV ", 100, 0, 100, 100, 0, 100); 

  TH2D* DiTauJet2D_Pass_jet = new TH2D ("DiTauJet2D_jet", "double tau + jet", 100, 0, 100, 100,30, 130);
  TH2D* Ratio_DiTauJet2D_jet = new TH2D ("Ratio_DiTauJet2D_jet", "Ratio - 2 taus, jet E_{T} > 30 GeV", 100, 0, 100, 100, 30, 130);
  TH2D* Rate_DiTauJet2D_jet = new TH2D ("Rate_DiTauJet2D_jet", "Rate - 2 taus, jet E_{T} > 30 GeV ", 100, 0, 100, 100, 30, 130); 
  TH1D* IsoTau = new TH1D ("IsoTau", "", 20, 0, 10);

  fOutTausBjet->cd();
  //DiTau + bjet
  TH2D* DiTauBJet2D_Pass = new TH2D ("DiTauJetB2D", "double tau + bjet", 100, 0, 100, 100,0, 100);
  TH2D* Ratio_DiTauBJet2D = new TH2D ("Ratio_DiTauBJet2D", "Ratio - 2 taus + bjet", 100, 0, 100, 100, 0, 100);
  TH2D* Rate_DiTauBJet2D = new TH2D ("Rate_DiTauBJet2D", "Rate - 2 taus + bjet ", 100, 0, 100, 100, 0, 100); 


  std::map<Int_t,Float_t> PU_per_LS;
  std::string str; 
  while (std::getline(PUFile, str))
    {
      TString temp(str);

      //temp.ReplaceAll("5412,283171,",""); //2016

             regex reg("6194,[0-9]+,");
          temp = regex_replace(str, reg, "");
      int pos_coma = temp.First(",");
      TString LS_str(temp,pos_coma);
      // cout<<LS_str<<endl;
      TString Replacing = LS_str ;
      Replacing += ",";
      temp.ReplaceAll(Replacing.Data(),"");
      TString PU_str = temp;
      // cout<<PU_str<<endl;
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


  std::vector<object> tau;
  std::vector<object> muon;
  std::vector<object> tauNoOverlap;
  std::vector<object> tau25noOverlap;
  std::vector<object> jet30;   
  std::vector<object> jet30rej;   
  std::vector<object> jet30noOverlap;   

  std::vector< std::tuple<double,int,int> > et_ditau_pass; //et of tau pair
  std::vector< std::tuple<double,int,int> > et_ditau_all; //et of tau pair  
  std::vector< std::tuple<double,int,int> > m_ditau_pass; //m of tau pair  
  std::vector< std::tuple<double,int,int> > mjj_pass; //VBF
  std::vector< std::tuple<double,int,int,double, double> > mjj_pass_sortPt; //VBF
  std::vector< std::tuple<double,int,int> > mjj_rej; //VBF
  std::vector< std::tuple<double,int,int,double, double> > mjj_rej_sortPt; //VBF


  double singletau = 0;
  double ditau25 = 0;  
  double ditaujet = 0;  
  double ditau = 0;  
  double ditau25jet = 0;  

  bool L1_DoubleIsoTau25er_Jet50= false;
  bool L1_DoubleIsoTauXXer[4];
  bool L1_DoubleIsoTau25er_PtTauTau70 = false;
  
  
  
  int Nfill = 0;  
  
  for (Long64_t iEv =0 ;iEv<30000; ++iEv){
    lumi = 0;			   
    run = 0;
    stage2_tauN=0;		   
    stage2_tauEt.clear(); 
    stage2_tauEta.clear();
    stage2_tauPhi.clear();
    stage2_tauIso.clear();    
		   
		   
    stage2_muonN=0;		   
    stage2_muonEt.clear();    
    stage2_muonEta.clear();
    stage2_muonPhi.clear();
    stage2_muonIso.clear();   
		   
    stage2_jetN =0;		   
    stage2_jetEt.clear(); 
    stage2_jetEta.clear();
    stage2_jetPhi.clear();
    int got = 0;
   
      lInput->GetEntry(iEv);
      got = cInput->GetEntry(iEv);
    
    if (got == 0) break;

   
  
    // if(lumi>781) continue;
    //if(lumi>150) continue;
    
    if (iEv%1000 == 0) cout << iEv << " / " << nEvents << endl;
   





    
    jet30.clear();
    jet30rej.clear();
    jet30noOverlap.clear();
    tauNoOverlap.clear();
    tau25noOverlap.clear();
    mjj_pass.clear();
    mjj_pass_sortPt.clear();
    mjj_rej.clear();
    mjj_rej_sortPt.clear();
    et_ditau_pass.clear();
    et_ditau_all.clear();
    m_ditau_pass.clear();
    tau.clear(); 
    muon.clear();


    
    //    if(lumi<56 || lumi>69) continue; //2016

    if(run>302674) continue;
     if(lumi<135 || lumi>264) continue; //2017
     
     if(PU_per_LS.find(lumi)==PU_per_LS.end()) continue;
     Float_t weight = 0;
    if (reweight){
      weight = PU_per_LS[56]/PU_per_LS[lumi];
    }else{
      weight =  1.;
    }

    nEventsPass ++;
    L1_DoubleIsoTau25er_PtTauTau70 = false;
    L1_DoubleIsoTau25er_Jet50= false;
    for (int xx= 0; xx<8; xx++){
      L1_DoubleIsoTauXXer[xx]= false;
    }


        for (long int iL1 = 0; iL1 < stage2_tauN; iL1++){ //loop on taus
    //    for (long int iL1 = 0; iL1 < stage2_tauEt.size(); iL1++){ //loop on taus
      // selections
      double tauEta  = stage2_tauEta.at(iL1);
      double tauIso  = stage2_tauIso.at(iL1);

      if(tauIso>0.5 && fabs(tauEta)<2.1) tau.push_back(object(stage2_tauEt.at(iL1),stage2_tauEta.at(iL1),stage2_tauPhi.at(iL1),stage2_tauIso.at(iL1))) ;
     
    }


    for (long int iL1 = 0; iL1 < stage2_jetEt.size(); iL1++){ //loop on jets
      // selections
      double jetPt  = stage2_jetEt.at(iL1);
      double jetEta  = fabs(stage2_jetEta.at(iL1));

      if(jetPt>30.) {
	jet30.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	//	if(jetEta>2.7 && jetEta<3.0){
	//  if(jetPt>60.) jet30rej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	
	if(jetEta<3){
	   jet30rej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	   //}else{
	   // jet30rej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	  }
      }

    }

    ///    for (long int iL1 = 0; iL1 < stage2_muonN; iL1++){ //loop on muons
      for (long int iL1 = 0; iL1 < stage2_muonEt.size(); iL1++){ //loop on muons
	muon.push_back(object(stage2_muonEt.at(iL1),stage2_muonEta.at(iL1),stage2_muonPhi.at(iL1),stage2_muonIso.at(iL1))) ;
 }

    std::sort (jet30.begin(),jet30.end());
    std::sort (jet30rej.begin(),jet30rej.end());
    std::sort (tau.begin(),tau.end());
    std::sort (muon.begin(),muon.end());

    //overap removal
    tauNoOverlap.push_back(object(tau[0].Et(),tau[0].Eta(),tau[0].Phi(),tau[0].Iso())) ;      
    tau25noOverlap.push_back(object(tau[0].Et(),tau[0].Eta(),tau[0].Phi(),tau[0].Iso())) ;      
    for (int iTau =1;iTau<tau.size();iTau++){
      if (tau[iTau].DeltaR(tau[0])>0.2) tauNoOverlap.push_back(object(tau[iTau].Et(),tau[iTau].Eta(),tau[iTau].Phi(),tau[iTau].Iso())) ;      
      if (tau[iTau].DeltaR(tau[0])>0.2 && tau[iTau].Et()>25) tau25noOverlap.push_back(object(tau[iTau].Et(),tau[iTau].Eta(),tau[iTau].Phi(),tau[iTau].Iso())) ;      
    }
    for (int iJet =0;iJet<jet30.size();iJet++){
      if ((jet30[iJet].DeltaR(tau[0])>0.2)&&(jet30[iJet].DeltaR(tau[1])>0.2)) jet30noOverlap.push_back(object(jet30[iJet].Et(),jet30[iJet].Eta(),jet30[iJet].Phi(),-999)) ;
        
    }
    
  
    
    std::sort (jet30noOverlap.begin(),jet30noOverlap.end());//not necessary
    std::sort (tauNoOverlap.begin(),tauNoOverlap.end());//not necessary
    std::sort (tau25noOverlap.begin(),tau25noOverlap.end());//not necessary

   
       
    //DiTau+PtDiTau
    if (tauNoOverlap.size() >= 2){
      int Ntau = tauNoOverlap.size();
      if (Ntau > 5) Ntau=5; 
      for (int iTau = 0; iTau <Ntau; iTau++){      
	for (int kTau = 0; kTau <Ntau; kTau++){      
	  if (kTau!=iTau) {
	    TLorentzVector itau;
	    itau.SetPtEtaPhiM(
			      tauNoOverlap[iTau].Et(),
			      tauNoOverlap[iTau].Eta(),
			      tauNoOverlap[iTau].Phi(),
			      0.);
	    TLorentzVector ktau;
	    ktau.SetPtEtaPhiM(
			      tauNoOverlap[kTau].Et(),
			      tauNoOverlap[kTau].Eta(),
			      tauNoOverlap[kTau].Phi(),
			      0.);
	    TLorentzVector tauPair = itau+ktau;
	    et_ditau_all.push_back(make_tuple(tauPair.Et(),iTau,kTau));

	     
	  }
	}
      }
      std::sort(et_ditau_all.begin(),et_ditau_all.end());
    }
    if (tau25noOverlap.size() >= 2){
      int Ntau = tau25noOverlap.size();
      if (Ntau > 5) Ntau=5; 
      for (int iTau = 0; iTau <Ntau; iTau++){      
	for (int kTau = 0; kTau <Ntau; kTau++){      
	  if (kTau!=iTau) {
	    TLorentzVector itau;
	    itau.SetPtEtaPhiM(
			      tau25noOverlap[iTau].Et(),
			      tau25noOverlap[iTau].Eta(),
			      tau25noOverlap[iTau].Phi(),
			      0.);
	    TLorentzVector ktau;
	    ktau.SetPtEtaPhiM(
			      tau25noOverlap[kTau].Et(),
			      tau25noOverlap[kTau].Eta(),
			      tau25noOverlap[kTau].Phi(),
			      0.);
	    TLorentzVector tauPair = itau+ktau;
	    et_ditau_pass.push_back(make_tuple(tauPair.Et(),iTau,kTau));
	    m_ditau_pass.push_back(make_tuple(tauPair.M(),iTau,kTau));
	     
	  }
	}
      }
      std::sort(et_ditau_pass.begin(),et_ditau_pass.end());
      std::sort(m_ditau_pass.begin(),m_ditau_pass.end());
      PtTauTau->Fill(std::get<0>(*(et_ditau_pass.rbegin())));
      MTauTau->Fill(std::get<0>(*(m_ditau_pass.rbegin())));


      for (int xx=0; xx<8; xx++){	
	if (tauNoOverlap[1].Et()>(float)(30+xx)) L1_DoubleIsoTauXXer[xx] = true;
      }
      
      if (tau25noOverlap[1].Et()>25){
	if(std::get<0>(*(et_ditau_pass.rbegin()))>70) L1_DoubleIsoTau25er_PtTauTau70 = true;

      }
      DiTau2D_Pass -> Fill ( std::get<0>(*(et_ditau_pass.rbegin())),tau25noOverlap[1].Et(),weight);
      for (int xx=0; xx<8; xx++){
	if(!L1_DoubleIsoTauXXer[xx]) {
	  	  PureDiTau2D_Pass[xx] -> Fill ( std::get<0>(*(et_ditau_pass.rbegin())),tau25noOverlap[1].Et(),weight);
	  if(tau25noOverlap[1].Et()>33){
	    // cout<<"plot "<<xx<<"; pass wrt"<<xx +30<<endl;
	  }

	}else{
	    PureDiTau2D_Pass[xx] -> Fill ( -1,-1);
	  
	}
      }
      
    }else{
      DiTau2D_Pass -> Fill (-1,-1);
      for (int xx=0; xx<8; xx++){
       PureDiTau2D_Pass[xx] -> Fill ( -1,-1);
      }
    }
      

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
     
      
      DiJet2D_Pass -> Fill ( std::get<0>(*(mjj_pass.rbegin())),jet30[0].Et(),weight);        
      if(mjj_pass_sortPt.size()>0) {

	if (std::get<0>(*(mjj_pass_sortPt.rbegin()))>620) DiJet2D_Sub_Pass -> Fill ( std::get<4>(*(mjj_pass_sortPt.rbegin())),jet30[0].Et(),weight);  
	if (std::get<4>(*(mjj_pass_sortPt.rbegin()))>35 && std::get<0>(*(mjj_pass.rbegin()))>620) VBF_lead->Fill(jet30[0].Et(),weight); 
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
     
      DiJet2D_Pass_rej -> Fill ( std::get<0>(*(mjj_rej.rbegin())),jet30rej[0].Et(),weight);  
      
      if(mjj_rej_sortPt.size()>0) {
	DiJet2D_Sub_Pass_rej -> Fill ( std::get<4>(*(mjj_rej_sortPt.rbegin())),jet30rej[0].Et(),weight);  
	
      }else{
       DiJet2D_Sub_Pass_rej->Fill(-1,-1);
      }

      
    } else {
      DiJet2D_Pass_rej -> Fill (-1, -1);
      DiJet2D_Sub_Pass_rej->Fill(-1,-1);
    }  
    

    //Taus and ditau+jet
    if (tauNoOverlap.size() >= 2 )
      {
	SubleadTauPt_Pass -> Fill(tauNoOverlap[1].Et(),weight);
	if(et_ditau_all.size()>0){
	  if(std::get<0>(*(et_ditau_all.rbegin()))>20){
	    SubleadTauPt_Boost20_Pass-> Fill(tauNoOverlap[1].Et(),weight);
	  }else{
	    SubleadTauPt_Boost20_Pass-> Fill(-1);
	  }
	  if(std::get<0>(*(et_ditau_all.rbegin()))>30){
	    SubleadTauPt_Boost30_Pass-> Fill(tauNoOverlap[1].Et(),weight);
	  }else{
	    SubleadTauPt_Boost30_Pass-> Fill(-1);
	  }
	  if(std::get<0>(*(et_ditau_all.rbegin()))>40) {
	    SubleadTauPt_Boost40_Pass-> Fill(tauNoOverlap[1].Et(),weight);
	  }else{
	    SubleadTauPt_Boost40_Pass-> Fill(-1);
	  }
	  
	  if(std::get<0>(*(et_ditau_all.rbegin()))>50) {
	    SubleadTauPt_Boost50_Pass-> Fill(tauNoOverlap[1].Et(),weight);
	  }else{
	    SubleadTauPt_Boost50_Pass-> Fill(-1);
	  }
	}else{
	  SubleadTauPt_Boost20_Pass-> Fill(-1);
	  SubleadTauPt_Boost30_Pass-> Fill(-1);
	  SubleadTauPt_Boost40_Pass-> Fill(-1);
	  SubleadTauPt_Boost50_Pass-> Fill(-1);
	}
	IsoTau -> Fill(tauNoOverlap[0].Iso());
	if(muon.size()>0){
	  DiTauBJet2D_Pass->Fill(tauNoOverlap[1].Et(),muon[0].Et(),weight); 
	}else{
	  DiTauBJet2D_Pass->Fill(-1,-1);
	}
	
	if (jet30noOverlap.size()>0){

	  DiTauJet2D_Pass->Fill( tauNoOverlap[0].Et(),tauNoOverlap[1].Et(),weight );
	  DiTauJet2D_Pass_jet->Fill( tauNoOverlap[1].Et(),jet30noOverlap[0].Et(),weight);

	}else{
	  SubleadTauPt_Pass -> Fill(-1);
	  SubleadTauPt_Boost20_Pass-> Fill(-1);
	  SubleadTauPt_Boost30_Pass-> Fill(-1);
	  SubleadTauPt_Boost40_Pass-> Fill(-1);
	  SubleadTauPt_Boost50_Pass-> Fill(-1);
	  DiTauJet2D_Pass->Fill( -1,-1 );
	  DiTauJet2D_Pass_jet->Fill( -1,-1);
	  DiTauBJet2D_Pass->Fill(-1,-1);
	}
      }else if(tauNoOverlap.size() == 1){
      SubleadTauPt_Pass -> Fill ( -1 );
      SubleadTauPt_Boost20_Pass-> Fill(-1);
      SubleadTauPt_Boost30_Pass-> Fill(-1);
      SubleadTauPt_Boost40_Pass-> Fill(-1);
      SubleadTauPt_Boost50_Pass-> Fill(-1);
      
      DiTauJet2D_Pass->Fill( -1,-1 );
      DiTauJet2D_Pass_jet->Fill( -1,-1 );
      DiTauBJet2D_Pass->Fill(-1,-1);
      
    }else{
      SubleadTauPt_Pass -> Fill ( -1 );        
      SubleadTauPt_Boost20_Pass-> Fill(-1);
      SubleadTauPt_Boost30_Pass-> Fill(-1);
      SubleadTauPt_Boost40_Pass-> Fill(-1);
      SubleadTauPt_Boost50_Pass-> Fill(-1);
      DiTauJet2D_Pass->Fill( -1,-1 );       
      DiTauJet2D_Pass_jet->Fill( -1,-1 );
      DiTauBJet2D_Pass->Fill(-1,-1);
    }       
    
  
  
  
  
  
  
  
  }

  // compute rate plots




  //cout << endl;
  cout << "Computing rates..." << endl; 
  //taus rates
  for (int i = 1; i <= SubleadTauPt_Pass->GetNbinsX(); i++){
    double binDiTauSublead = 1.*(SubleadTauPt_Pass->Integral(i,SubleadTauPt_Pass->GetNbinsX()+1))/nEventsPass;
    Ratio_diTauSublead -> SetBinContent (i, binDiTauSublead);

    binDiTauSublead *=scale;
    Rate_diTauSublead -> SetBinContent (i, binDiTauSublead);
    binDiTauSublead = 1.*(SubleadTauPt_Boost20_Pass->Integral(i,SubleadTauPt_Boost20_Pass->GetNbinsX()+1))/nEventsPass;
    Ratio_diTauSublead_Boost20 -> SetBinContent (i, binDiTauSublead);
    binDiTauSublead *=scale;
    Rate_diTauSublead_Boost20 -> SetBinContent (i, binDiTauSublead);
    binDiTauSublead = 1.*(SubleadTauPt_Boost30_Pass->Integral(i,SubleadTauPt_Boost30_Pass->GetNbinsX()+1))/nEventsPass;
    Ratio_diTauSublead_Boost30 -> SetBinContent (i, binDiTauSublead);
    binDiTauSublead *=scale;
    Rate_diTauSublead_Boost30 -> SetBinContent (i, binDiTauSublead);
    binDiTauSublead = 1.*(SubleadTauPt_Boost40_Pass->Integral(i,SubleadTauPt_Boost40_Pass->GetNbinsX()+1))/nEventsPass;
    Ratio_diTauSublead_Boost40 -> SetBinContent (i, binDiTauSublead);
    binDiTauSublead *=scale;
    Rate_diTauSublead_Boost40 -> SetBinContent (i, binDiTauSublead);
    binDiTauSublead = 1.*(SubleadTauPt_Boost50_Pass->Integral(i,SubleadTauPt_Boost50_Pass->GetNbinsX()+1))/nEventsPass;
    Ratio_diTauSublead_Boost50 -> SetBinContent (i, binDiTauSublead);        
    binDiTauSublead *=scale;
    Rate_diTauSublead_Boost50 -> SetBinContent (i, binDiTauSublead);        
    
    

  }

  

  //VBF
  for (int i = 1; i <=DiJet2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiJet2D_Pass->GetNbinsY();j++ ){

      double binDiJet2D = 1.*(DiJet2D_Pass->Integral(i, DiJet2D_Pass->GetNbinsX()+1,j,DiJet2D_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      double addTerm = 0.;
     
      binDiJet2D *=scale;
      Rate_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      //rej
      binDiJet2D = 1.*(DiJet2D_Pass_rej->Integral(i, DiJet2D_Pass_rej->GetNbinsX()+1,j,DiJet2D_Pass_rej->GetNbinsY()+1))/nEventsPass;
      Ratio_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        
      addTerm=binDiJet2D*add;
      binDiJet2D =(binDiJet2D+addTerm)*scale;
      Rate_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        
    }
  }

    for (int j = 1; j<= VBF_lead->GetNbinsX();j++ ){
      double addTerm = 0.;
      double binDiJet = 1.*(VBF_lead->Integral(j, VBF_lead->GetNbinsX()+1))/nEventsPass;
      Ratio_VBF_lead -> SetBinContent (j, binDiJet);        
      addTerm=binDiJet*add;
      binDiJet =(binDiJet+addTerm)*scale;

      Rate_VBF_lead ->  SetBinContent (j, binDiJet);
    }
  
    for (int i = 1; i <=DiJet2D_Sub_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiJet2D_Sub_Pass->GetNbinsY();j++ ){
      double addTerm = 0.;
      double binDiJet2D = 1.*(DiJet2D_Sub_Pass->Integral(i, DiJet2D_Sub_Pass->GetNbinsX()+1,j,DiJet2D_Sub_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_Sub_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      addTerm=binDiJet2D*add;

      binDiJet2D =(binDiJet2D+addTerm)*scale;

      
      Rate_Sub_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      //rej
      
      binDiJet2D = 1.*(DiJet2D_Sub_Pass_rej->Integral(i, DiJet2D_Sub_Pass_rej->GetNbinsX()+1,j,DiJet2D_Sub_Pass_rej->GetNbinsY()+1))/nEventsPass;
      Ratio_Sub_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        

      addTerm=binDiJet2D*add;
      binDiJet2D =(binDiJet2D+addTerm)*scale;

      Rate_Sub_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        
    }
  }

  //ditau + ptditau
  for (int i = 1; i <=DiTau2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiTau2D_Pass->GetNbinsY();j++ ){
      double binDiTau2D = 1.*(DiTau2D_Pass->Integral(i, DiTau2D_Pass->GetNbinsX()+1,j,DiTau2D_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_DiTau2D -> SetBinContent (i, j, binDiTau2D);        
      binDiTau2D *=scale;
      Rate_DiTau2D -> SetBinContent (i, j, binDiTau2D);        
      for (int xx=0; xx<8; xx++){
	PureRatio_DiTau2D[xx]->SetName(Form("PureRatio_DiTau2D_wrt%d",30+xx));
	PureRate_DiTau2D[xx]->SetName(Form("PureRate_DiTau2D_wrt%d",30+xx));
	binDiTau2D = 1.*(PureDiTau2D_Pass[xx]->Integral(i, PureDiTau2D_Pass[xx]->GetNbinsX()+1,j,PureDiTau2D_Pass[xx]->GetNbinsY()+1))/nEventsPass;
	PureRatio_DiTau2D[xx] -> SetBinContent (i, j, binDiTau2D);        
	binDiTau2D *=scale;
	PureRate_DiTau2D[xx] -> SetBinContent (i, j, binDiTau2D);
	}
    }
  }
    
  //diTau+jet
  for (int i = 1; i <=DiTauJet2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiTauJet2D_Pass->GetNbinsY();j++ ){
      double binDiTauJet2D = 1.*(DiTauJet2D_Pass->Integral(i, DiTauJet2D_Pass->GetNbinsX()+1,j,DiTauJet2D_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_DiTauJet2D -> SetBinContent (i, j, binDiTauJet2D);        
      binDiTauJet2D *=scale;
      Rate_DiTauJet2D -> SetBinContent (i, j, binDiTauJet2D);        

      
    }
  }
  
  for (int i = 1; i <=DiTauJet2D_Pass_jet->GetNbinsX(); i++){
    for (int j = 1; j<= DiTauJet2D_Pass_jet->GetNbinsY();j++ ){
      double binDiTauJet2D_jet = 1.*(DiTauJet2D_Pass_jet->Integral(i, DiTauJet2D_Pass_jet->GetNbinsX()+1,j,DiTauJet2D_Pass_jet->GetNbinsY()+1))/nEventsPass;
      Ratio_DiTauJet2D_jet -> SetBinContent (i, j, binDiTauJet2D_jet);        
      binDiTauJet2D_jet *=scale;
      Rate_DiTauJet2D_jet -> SetBinContent (i, j, binDiTauJet2D_jet);        
      //      binDiTauJet2D_jet = 1.*(PureDiTauJet2D_Pass_jet->Integral(i, PureDiTauJet2D_Pass_jet->GetNbinsX()+1,j,PureDiTauJet2D_Pass_jet->GetNbinsY()+1))/nEventsPass;
      // PureRatio_DiTauJet2D_jet -> SetBinContent (i, j, binDiTauJet2D_jet);        
      // binDiTauJet2D_jet *=scale;
      //PureRate_DiTauJet2D_jet -> SetBinContent (i, j, binDiTauJet2D_jet);              
    }
  }
  
  //diTau+bjet
  for (int i = 1; i <=DiTauBJet2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiTauBJet2D_Pass->GetNbinsY();j++ ){
      double binDiTauBJet2D = 1.*(DiTauBJet2D_Pass->Integral(i, DiTauBJet2D_Pass->GetNbinsX()+1,j,DiTauBJet2D_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_DiTauBJet2D -> SetBinContent (i, j, binDiTauBJet2D);        
      binDiTauBJet2D *=scale;
      Rate_DiTauBJet2D -> SetBinContent (i, j, binDiTauBJet2D);        

      
    }
  }
  
  int xbin = Rate_DiTau2D->GetXaxis()->FindBin(70.0);
  int ybin = Rate_DiTau2D->GetYaxis()->FindBin(25.0);
  cout<<"Rate L1_DoubleIsoTau25er_PtTauTau70   "<<Rate_DiTau2D->GetBinContent(xbin,ybin)<<" kHz"<<endl; 

  
  xbin = PureRate_DiTau2D[0]->GetXaxis()->FindBin(70.0);
  ybin = PureRate_DiTau2D[0]->GetYaxis()->FindBin(25.0);
  cout<<"Pure Rate L1_DoubleIsoTau25er_PtTauTau70  wrt 30 "<<PureRate_DiTau2D[0]->GetBinContent(xbin,ybin)<<" kHz"<<endl; 

  cout<<"old "<<endl;
  xbin = Rate_DiJet2D->GetXaxis()->FindBin(620.0);
  ybin = Rate_DiJet2D->GetYaxis()->FindBin(90.0);
  cout<<" Rate VBFseed 90 30 620  "<<Rate_DiJet2D->GetBinContent(xbin,ybin)<<" kHz"<<endl;

  
  cout<<"new "<<endl;
  xbin = Rate_Sub_DiJet2D->GetXaxis()->FindBin(35.0);
  ybin = Rate_Sub_DiJet2D->GetYaxis()->FindBin(110.0);
  cout<<" Rate VBFseed 110 35 620  "<<Rate_Sub_DiJet2D->GetBinContent(xbin,ybin)<<" kHz"<<endl;
  xbin = Rate_Sub_DiJet2D->GetXaxis()->FindBin(30.0);
  ybin = Rate_Sub_DiJet2D->GetYaxis()->FindBin(90.0);
  cout<<" Rate VBFseed 90 30 620  "<<Rate_Sub_DiJet2D->GetBinContent(xbin,ybin)<<" kHz"<<endl; 
 

  fOutTaus->cd();
   
  fOutTaus -> Write();
  fOutVBF->cd();
  fOutVBF -> Write();
  fOutXtrigger->cd();
  fOutXtrigger->Write();
  fOutTausBjet->cd();
  fOutTausBjet->Write();
}

