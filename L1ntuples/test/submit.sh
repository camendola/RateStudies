#ZeroBias Events, PU55

./submitOnTier3_L1Ntuplizer.py -o /data_CMS/cms/amendola/RateStudiesL1Ntuples/ -t ZeroBias28Feb2017HighPU_ibx0_BunchTrain1_2016H9Nov -s fileLists/Data_ZeroBiasBunchTrain1_2016H9Nov.py -w True -n 20
./submitOnTier3_L1Ntuplizer.py -o /data_CMS/cms/amendola/RateStudiesL1Ntuples/ -t ZeroBias28Feb2017HighPU_ibx0_BunchTrain2_2016H9Nov -s fileLists/Data_ZeroBiasBunchTrain2_2016H9Nov.py -w True -n 20
./submitOnTier3_L1Ntuplizer.py -o /data_CMS/cms/amendola/RateStudiesL1Ntuples/ -t ZeroBias28Feb2017HighPU_ibx0_BunchTrain3_2016H9Nov -s fileLists/Data_ZeroBiasBunchTrain3_2016H9Nov.py -w True -n 20
./submitOnTier3_L1Ntuplizer.py -o /data_CMS/cms/amendola/RateStudiesL1Ntuples/ -t ZeroBias28Feb2017HighPU_ibx0_BunchTrain4_2016H9Nov -s fileLists/Data_ZeroBiasBunchTrain4_2016H9Nov.py -w True -n 20
./submitOnTier3_L1Ntuplizer.py -o /data_CMS/cms/amendola/RateStudiesL1Ntuples/ -t ZeroBias28Feb2017HighPU_ibx0_BunchTrain5_2016H9Nov -s fileLists/Data_ZeroBiasBunchTrain5_2016H9Nov.py -w True -n 20

#mkdir /data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias28Feb2017HighPU_ibx0_BunchTrain1-5_2016H9Nov_-1Events
#hadd /data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias28Feb2017HighPU_ibx0_BunchTrain1-5_2016H9Nov_-1Events/L1total.root /data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias28Feb2017HighPU_ibx0_BunchTrain*_2016H9Nov_-1Events/L1ntuples_*.root


