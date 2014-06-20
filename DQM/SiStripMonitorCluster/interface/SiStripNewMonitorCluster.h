#ifndef SiStripMonitorCluster_SiStripNewMonitorCluster_h
#define SiStripMonitorCluster_SiStripNewMonitorCluster_h
// -*- C++ -*-
// Package:     SiStripMonitorCluster
// Class  :     SiStripNewMonitorCluster

// Original Author:  hbakhshi
//         Created:  Sat Oct  26 12:37:04 CET 2013

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQM/SiStripCommon/interface/SiStripNewFolderOrganizer.h"

#include <vector>

class DQMStore;
class SiStripDetCabling;
class SiStripCluster;
class SiStripDCSStatus;
class GenericTriggerEventFlag;

class SiStripNewMonitorCluster : public edm::EDAnalyzer {
 public:
  explicit SiStripNewMonitorCluster(const edm::ParameterSet&);
  ~SiStripNewMonitorCluster();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void beginJob() ;
  virtual void endJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);

  struct ModMEs{ // MEs for one single detector module

    MonitorElement* NumberOfClusters = 0;
    MonitorElement* ClusterPosition = 0;
    MonitorElement* ClusterDigiPosition = 0;
    MonitorElement* ClusterWidth = 0;
    MonitorElement* ClusterCharge = 0;
    MonitorElement* ClusterNoise = 0;
    MonitorElement* ClusterSignalOverNoise = 0;
    MonitorElement* ClusterSignalOverNoiseVsPos = 0;
    MonitorElement* ModuleLocalOccupancy = 0;
    MonitorElement* NrOfClusterizedStrips = 0; // can be used at client level for occupancy calculations
  };

  struct ClusterProperties { // Cluster Properties
    float charge;
    float position;
    short start;
    short width;
    float noise;
  };


 private:

  void createMEs(const edm::EventSetup& es);
  void createModuleMEs(ModMEs& mod_single, uint32_t detid);
  //  int FindRegion(int nstrip,int npixel);
  void fillModuleMEs(ModMEs& mod_mes, ClusterProperties& cluster);


  inline void fillME(MonitorElement* ME,float value1){if (ME!=0)ME->Fill(value1);}
  inline void fillME(MonitorElement* ME,float value1,float value2){if (ME!=0)ME->Fill(value1,value2);}
  inline void fillME(MonitorElement* ME,float value1,float value2,float value3){if (ME!=0)ME->Fill(value1,value2,value3);}
  inline void fillME(MonitorElement* ME,float value1,float value2,float value3,float value4){if (ME!=0)ME->Fill(value1,value2,value3,value4);}
  MonitorElement * bookMETrend(const char*, const char*);
  MonitorElement* bookME1D(const char* ParameterSetLabel, const char* HistoName);

 private:
  DQMStore* dqmStore_;
  unsigned long long m_cacheID_;
  SiStripNewFolderOrganizer folder_organizer;
  edm::ParameterSet conf_;
  // flags


  edm::ESHandle<SiStripDetCabling> SiStripDetCabling_;

  int runNb, eventNb;


  edm::InputTag clusterProducerStrip_;

  std::string qualityLabel_;

  bool   applyClusterQuality_;
  double sToNLowerLimit_;  
  double sToNUpperLimit_;  
  double widthLowerLimit_;
  double widthUpperLimit_;

  edm::ParameterSet Parameters;

  // add for selecting on ZeroBias events in the MinimumBias PD
  GenericTriggerEventFlag* genTriggerEventFlagBPTXfilter_;
  GenericTriggerEventFlag* genTriggerEventFlagPixelDCSfilter_;
  GenericTriggerEventFlag* genTriggerEventFlagStripDCSfilter_;

  bool passBPTXfilter_;
  bool passPixelDCSfilter_;
  bool passStripDCSfilter_;
};
#endif
