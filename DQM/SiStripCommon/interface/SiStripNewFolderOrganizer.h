#ifndef SiStripCommon_SiStripNewFolderOrganizer_h
#define SiStripCommon_SiStripNewFolderOrganizer_h
// -*- C++ -*-
//
// Package:     SiStripCommon
// Class  :     SiStripNewFolderOrganizer
// 
/**\class SiStripNewFolderOrganizer SiStripNewFolderOrganizer.h DQM/SiStripCommon/interface/SiStripNewFolderOrganizer.h

   Description: <Organizes the folders for the monitoring elements of the SiStrip Tracker. Its methods return strings with names of folders to be created and used.>

   Usage:
   <usage>

*/
//
// Original Author:  hbakhshi
//         Created:  Thu Oct 17 12:11:45 CET 2013

//
#include "CommonTools/UtilAlgos/interface/DetIdSelector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>

using namespace std;
using namespace edm;

class DQMStore;
class MonitorElement;


class SiStripNewFolderOrganizer
{
public:
  SiStripNewFolderOrganizer( edm::ParameterSet parent_dir , const SiStripNewFolderOrganizer* TheParent = NULL );

  void NavigateHere();

  bool Contains(const unsigned int& detid) const; //search only this directory
  SiStripNewFolderOrganizer* SearchInSubDirectories( const unsigned int& detid ); 
  
  void AddME( string name , MonitorElement* );
  void AddME( const unsigned int& detid , string name , MonitorElement* me);

  MonitorElement* GetME( string name );
  MonitorElement* GetME( const unsigned int& detid , string name);



  string GetTabs() const;
  string GetFullName() const;
  void Print() const;

private:
  std::string Name;
  std::vector<std::string> DetIdRanges;
  DetIdSelector detIdSelector;

  edm::VParameterSet SubDirInfo;
  std::vector< SiStripNewFolderOrganizer* > SubDirs;

  const SiStripNewFolderOrganizer* Parent;

  DQMStore* dbe_;

  map< string , MonitorElement* > SummaryMEs;
  map< unsigned int ,  map< string , MonitorElement* > > ModuleMEs;
};

/*
  public:
  static unsigned short const all_ = 65535;

  SiStripNewFolderOrganizer();
  virtual ~SiStripNewFolderOrganizer();

  // top folder
  void setSiStripNewFolderName(std::string name);
  std::string getSiStripNewFolder();
  void setSiStripNewFolder();

  // control folder
  std::string getSiStripTopControlFolder();
  void setSiStripTopControlFolder();
  std::string getSiStripControlFolder(
  // unsigned short crate,
  unsigned short slot = all_,
  unsigned short ring = all_,
  unsigned short addr = all_,
  unsigned short chan = all_
  // unsigned short i2c
  );
  void setSiStripControlFolder(
  // unsigned short crate,
  unsigned short slot = all_,
  unsigned short ring = all_,
  unsigned short addr = all_,
  unsigned short chan = all_
  // unsigned short i2c
  );

  std::pair<std::string,int32_t> GetSubDetAndLayer(const uint32_t& detid, const TrackerTopology* tTopo, bool ring_flag = 0);
  // detector folders
  void setDetectorFolder(uint32_t rawdetid, const TrackerTopology* tTopo);
  void getFolderName(int32_t rawdetid, const TrackerTopology* tTopo, std::string& lokal_folder);
  void getFolderName(int32_t rawdetid, std::string& lokal_folder);  // deprecated version, still needed for now

  // layer folders
  void setLayerFolder(uint32_t rawdetid,const TrackerTopology* tTopo,int32_t layer=0,bool ring_flag = 0);
  void getLayerFolderName(std::stringstream& ss, uint32_t rawdetid, const TrackerTopology* tTopo, bool ring_flag = 0);
  void getSubDetLayerFolderName(std::stringstream& ss, SiStripDetId::SubDetector subDet, uint32_t layer, uint32_t side=0);
  // SubDetector Folder
  void getSubDetFolder(const uint32_t& detid, const TrackerTopology* tTopo, std::string& folder_name);
  std::pair<std::string, std::string> getSubDetFolderAndTag(const uint32_t& detid, const TrackerTopology* tTopo);
  private:
  SiStripNewFolderOrganizer(const SiStripNewFolderOrganizer&); // stop default
  const SiStripNewFolderOrganizer& operator=(const SiStripNewFolderOrganizer&); // stop default

  private:
  std::string TopFolderName;
  DQMStore* dbe_;
  };
*/
#endif
