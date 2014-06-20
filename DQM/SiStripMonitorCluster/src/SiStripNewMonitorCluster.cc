#include <vector>
#include <numeric>
#include <fstream>
#include <math.h>
#include "TNamed.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
#include "CalibFormats/SiStripObjects/interface/SiStripGain.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DQM/SiStripCommon/interface/SiStripNewFolderOrganizer.h"
#include "DQM/SiStripCommon/interface/SiStripHistoId.h"
#include "DQM/SiStripMonitorCluster/interface/SiStripNewMonitorCluster.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "CalibTracker/SiStripCommon/interface/SiStripDCSStatus.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DPGAnalysis/SiStripTools/interface/APVCyclePhaseCollection.h"
#include "DPGAnalysis/SiStripTools/interface/EventWithHistory.h"

#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"

#include "TMath.h"
#include <iostream>

//--------------------------------------------------------------------------------------------
SiStripNewMonitorCluster::SiStripNewMonitorCluster(const edm::ParameterSet& iConfig):
  dqmStore_(edm::Service<DQMStore>().operator->()),
  m_cacheID_(0),
  folder_organizer( iConfig.getParameter< edm::ParameterSet >("DirectoryStructure") ),
  conf_(iConfig)
{
  genTriggerEventFlagBPTXfilter_     = new GenericTriggerEventFlag(iConfig.getParameter<edm::ParameterSet>("BPTXfilter")     );
  genTriggerEventFlagPixelDCSfilter_ = new GenericTriggerEventFlag(iConfig.getParameter<edm::ParameterSet>("PixelDCSfilter") );
  genTriggerEventFlagStripDCSfilter_ = new GenericTriggerEventFlag(iConfig.getParameter<edm::ParameterSet>("StripDCSfilter") );

  eventNb = 0;


  qualityLabel_  = conf_.getParameter<std::string>("StripQualityLabel");
  clusterProducerStrip_ = conf_.getParameter<edm::InputTag>("ClusterProducerStrip");

  edm::ParameterSet cluster_condition = conf_.getParameter<edm::ParameterSet>("ClusterConditions");
  applyClusterQuality_ = cluster_condition.getParameter<bool>("On");
  sToNLowerLimit_      = cluster_condition.getParameter<double>("minStoN");
  sToNUpperLimit_      = cluster_condition.getParameter<double>("maxStoN");
  widthLowerLimit_     = cluster_condition.getParameter<double>("minWidth"); 
  widthUpperLimit_     = cluster_condition.getParameter<double>("maxWidth"); 

  cout << "folder organizer" << endl;
  folder_organizer.Print();
  cout << "folder organizer" << endl;

} 

SiStripNewMonitorCluster::~SiStripNewMonitorCluster() { 
  if (genTriggerEventFlagBPTXfilter_    ) delete genTriggerEventFlagBPTXfilter_;
  if (genTriggerEventFlagPixelDCSfilter_) delete genTriggerEventFlagPixelDCSfilter_;
  if (genTriggerEventFlagStripDCSfilter_) delete genTriggerEventFlagStripDCSfilter_;
}

//--------------------------------------------------------------------------------------------
void SiStripNewMonitorCluster::beginRun(const edm::Run& run, const edm::EventSetup& es){

  // Initialize the GenericTriggerEventFlag
  if ( genTriggerEventFlagBPTXfilter_->on() )
    genTriggerEventFlagBPTXfilter_->initRun( run, es );
  if ( genTriggerEventFlagPixelDCSfilter_->on() )
    genTriggerEventFlagPixelDCSfilter_->initRun( run, es );
  if ( genTriggerEventFlagStripDCSfilter_->on() )
    genTriggerEventFlagStripDCSfilter_->initRun( run, es );

  unsigned long long cacheID = es.get<SiStripDetCablingRcd>().cacheIdentifier();
  if (m_cacheID_ != cacheID) {
    m_cacheID_ = cacheID;       
    createMEs(es);
  }
}

//--------------------------------------------------------------------------------------------
void SiStripNewMonitorCluster::createMEs(const edm::EventSetup& es){

  // take from eventSetup the SiStripDetCabling object - here will use SiStripDetControl later on
  es.get<SiStripDetCablingRcd>().get(SiStripDetCabling_);
    
  // get list of active detectors from SiStripDetCabling 
  std::vector<uint32_t> activeDets;
  SiStripDetCabling_->addActiveDetectorsRawIds(activeDets);
    

  // loop over detectors and book MEs
  edm::LogInfo("SiStripTkDQM|SiStripNewMonitorCluster")<<"nr. of activeDets:  "<<activeDets.size();
  for(std::vector<uint32_t>::iterator detid_iterator = activeDets.begin(); detid_iterator!=activeDets.end(); detid_iterator++){
    uint32_t detid = (*detid_iterator);
    // remove any eventual zero elements - there should be none, but just in case
    if(detid == 0) {
      activeDets.erase(detid_iterator);
      continue;
    }
      
    ModMEs mod_single;
    SiStripNewFolderOrganizer* folder = folder_organizer.SearchInSubDirectories( detid );
    if( folder != NULL ){
      folder->NavigateHere() ;
      createModuleMEs(mod_single, detid);
      folder->AddME( detid ,  "NumberOfClusters" , mod_single.NumberOfClusters);
      folder->AddME( detid ,  "ClusterPosition" , mod_single.ClusterPosition);
      folder->AddME( detid ,  "ClusterDigiPosition" , mod_single.ClusterDigiPosition);
      folder->AddME( detid ,  "ClusterWidth" , mod_single.ClusterWidth);
      folder->AddME( detid ,  "ClusterCharge" , mod_single.ClusterCharge);
      folder->AddME( detid ,  "ClusterNoise" , mod_single.ClusterNoise);
      folder->AddME( detid ,  "ClusterSignalOverNoise" , mod_single.ClusterSignalOverNoise);
      folder->AddME( detid ,  "ClusterSignalOverNoiseVsPos" , mod_single.ClusterSignalOverNoiseVsPos);
      folder->AddME( detid ,  "ModuleLocalOccupancy" , mod_single.ModuleLocalOccupancy);
      folder->AddME( detid ,  "NrOfClusterizedStrips" , mod_single.NrOfClusterizedStrips);
    }//end of loop over detectors
  }    
}//end of method

//--------------------------------------------------------------------------------------------
void SiStripNewMonitorCluster::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Filter out events if Trigger Filtering is requested
  passBPTXfilter_     = ( iEvent.isRealData() and genTriggerEventFlagBPTXfilter_->on()     ) ? genTriggerEventFlagBPTXfilter_->accept( iEvent, iSetup)     : true;
  passPixelDCSfilter_ = ( iEvent.isRealData() and genTriggerEventFlagPixelDCSfilter_->on() ) ? genTriggerEventFlagPixelDCSfilter_->accept( iEvent, iSetup) : true;
  passStripDCSfilter_ = ( iEvent.isRealData() and genTriggerEventFlagStripDCSfilter_->on() ) ? genTriggerEventFlagStripDCSfilter_->accept( iEvent, iSetup) : true;

  runNb   = iEvent.id().run();
  eventNb++;

  //int NStripClusters=0;

  edm::ESHandle<SiStripNoises> noiseHandle;
  iSetup.get<SiStripNoisesRcd>().get(noiseHandle);

  edm::ESHandle<SiStripGain> gainHandle;
  iSetup.get<SiStripGainRcd>().get(gainHandle);

  edm::ESHandle<SiStripQuality> qualityHandle;
  iSetup.get<SiStripQualityRcd>().get(qualityLabel_,qualityHandle);

  //iSetup.get<SiStripDetCablingRcd>().get(SiStripDetCabling_);

  // get collection of DetSetVector of clusters from Event
  edm::Handle< edmNew::DetSetVector<SiStripCluster> > cluster_detsetvektor;
  iEvent.getByLabel(clusterProducerStrip_, cluster_detsetvektor);

  if (!cluster_detsetvektor.isValid()) return;
  
  //const edmNew::DetSetVector<SiStripCluster> * StrC= cluster_detsetvektor.product();
  //NStripClusters= StrC->data().size(); 


  iSetup.get<SiStripDetCablingRcd>().get(SiStripDetCabling_);  
  // get list of active detectors from SiStripDetCabling 
  std::vector<uint32_t> activeDets;
  SiStripDetCabling_->addActiveDetectorsRawIds(activeDets);

  ModMEs mod_single;

  for(std::vector<uint32_t>::iterator detid_iterator = activeDets.begin(); detid_iterator!=activeDets.end(); detid_iterator++){
    uint32_t detid = (*detid_iterator);
    // remove any eventual zero elements - there should be none, but just in case
    if(detid == 0) {
      activeDets.erase(detid_iterator);
      continue;
    }

    SiStripNewFolderOrganizer* folder = folder_organizer.SearchInSubDirectories( detid );
    if( folder != NULL ){
      folder->NavigateHere() ;
      mod_single.NumberOfClusters = folder->GetME( detid ,  "NumberOfClusters" ) ;
      mod_single.ClusterPosition = folder->GetME( detid ,  "ClusterPosition");
      mod_single.ClusterDigiPosition = folder->GetME( detid ,  "ClusterDigiPosition");
      mod_single.ClusterWidth = folder->GetME( detid ,  "ClusterWidth" ); 
      mod_single.ClusterCharge = folder->GetME( detid ,  "ClusterCharge"); 
      mod_single.ClusterNoise = folder->GetME( detid ,  "ClusterNoise"); 
      mod_single.ClusterSignalOverNoise = folder->GetME( detid ,  "ClusterSignalOverNoise");
      mod_single.ClusterSignalOverNoiseVsPos = folder->GetME( detid ,  "ClusterSignalOverNoiseVsPos");
      mod_single.ModuleLocalOccupancy = folder->GetME( detid ,  "ModuleLocalOccupancy"); 
      mod_single.NrOfClusterizedStrips = folder->GetME( detid ,  "NrOfClusterizedStrips"); 

      edmNew::DetSetVector<SiStripCluster>::const_iterator isearch = cluster_detsetvektor->find(detid); // search  clusters of detid
      
      if(isearch==cluster_detsetvektor->end())
	continue; // no clusters for this detid => jump to next step of loop
      
      //cluster_detset is a structure, cluster_detset.data is a std::vector<SiStripCluster>, cluster_detset.id is uint32_t
      edmNew::DetSet<SiStripCluster> cluster_detset = (*cluster_detsetvektor)[detid]; // the statement above makes sure there exists an element with 'detid'
      
      if(mod_single.NumberOfClusters != NULL){ // nr. of clusters per module
	(mod_single.NumberOfClusters)->Fill(static_cast<float>(cluster_detset.size()));
      }
      

      short total_clusterized_strips = 0;

      SiStripNoises::Range detNoiseRange = noiseHandle->getRange(detid);
      SiStripApvGain::Range detGainRange = gainHandle->getRange(detid); 
      SiStripQuality::Range qualityRange = qualityHandle->getRange(detid);
      

      for(edmNew::DetSet<SiStripCluster>::const_iterator clusterIter = cluster_detset.begin(); clusterIter!= cluster_detset.end(); clusterIter++){

	const std::vector<uint8_t>& ampls = clusterIter->amplitudes();
	// cluster position
	float cluster_position = clusterIter->barycenter();
	// start defined as nr. of first strip beloning to the cluster
	short cluster_start    = clusterIter->firstStrip();
	// width defined as nr. of strips that belong to cluster
	short cluster_width    = ampls.size(); 
	// add nr of strips of this cluster to total nr. of clusterized strips
	total_clusterized_strips = total_clusterized_strips + cluster_width; 
	
	// cluster signal and noise from the amplitudes
	float cluster_signal = 0.0;
	float cluster_noise  = 0.0;
	int nrnonzeroamplitudes = 0;
	float noise2 = 0.0;
	float noise  = 0.0;
	for(uint iamp=0; iamp<ampls.size(); iamp++){
	  if(ampls[iamp]>0){ // nonzero amplitude
	    cluster_signal += ampls[iamp];
	    if(!qualityHandle->IsStripBad(qualityRange, clusterIter->firstStrip()+iamp)){
	      noise = noiseHandle->getNoise(clusterIter->firstStrip()+iamp,detNoiseRange)/gainHandle->getStripGain(clusterIter->firstStrip()+iamp, detGainRange);
	    }
	    noise2 += noise*noise;
	    nrnonzeroamplitudes++;
	  }
	}// End loop over cluster amplitude
	
	if (nrnonzeroamplitudes > 0) cluster_noise = sqrt(noise2/nrnonzeroamplitudes);
	
	if( applyClusterQuality_ &&
	    (cluster_signal/cluster_noise < sToNLowerLimit_ ||
	     cluster_signal/cluster_noise > sToNUpperLimit_ ||
	     cluster_width < widthLowerLimit_ ||
	     cluster_width > widthUpperLimit_) ) continue;  
	
	ClusterProperties cluster_properties;
	cluster_properties.charge    = cluster_signal;
	cluster_properties.position  = cluster_position;
	cluster_properties.start     = cluster_start;
	cluster_properties.width     = cluster_width;
	cluster_properties.noise     = cluster_noise;
	
	// Fill Module Level MEs
	fillModuleMEs(mod_single, cluster_properties);
      
      }
    }

  }
  
}
      

// -- EndJob
//
void SiStripNewMonitorCluster::endJob(void){
  bool outputMEsInRootFile = conf_.getParameter<bool>("OutputMEsInRootFile");
  std::string outputFileName = conf_.getParameter<std::string>("OutputFileName");
 
  // save histos in a file
  if(outputMEsInRootFile) dqmStore_->save(outputFileName);
}
//
// -- Create Module Level MEs
//
void SiStripNewMonitorCluster::createModuleMEs(ModMEs& mod_single, uint32_t detid) {

  // use SistripHistoId for producing histogram id (and title)
  SiStripHistoId hidmanager;
  std::string hid;
  
  //nr. of clusters per module
  hid = hidmanager.createHistoId("NumberOfClusters","det",detid);
  mod_single.NumberOfClusters = bookME1D("TH1nClusters", hid.c_str());
  dqmStore_->tag(mod_single.NumberOfClusters, detid);
  mod_single.NumberOfClusters->setAxisTitle("number of clusters in one detector module");
  mod_single.NumberOfClusters->getTH1()->StatOverflows(kTRUE);  // over/underflows in Mean calculation

  //ClusterPosition
  short total_nr_strips = SiStripDetCabling_->nApvPairs(detid) * 2 * 128; // get correct # of avp pairs    
  hid = hidmanager.createHistoId("ClusterPosition","det",detid);
  mod_single.ClusterPosition = dqmStore_->book1D(hid, hid, total_nr_strips, 0.5, total_nr_strips+0.5);
  dqmStore_->tag(mod_single.ClusterPosition, detid);
  mod_single.ClusterPosition->setAxisTitle("cluster position [strip number +0.5]");

  //ClusterDigiPosition
  //    short total_nr_strips = SiStripDetCabling_->nApvPairs(detid) * 2 * 128; // get correct # of avp pairs    
    hid = hidmanager.createHistoId("ClusterDigiPosition","det",detid);
    mod_single.ClusterDigiPosition = dqmStore_->book1D(hid, hid, total_nr_strips, 0.5, total_nr_strips+0.5);
    dqmStore_->tag(mod_single.ClusterDigiPosition, detid);
    mod_single.ClusterDigiPosition->setAxisTitle("digi in cluster position [strip number +0.5]");
  
  //ClusterWidth
    hid = hidmanager.createHistoId("ClusterWidth","det",detid);
    mod_single.ClusterWidth = bookME1D("TH1ClusterWidth", hid.c_str());
    dqmStore_->tag(mod_single.ClusterWidth, detid);
    mod_single.ClusterWidth->setAxisTitle("cluster width [nr strips]");
  
  //ClusterCharge
    hid = hidmanager.createHistoId("ClusterCharge","det",detid);
    mod_single.ClusterCharge = bookME1D("TH1ClusterCharge", hid.c_str());
    dqmStore_->tag(mod_single.ClusterCharge, detid);
    mod_single.ClusterCharge->setAxisTitle("cluster charge [ADC]");

  
  //ClusterNoise
    hid = hidmanager.createHistoId("ClusterNoise","det",detid);
    mod_single.ClusterNoise = bookME1D("TH1ClusterNoise", hid.c_str());
    dqmStore_->tag(mod_single.ClusterNoise, detid);
    mod_single.ClusterNoise->setAxisTitle("cluster noise");

  
  //ClusterSignalOverNoise
    hid = hidmanager.createHistoId("ClusterSignalOverNoise","det",detid);
    mod_single.ClusterSignalOverNoise = bookME1D("TH1ClusterStoN", hid.c_str());
    dqmStore_->tag(mod_single.ClusterSignalOverNoise, detid);
    mod_single.ClusterSignalOverNoise->setAxisTitle("ratio of signal to noise for each cluster");


  //ClusterSignalOverNoiseVsPos
    hid = hidmanager.createHistoId("ClusterSignalOverNoiseVsPos","det",detid);
    Parameters =  conf_.getParameter<edm::ParameterSet>("TH1ClusterStoNVsPos");
    mod_single.ClusterSignalOverNoiseVsPos= dqmStore_->bookProfile(hid.c_str(),hid.c_str(),
								   Parameters.getParameter<int32_t>("Nbinx"),
								   Parameters.getParameter<double>("xmin"),
								   Parameters.getParameter<double>("xmax"),
								   Parameters.getParameter<int32_t>("Nbiny"),
								   Parameters.getParameter<double>("ymin"),
								   Parameters.getParameter<double>("ymax")
								   );
    dqmStore_->tag(mod_single.ClusterSignalOverNoiseVsPos, detid);
    mod_single.ClusterSignalOverNoiseVsPos->setAxisTitle("pos");

  
  //ModuleLocalOccupancy
    hid = hidmanager.createHistoId("ClusterLocalOccupancy","det",detid);
    mod_single.ModuleLocalOccupancy = bookME1D("TH1ModuleLocalOccupancy", hid.c_str());
    dqmStore_->tag(mod_single.ModuleLocalOccupancy, detid);
    mod_single.ModuleLocalOccupancy->setAxisTitle("module local occupancy [% of clusterized strips]");

  
  //NrOfClusterizedStrips
    hid = hidmanager.createHistoId("NrOfClusterizedStrips","det",detid);
    mod_single.NrOfClusterizedStrips = bookME1D("TH1NrOfClusterizedStrips", hid.c_str());
    dqmStore_->tag(mod_single.NrOfClusterizedStrips, detid);
    mod_single.NrOfClusterizedStrips->setAxisTitle("number of clusterized strips");

}  

//
// -- Fill Module Level Histograms
//
void SiStripNewMonitorCluster::fillModuleMEs(ModMEs& mod_mes, ClusterProperties& cluster) {
  
  if((mod_mes.ClusterPosition)) // position of cluster
    (mod_mes.ClusterPosition)->Fill(cluster.position);
  
  // position of digis in cluster
  if((mod_mes.ClusterDigiPosition)) {
    for(int ipos=cluster.start+1; ipos<=cluster.start+cluster.width; ipos++){
      (mod_mes.ClusterDigiPosition)->Fill(ipos);
    }
  }

  if((mod_mes.ClusterWidth)) // width of cluster
    (mod_mes.ClusterWidth)->Fill(static_cast<float>(cluster.width));
 
  if((mod_mes.ClusterSignalOverNoise)) {// SignalToNoise
    if (cluster.noise > 0) 
      (mod_mes.ClusterSignalOverNoise)->Fill(cluster.charge/cluster.noise);
  }

  if((mod_mes.ClusterSignalOverNoiseVsPos)) {// SignalToNoise
    if (cluster.noise > 0) 
      (mod_mes.ClusterSignalOverNoiseVsPos)->Fill(cluster.position,cluster.charge/cluster.noise);
  }

  if((mod_mes.ClusterNoise))  // Noise
    (mod_mes.ClusterNoise)->Fill(cluster.noise);

  if((mod_mes.ClusterCharge)) // charge of cluster
    (mod_mes.ClusterCharge)->Fill(cluster.charge);
  
} 
//------------------------------------------------------------------------------------------
MonitorElement* SiStripNewMonitorCluster::bookMETrend(const char* ParameterSetLabel, const char* HistoName)
{
  Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
  edm::ParameterSet ParametersTrend =  conf_.getParameter<edm::ParameterSet>("Trending");
  MonitorElement* me = dqmStore_->bookProfile(HistoName,HistoName,
					      ParametersTrend.getParameter<int32_t>("Nbins"),
					      // 					      0,
					      ParametersTrend.getParameter<double>("xmin"),
					      ParametersTrend.getParameter<double>("xmax"),
					      // 					      ParametersTrend.getParameter<int32_t>("Nbins"),
					      100, //that parameter should not be there !?
					      ParametersTrend.getParameter<double>("ymin"),
					      ParametersTrend.getParameter<double>("ymax"),
					      "" );
  if(!me) return me;
  me->setAxisTitle("Event Time in Seconds",1);
  if (me->kind() == MonitorElement::DQM_KIND_TPROFILE) me->getTH1()->SetBit(TH1::kCanRebin);
  return me;
}

//------------------------------------------------------------------------------------------
MonitorElement* SiStripNewMonitorCluster::bookME1D(const char* ParameterSetLabel, const char* HistoName)
{
  Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
  return dqmStore_->book1D(HistoName,HistoName,
			   Parameters.getParameter<int32_t>("Nbinx"),
			   Parameters.getParameter<double>("xmin"),
			   Parameters.getParameter<double>("xmax")
			   );
}

// int SiStripNewMonitorCluster::FindRegion(int nstrip,int npix){
  
//   double kplus= k0*(1+dk0/100);
//   double kminus=k0*(1-dk0/100);
//   int region=0;
  
//   if (nstrip!=0 && npix >= (nstrip*kminus-q0) && npix <=(nstrip*kplus+q0)) region=1; 
//   else if (nstrip!=0 && npix < (nstrip*kminus-q0) &&  nstrip <= maxClus) region=2;
//   else if (nstrip!=0 && npix < (nstrip*kminus-q0) &&  nstrip > maxClus) region=3;
//   else if (nstrip!=0 && npix > (nstrip*kplus+q0)) region=4;
//   else if (npix > minPix && nstrip==0) region=5;
//   return region;

// }


    
