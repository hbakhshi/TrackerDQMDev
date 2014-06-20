#include "../interface/SiStripNewFolderOrganizer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

SiStripNewFolderOrganizer::SiStripNewFolderOrganizer( edm::ParameterSet parent_dir , const SiStripNewFolderOrganizer* TheParent )
  :detIdSelector( parent_dir.getParameter< std::vector<std::string> >("DetIds")){
  this->Name = parent_dir.getParameter< string >("Name");
  //  this->DetIdRanges = 
  this->Parent = TheParent ;

  SubDirInfo = parent_dir.getParameterSetVector("SubDirs");
  SubDirs.reserve(SubDirInfo.size());

  int subdir_id = 0;
  for( edm::VParameterSet::const_iterator subdir = SubDirInfo.begin() ; subdir != SubDirInfo.end() ; subdir++){
    SiStripNewFolderOrganizer* subfolder = new SiStripNewFolderOrganizer( *subdir , this ); //(&SubDirs[subdir_id])
    SubDirs.push_back( subfolder );
    subdir_id ++;
  }

  dbe_  = edm::Service<DQMStore>().operator->();
  dbe_->setCurrentFolder( this->GetFullName() );
}

void SiStripNewFolderOrganizer::NavigateHere(){
  dbe_->setCurrentFolder( this->GetFullName() );
}

bool SiStripNewFolderOrganizer::Contains( const unsigned int& detid) const{
  return this->detIdSelector( detid );
}

SiStripNewFolderOrganizer* SiStripNewFolderOrganizer::SearchInSubDirectories( const unsigned int& detid ){
  if( this->Contains( detid ) )
    return this;
  
  for( std::vector< SiStripNewFolderOrganizer* >::const_iterator subs = SubDirs.begin() ; subs != SubDirs.end() ; subs++ ){
    SiStripNewFolderOrganizer* ret = (*subs)->SearchInSubDirectories( detid );
    if( ret != NULL )
      return ret;
  }

  return  NULL;
}


void SiStripNewFolderOrganizer::AddME( string name , MonitorElement* me ){
  SummaryMEs[ name ] = me ; 
}

void SiStripNewFolderOrganizer::AddME( const unsigned int& detid , string name , MonitorElement* me){
  map< string , MonitorElement* >* all_module_mes = &( ModuleMEs[ detid ] ) ;
  (*all_module_mes)[ name ] = me; 
}


MonitorElement* SiStripNewFolderOrganizer::GetME( string name ){
  map<string, MonitorElement*>::const_iterator value = this->SummaryMEs.find( name );

  if( value != SummaryMEs.end() )
    return value->second ;
  else
    return NULL;
}

MonitorElement* SiStripNewFolderOrganizer::GetME( const unsigned int& detid , string name){

  map< string , MonitorElement* >* all_module_mes = &( ModuleMEs[ detid ] ) ;
  map<string, MonitorElement*>::const_iterator value = all_module_mes->find( name );

  if( value != all_module_mes->end() )
    return value->second ;
  else
    return NULL;

}

string SiStripNewFolderOrganizer::GetTabs() const{
  string ret;
  if(Parent == NULL)
    ret = "";
  else
    ret = (Parent)->GetTabs() + "\t";

  return ret;
}

string SiStripNewFolderOrganizer::GetFullName() const{
  string ret;
  if(Parent == NULL)
    ret = this->Name ;
  else
    ret = (Parent)->GetFullName() + "/" + this->Name ;

  return ret;
}

void SiStripNewFolderOrganizer::Print() const{
  cout << GetFullName() << endl;
  for( std::vector< SiStripNewFolderOrganizer* >::const_iterator subs = SubDirs.begin() ; subs != SubDirs.end() ; subs++ )
    (*subs)->Print();
}
