// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#define GSL_THROW_ON_CONTRACT_VIOLATION
#include <gsl/gsl>

#include "DataFormatsITSMFT/Cluster.h"
#include "ITStracking/Cluster.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPObject.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCTrack.h"

#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/ROframe.h"

#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "ITStracking/IOUtils.h"
#include "ITStracking/Vertexer.h"
#include "ITStracking/Constants.h"

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

void plot_efficiency_Pt(const int startRof = 0,
                        const int nRof = -1,
                        const unsigned char muted = true,
                        const std::string path = "./",
                        const std::string inputClustersITS = "o2clus_its.root",
                        const std::string simfilename = "o2sim.root",
                        const std::string paramfilename = "O2geometry.root"){



  TChain itsClusters("o2sim");
  itsClusters.AddFile((path + inputClustersITS).data());

  // Setup Runtime DB
  TFile paramFile((path + paramfilename).data());
  paramFile.Get("FAIRGeom");
  auto gman = o2::its::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2L, o2::TransformType::T2GRot,
                                            o2::TransformType::L2G)); // request cached transforms

  // Get event header
  TChain mcHeaderTree("o2sim");
  mcHeaderTree.AddFile((path + simfilename).data());
  o2::dataformats::MCEventHeader* mcHeader = nullptr;
  if (!mcHeaderTree.GetBranch("MCEventHeader.")) {
    LOG(FATAL) << "Did not find MC event header in the input header file." << FairLogger::endl;
  }
  mcHeaderTree.SetBranchAddress("MCEventHeader.", &mcHeader);

  // get clusters
  std::vector<o2::itsmft::Cluster>* clusters = nullptr;
  itsClusters.SetBranchAddress("ITSCluster", &clusters);

  TChain itsClustersROF("ITSClustersROF");
  itsClustersROF.AddFile((path + inputClustersITS).data());

  if (!itsClustersROF.GetBranch("ITSClustersROF")) {
    LOG(FATAL) << "Did not find ITS clusters branch ITSClustersROF in the input tree" << FairLogger::endl;
  }
  std::vector<o2::itsmft::ROFRecord>* rofs = nullptr;
  itsClustersROF.SetBranchAddress("ITSClustersROF", &rofs);
  itsClustersROF.GetEntry(0);

  //get labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels = nullptr;
  itsClusters.SetBranchAddress("ITSClusterMCTruth", &labels);

  //getting the clusters ID etc
  std::string file = "label2Track0.root";
  TFile * labelFile = new TFile(file.data());
  TTree * labelsData ;
  labelFile->GetObject("Labels2Tracks", labelsData);
  
  std::vector<o2::MCCompLabel>  * MCLabels = nullptr;
	TBranch * MCLabelsBranch = labelsData->GetBranch("MCLabels");
	MCLabelsBranch->SetAddress(& MCLabels);
 
  std::vector<o2::MCTrack>   * Tracks = nullptr;
	TBranch * TracksBranch = labelsData->GetBranch("Tracks");
	TracksBranch->SetAddress(& Tracks);

  std::vector<double> BinsLimits; 
  for(int i=0; i<20; i++){
      BinsLimits.push_back(0.05*i);
  }
  //BinsLimits.push_back(0.25);
  BinsLimits.push_back(1);
  BinsLimits.push_back(2);
  BinsLimits.push_back(5);
  BinsLimits.push_back(10);
  BinsLimits.push_back(20);

  int NbBins = BinsLimits.size()-1;
  std::cout<< "Number of bins : "<< NbBins << std::endl;

  std::map<double,int> PtBins ;

  for(auto & binlimit : BinsLimits){
    PtBins.insert(std::pair<double,int>(binlimit,0));
  }

  std::uint32_t roFrame = 0;
  o2::its::ROframe frame(-123);
  o2::its::VertexerTraits* traits = nullptr;
  traits = o2::its::createVertexerTraits();
 
  o2::its::Vertexer vertexer(traits);

  const int stopAt = (nRof == -1) ? rofs->size() : startRof + nRof;

  itsClusters.GetEntry(0);
  mcHeaderTree.GetEntry(0);

  for (size_t iROfCount{ static_cast<size_t>(startRof) }; iROfCount < static_cast<size_t>(stopAt); ++iROfCount) {

    
    auto& rof = (*rofs)[iROfCount];
    //std::cout << "  ROframe: " << iROfCount << std::endl;
     
    int nclUsed = o2::its::ioutils::loadROFrameData(rof, frame, clusters, labels);

    std::array<float, 3> total{ 0.f, 0.f, 0.f };
    o2::its::ROframe* eventptr = &frame;
    
    labelsData->GetEntry(iROfCount);

    std::map<o2::MCCompLabel,o2::its::label>  LabelMap ;
    

    for(unsigned int i=0; i< Tracks->size(); i++){
      o2::MCCompLabel tmpMC = (*MCLabels)[i];
      o2::MCTrack tmpTrack = (*Tracks)[i];
      o2::its::label tmpLabel = {tmpMC.getTrackID(),  tmpTrack.getMotherTrackId(), tmpMC.getEventID(), static_cast<float>( tmpTrack.GetPt())};
      LabelMap.insert( std::pair <o2::MCCompLabel,o2::its::label> (tmpMC, tmpLabel));
    }

    vertexer.initialiseVertexer(eventptr);
    vertexer.findTracklets(true);
    
    int count =0;

    auto tracklets01 = vertexer.getTracklets01();
    for(auto & tracklet : tracklets01){
      auto Cluster0 = vertexer.getClusters()[0][tracklet.firstClusterIndex] ;
      const auto MCLabel0 = frame.getClusterLabels(0, Cluster0.clusterId);
      if(LabelMap.find(MCLabel0)!=LabelMap.end()){ //the cluster is in the LabelMap : it was reconstructed, we have to get the Pt
        //std::cout << "Found a Cluster \n";
        double Pt = LabelMap[MCLabel0].Pt;
        std::cout<< "Pt of the Cluster : "<< Pt << std::endl;
        //std::map<double,int>::iterator it;
        for ( std::map<double,int>::iterator it = PtBins.begin(); it != PtBins.end(); it++ ){
          //std::cout << "Bin : " << it->first << std::endl; //does not work
          if(it->first > Pt){ //the bin beginning is greater than the Pt : it belongs to the previous bin
            it--;
            PtBins[(it)->first]++;
            break;
          }
        }

      }
      //std::cout << "Cluster 0 track ID : " << Label0.getTrackID()<< std::endl;
    }

    //for (auto & label : LabelMap){
    //  std::cout << "track ID in LabelMap: " << label.TrackId << std::endl;
    //}

    //vertexer.findTrivialMCTracklets(LabelMap);

    std::cout<< "Pt bins content : \n";
    for( auto const& [key, val] : PtBins ){
    std::cout << key         // string (key)
              << ':'  
              << val        // string's value
              << std::endl ;
}

    vertexer.findVertices();
    vertexer.processLines();
    std::vector<Vertex> vertITS = vertexer.exportVertices();

    if ( vertITS.size() > 0 ) {
     std::cout<<"Found Vertex : x : "<<vertITS[0].getX()<<" y : "<<vertITS[0].getY()<<" z "<<vertITS[0].getZ()<<std::endl;
    }
    
    
  }


}
#endif