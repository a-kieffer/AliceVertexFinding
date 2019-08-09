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

#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/ROframe.h"

#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"

void event_map(
    const int startRof = 0,
    const int nRof = -1,
    const unsigned char muted = true,
    const std::string path = "./",
    const std::string inputClustersITS = "o2clus_its.root",
    const std::string simfilename = "o2sim.root",
    const std::string paramfilename = "O2geometry.root")
{
  TString outputfile;
  outputfile.Form("label2Track%d.root", startRof);

  // Setup Runtime DB
  TFile paramFile((path + paramfilename).data());
  paramFile.Get("FAIRGeom");
  auto gman = o2::its::GeometryTGeo::Instance();
  gman->fillMatrixCache(
      o2::utils::bit2Mask(o2::TransformType::T2L, o2::TransformType::T2GRot,
                          o2::TransformType::L2G)); // request cached transforms

  // ROframe data
  TChain itsClustersROF("ITSClustersROF");
  TChain itsClusters("o2sim");
  itsClusters.AddFile((path + inputClustersITS).data());
  itsClustersROF.AddFile((path + inputClustersITS).data());
  std::vector<o2::itsmft::Cluster> *clusters = nullptr;
  std::vector<o2::itsmft::ROFRecord> *rofs = nullptr; // only one rofrecord
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *labels = nullptr;
  itsClusters.SetBranchAddress("ITSCluster", &clusters);
  itsClustersROF.SetBranchAddress("ITSClustersROF", &rofs);
  itsClusters.SetBranchAddress("ITSClusterMCTruth", &labels);
  itsClustersROF.GetEntry(0);
  itsClusters.GetEntry(0);

  // Simulation data
  TChain mcHeaderTree("o2sim");
  mcHeaderTree.AddFile((path + simfilename).data());
  std::vector<o2::MCTrack> *tracks = nullptr;
  mcHeaderTree.SetBranchAddress("MCTrack", &tracks);

  // Output data
  TTree outTree("Labels2Tracks", "labels2tracks");
  std::vector<o2::MCCompLabel> *labelsToSave = new std::vector<o2::MCCompLabel>;
  std::vector<o2::MCTrack> *tracksToSave = new std::vector<o2::MCTrack>;
  outTree.Branch("MCLabels", &labelsToSave);
  outTree.Branch("Tracks", &tracksToSave);

  std::vector <int> ROFNumbers ;

  std::vector <int>  EventNumbers ;


  o2::its::GeometryTGeo *geom = o2::its::GeometryTGeo::Instance();
  geom->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2GRot));
  int stopAt = nRof == -1 ? (*rofs).size() : static_cast<int>(std::min(static_cast<int>((*rofs).size()), startRof + nRof));

  //std::cout << "Rofs size : " << (*rofs).size() << std::endl;

  for (auto iRof{startRof}; iRof < stopAt; ++iRof)
  {
    auto &rofRec = (*rofs)[iRof];
    int clusterId{0};
    int first = static_cast<int>(rofRec.getROFEntry().getIndex());
    int number = static_cast<int>(rofRec.getNROFEntries());
    std::vector<o2::MCCompLabel> roLabels;
    std::vector<o2::MCTrack> roTracks;
    if (!muted)
      std::cerr << "ROframe: " << iRof << std::endl;
    auto clusters_in_frame = gsl::make_span(&(*clusters)[first], number);

    int eventID = -1;

    for (auto &c : clusters_in_frame){
      if (geom->getLayer(c.getSensorID()) == 0) // Select clusters on Layer 0
      {
        const auto &lbl = *(labels->getLabels(first + clusterId).begin());
        
        if (lbl.isValid() && lbl.getSourceID() == 0)
        {
          if (mcHeaderTree.GetReadEntry() != lbl.getEventID())
          {
            // to avoid repetitive reading of the same event
            mcHeaderTree.GetEntry(lbl.getEventID());
            eventID = lbl.getEventID();
          }

          
        }
      }
      clusterId++;
    }
    ROFNumbers.push_back(iRof);
    EventNumbers.push_back(eventID);
  }
  for (int i =0; i<ROFNumbers.size(); i++){
      if(EventNumbers[i-startRof] != -1){
           std::cout << " ROFRame "<< i << " : Event " << EventNumbers[i-startRof] << std::endl;
      }
     
  }
 // std::cout<< "number of events :" << EventNumbers.size() << std::endl;
 
}
#endif