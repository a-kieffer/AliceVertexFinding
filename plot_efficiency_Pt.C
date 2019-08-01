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
#include "ITStracking/Tracklet.h"
#include "ITStracking/Constants.h"

#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

void FillPtBins(std::vector<o2::its::Tracklet> tracklets, int layer, std::map<double,int> * PtBins,  std::map<o2::MCCompLabel,o2::its::label>  LabelMap, 
                o2::its::Vertexer* vertexer, o2::its::ROframe frame);

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

  std::vector<float> BinsLimits; 
  for(int i=0; i<80; i++){
      BinsLimits.push_back(0.05*i);
  }
  for(int i=0; i<12; i++){
      BinsLimits.push_back(4+0.5*i);
  }
  for(int i=0; i<5; i++){
      BinsLimits.push_back(10+2*i);
  }

  int NbBins = BinsLimits.size()-1;
  std::cout<< "Number of bins : "<< NbBins << std::endl;

  std::map<double,int> PtBins01 ;
  std::map<double,int> RealPtBins01 ;
  std::map<double,int> PtBins12 ;
  std::map<double,int> RealPtBins12 ;

  for(auto & binlimit : BinsLimits){
    PtBins01.insert(std::pair<double,int>(binlimit,0));
    RealPtBins01.insert(std::pair<double,int>(binlimit,0));
    PtBins12.insert(std::pair<double,int>(binlimit,0));
    RealPtBins12.insert(std::pair<double,int>(binlimit,0));
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
    std::cout << "ROframe: " << iROfCount << std::endl;
     
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

    auto tracklets01 = vertexer.getTracklets01();
    auto tracklets12 = vertexer.getTracklets12();

    FillPtBins(tracklets01, 0, &PtBins01,   LabelMap, &vertexer,  frame);
    FillPtBins(tracklets12, 1, &PtBins12,   LabelMap, &vertexer,  frame);

    vertexer.initialiseVertexer(&frame);
    vertexer.findTrivialMCTracklets(LabelMap);
    tracklets01 = vertexer.getTracklets01();
    tracklets12 = vertexer.getTracklets12();

    FillPtBins(tracklets01, 0, &RealPtBins01,   LabelMap, &vertexer,  frame);
    FillPtBins(tracklets12, 1, &RealPtBins12,   LabelMap, &vertexer,  frame);

    
  }


  std::cout<< "Pt bins 01 content : \n";
    for( auto const& [key, val] : PtBins01 ){
      std::cout << key    // string (key)
              << ':'  << val        // string's value
              << std::endl ;
      }

/*
    std::cout<< "Real Pt bins 12 content : \n";
    for( auto const& [key, val] : RealPtBins12 ){
      std::cout << key    // string (key)
              << ':'  << val        // string's value
              << std::endl ;
      }

 */

  //plotting histograms 
  /* 
  TCanvas *c1= new TCanvas ("c1", "Histograms", 1600, 900);
  c1->SetLogy();
  c1->SetLogy();
  TH1F *Hist01 = new TH1F("hist01", "Reconstructed Tracklets 01 vs Pt", NbBins, BinsLimits.data());
  TH1F *HistReal01 = new TH1F("realhist01", "Real Tracklets 01 vs Pt", NbBins, BinsLimits.data());

  for( auto const& [key, val] : PtBins01 ){
    for(int i =0; i< val;i++){
      //std::cout << "key :" <<key << " i :" << i << std::endl;
      Hist01->Fill(key);
    }
  }

  for( auto const& [key, val] : RealPtBins01 ){
    for(int i =0; i< val;i++){
      //std::cout << "key :" <<key << " i :" << i << std::endl;
      HistReal01->Fill(key);
    }
  }

  std::cout << "Entries in the histogram : "<< Hist01->GetEntries() << std::endl;

  Hist01->Draw();
  HistReal01->SetLineColor(kOrange + 8);
  HistReal01->Draw("same");

*/

  //plotting the graphs
  std::vector<float> Efficiency01;

  for( auto const& [key, val] : PtBins01 ){
    Efficiency01.push_back(((double)val)/((double)RealPtBins01[key]));
  }
  
  std::vector<float> Efficiency12;

  for( auto const& [key, val] : PtBins12 ){
    Efficiency12.push_back(((double)val)/((double)RealPtBins12[key]));
  }


  

  BinsLimits.pop_back();
  Efficiency01.pop_back();
  Efficiency12.pop_back();
  TCanvas *c= new TCanvas ("c", "Efficiency graphs", 1600, 900);
  c->Divide(2,1);
  c->cd(1);
  c->cd(1)->SetLogx();
  TGraph * graphEff01= new TGraph(NbBins, BinsLimits.data(),Efficiency01.data());
  graphEff01->SetMarkerStyle(8);
  graphEff01->GetXaxis()->SetTitle("Pt ( GeV) ");
  graphEff01->GetYaxis()->SetTitle("Reco/Total");
  graphEff01->SetTitle("Tracklets 01 reconstruction efficiency vs Pt");
  graphEff01->Draw("APC");

  c->cd(2);
  c->cd(2)->SetLogx();
  TGraph * graphEff12= new TGraph(NbBins, BinsLimits.data(),Efficiency12.data());
  graphEff12->SetMarkerStyle(8);
  graphEff12->GetXaxis()->SetTitle("Pt ( GeV) ");
  graphEff12->GetYaxis()->SetTitle("Reco/Total");
  graphEff12->SetTitle("Tracklets 12 reconstruction efficiency vs Pt");
  graphEff12->Draw("APC");

  std::cout << "To plot : \n";
  for(unsigned int i =0; i< Efficiency01.size(); i++) {
    std::cout<< "Bin : "<< BinsLimits[i] << " efficiency : "<< Efficiency01[i]<< std::endl;
  }
}


 
void FillPtBins(std::vector<o2::its::Tracklet> tracklets, int layer, std::map<double,int> * PtBins,  std::map<o2::MCCompLabel,o2::its::label>  LabelMap, 
                o2::its::Vertexer * vertexer, o2::its::ROframe frame){

  for(auto & tracklet : tracklets){
    auto Cluster = vertexer->getClusters()[layer][tracklet.firstClusterIndex] ;
    const auto MCLabel = frame.getClusterLabels(layer, Cluster.clusterId);
    if(LabelMap.find(MCLabel)!=LabelMap.end()){ //the cluster is in the LabelMap : it was reconstructed, we have to get the Pt
      //std::cout << "Found a Cluster \n";
      double Pt = LabelMap[MCLabel].Pt;
      //std::cout<< "Pt of the Cluster : "<< Pt << std::endl;
      for ( std::map<double,int>::iterator it = PtBins->begin(); it != PtBins->end(); it++ ){
        //std::cout << "Bin : " << it->first << std::endl; //does not work
        if(it->first > Pt){ //the bin beginning is greater than the Pt : it belongs to the previous bin
          it--;
          PtBins->at((it)->first)++;
          break;
          }
        }
      }
    }
}



#endif