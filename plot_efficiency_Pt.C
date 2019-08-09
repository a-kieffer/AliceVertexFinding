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
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
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


#include "ITStracking/IOUtils.h"
#include "ITStracking/Vertexer.h"
#include "ITStracking/Tracklet.h"
#include "ITStracking/Constants.h"



using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

void FillPtBins(std::vector<o2::its::Tracklet> tracklets, int layer, std::map<double,int> * PtBins, std::vector<double> * RawPt, 
                std::map<o2::MCCompLabel,o2::its::ClusterMCLabelInfo>  LabelMap, o2::its::Vertexer* vertexer, o2::its::ROframe frame);

void plot_efficiency_Pt(const int startRof = 0,
                        const int nRof = -1,
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

  //output file to save the plots

  TFile * outputfile = new TFile("efficiency_pt_plots.root", "recreate");

  //creating all the needed vectors 

  std::vector<float> BinsLimits; //vector that contains the limit values for all the pT bins for the histogram

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


  // Pt bins to plot the graphs
  std::map<double,int> PtBins01 ; //number of reconstructed tracklets corresponding to one pT bin
  // 15 tracklets in bin 0.2 and x in bin 0.25 means that there are 15 tracklets with a bin between 0.2 and 0.25
  std::map<double,int> RealPtBins01 ; //number of real tracklers
  std::map<double,int> PtBins12 ;
  std::map<double,int> RealPtBins12 ;
  // all the vectors corresponding to lines are commented out because I did not find a way to get the lines
  //std::map<double,int> PtBinsLines ;
  //std::map<double,int> RealPtBinsLines ;

  //set the number of tracklets for each bin to 0
  for(auto & binlimit : BinsLimits){
    PtBins01.insert(std::pair<double,int>(binlimit,0));
    RealPtBins01.insert(std::pair<double,int>(binlimit,0));
    PtBins12.insert(std::pair<double,int>(binlimit,0));
    RealPtBins12.insert(std::pair<double,int>(binlimit,0));
    //PtBinsLines.insert(std::pair<double,int>(binlimit,0));
    //RealPtBinsLines.insert(std::pair<double,int>(binlimit,0));
  }


  //vectors to fill the histograms : all the values of the Pt, no bins
  std::vector<double> RawPt01;
  std::vector<double> RawPt12;
  //std::vector<double> RawPtLines;
  std::vector<double> RawRealPt01;
  std::vector<double> RawRealPt12;
  //std::vector<double> RawRealPtLines;


  std::uint32_t roFrame = 0;
  o2::its::ROframe frame(-123);
  o2::its::VertexerTraits* traits = nullptr;
  traits = o2::its::createVertexerTraits();
 
  o2::its::Vertexer vertexer(traits);

  o2::its::VertexingParameters parameters;

  //here you can set the value of the cut 
  parameters.phiCut = 0.02f;
  vertexer.setParameters(parameters);

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

    std::map<o2::MCCompLabel,o2::its::ClusterMCLabelInfo>  LabelMap ;
    
    for(unsigned int i=0; i< Tracks->size(); i++){
      o2::MCCompLabel tmpMC = (*MCLabels)[i];
      o2::MCTrack tmpTrack = (*Tracks)[i];
      o2::its::ClusterMCLabelInfo tmpLabel = {tmpMC.getTrackID(),  tmpTrack.getMotherTrackId(), tmpMC.getEventID(), static_cast<float>( tmpTrack.GetPt())};
      LabelMap.insert( std::pair <o2::MCCompLabel,o2::its::ClusterMCLabelInfo> (tmpMC, tmpLabel));
    }

    vertexer.initialiseVertexer(eventptr);
    vertexer.findTracklets(true);

    auto tracklets01 = vertexer.getTracklets01();
    auto tracklets12 = vertexer.getTracklets12();

    std::cout<< "Tracklets 01 : "<< tracklets01.size() << " tracklets 12 : "<< tracklets12.size() << std::endl;
    //auto lines = vertexer.getLines();

    FillPtBins(tracklets01, 0, &PtBins01, &RawPt01,  LabelMap, &vertexer,  frame);
    FillPtBins(tracklets12, 1, &PtBins12, &RawPt12,  LabelMap, &vertexer,  frame);

    vertexer.initialiseVertexer(eventptr);
    vertexer.findTrivialMCTracklets(/*LabelMap */);
    tracklets01 = vertexer.getTracklets01();
    tracklets12 = vertexer.getTracklets12();
    //auto lines = vertexer.getLines();

    FillPtBins(tracklets01, 0, &RealPtBins01, &RawRealPt01,  LabelMap, &vertexer,  frame);
    FillPtBins(tracklets12, 1, &RealPtBins12, &RawRealPt12,  LabelMap, &vertexer,  frame);


    
  }

/*
  std::cout<< "Pt bins 01 content : \n";
    for( auto const& [key, val] : PtBins01 ){
      std::cout << key    // string (key)
              << ':'  << val        // string's value
              << std::endl ;
      }


    std::cout<< "Real Pt bins 12 content : \n";
    for( auto const& [key, val] : RealPtBins12 ){
      std::cout << key    // string (key)
              << ':'  << val        // string's value
              << std::endl ;
      }

 */

  //plotting histograms 

  int nbBinsHist = 400;
   
  TCanvas *cHist= new TCanvas ("CanvasHist", "Reconstruction histograms", 1600, 900);
  cHist->Divide(2,1);

  TH1F *Hist01 = new TH1F("hist01", " Real reconstructed Tracklets 01 vs p_{T};p_{T} (GeV);Reconstructed tracklets",nbBinsHist,0,25 );
  TH1F *HistReal01 = new TH1F("realhist01", "Real reconstructed Tracklets 01 vs p_{T};p_{T} (GeV);Reconstructed tracklets", nbBinsHist,0,25);

  for( auto PtValue : RawPt01 ){
      Hist01->Fill(PtValue);
  }
  for( auto PtValue : RawRealPt01 ){
      HistReal01->Fill(PtValue);
  }

  std::cout << "Entries in the histogram : "<< Hist01->GetEntries() << std::endl;

  cHist->cd(1);
  cHist->cd(1)->SetLogx();
  HistReal01->SetLineColor(kOrange + 8);
  HistReal01->Draw();
  Hist01->Draw("sames");

  

  TH1F *Hist12 = new TH1F("hist12", "Real reconstructed Tracklets 12 vs p_{T}; p_{T} (GeV);Reconstructed tracklet",nbBinsHist,0,25 );
  TH1F *HistReal12 = new TH1F("realhist12", "Real reconstructed Tracklets 12 vs p_{T};p_{T} (GeV);Reconstructed tracklet", nbBinsHist,0,25);

  for( auto PtValue : RawPt12 ){
      Hist12->Fill(PtValue);
  }
  for( auto PtValue : RawRealPt12 ){
      HistReal12->Fill(PtValue);
  }

  std::cout << "Entries in the histogram : "<< Hist01->GetEntries() << std::endl;

  cHist->cd(2);
  cHist->cd(2)->SetLogx();
  HistReal12->SetLineColor(kOrange + 8);
  HistReal12->Draw();
  Hist12->Draw("sames");

/* 
  TCanvas *cHist2= new TCanvas ("CanvasHist", "Selection histogram", 1600, 900);
  TH1F *HistLines = new TH1F("histlines", " Real selected lines vs Pt",100,0,25 );
  TH1F *HistRealLines = new TH1F("reallines", "Real lines vs Pt", 100,0,25);

  for( auto PtValue : RawPtLines ){
      HistLines->Fill(PtValue);
  }
  for( auto PtValue : RawRealPtLines ){
      HistRealLines->Fill(PtValue);
  }

  std::cout << "Entries in the histogram : "<< HistLines->GetEntries() << std::endl;

  HistRealLines->SetLineColor(kOrange + 8);
  HistRealLines->Draw();
  HistLines->Draw("sames");
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
  
  /* 
  std::vector<float> EfficiencyLines;
  for( auto const& [key, val] : PtBinsLines ){
    EfficiencyLines.push_back(((double)val)/((double)RealPtBinsLines[key]));
  }
*/
  //getting rid of the last value because the bin is empty
  BinsLimits.pop_back();
  Efficiency01.pop_back();
  Efficiency12.pop_back();
  //EfficiencyLines.pop_back();
  TCanvas *cGraph= new TCanvas ("cGraph", "Efficiency graphs", 1600, 900);
  cGraph->Divide(2,1);
  cGraph->cd(1);
  cGraph->cd(1)->SetLogx();
  TGraph * graphEff01= new TGraph(NbBins, BinsLimits.data(),Efficiency01.data());
  graphEff01->SetMarkerStyle(8);
  graphEff01->GetXaxis()->SetTitle("p_{T} ( GeV) ");
  graphEff01->GetYaxis()->SetTitle("Real reco/Total");
  graphEff01->SetTitle("Real tracklets 01 reconstruction efficiency vs p_{T}");
  graphEff01->Draw("APL");

  cGraph->cd(2);
  cGraph->cd(2)->SetLogx();
  TGraph * graphEff12= new TGraph(NbBins, BinsLimits.data(),Efficiency12.data());
  graphEff12->SetMarkerStyle(8);
  graphEff12->GetXaxis()->SetTitle("p_{T} ( GeV) ");
  graphEff12->GetYaxis()->SetTitle("Real reco/Total");
  graphEff12->SetTitle("Real tracklets 12 reconstruction efficiency vs p_{T}");
  graphEff12->Draw("APL");

/* 
  TCanvas *cGraph2= new TCanvas ("graphs2", "Efficiency graph for lines", 1600, 900);
  cGraph->cd(1)->SetLogx();
  TGraph * graphEffLines= new TGraph(NbBins, BinsLimits.data(),EfficiencyLines.data());
  graphEffLines->SetMarkerStyle(8);
  graphEffLines->GetXaxis()->SetTitle("Pt ( GeV) ");
  graphEffLines->GetYaxis()->SetTitle("Real selected/Total");
  graphEffLines->SetTitle("Real lines selection efficiency vs Pt");
  graphEffLines->Draw("APC");
*/
/* 
  std::cout << "To plot : \n";
  for(unsigned int i =0; i< Efficiency01.size(); i++) {
    std::cout<< "Bin : "<< BinsLimits[i] << " efficiency : "<< Efficiency01[i]<< std::endl;
  }
  */

  //write the plots into a file
  
  outputfile->WriteTObject(Hist01);
  outputfile->WriteTObject(HistReal01);
  outputfile->WriteTObject(Hist12);
  outputfile->WriteTObject(HistReal12);
  //outputfile->WriteTObject(HistLines);
  //outputfile->WriteTObject(HistRealLines);
  outputfile->WriteTObject(graphEff01, "GraphEff01");
  outputfile->WriteTObject(graphEff12, "GraphEff12");
  //outputfile->WriteTObject(graphEffLines, "GraphEffLines");

  //outputfile->Close();
}


 
void FillPtBins(std::vector<o2::its::Tracklet> tracklets, int layer, std::map<double,int> * PtBins, std::vector<double> * RawPt,
                std::map<o2::MCCompLabel,o2::its::ClusterMCLabelInfo>  LabelMap,  o2::its::Vertexer * vertexer, o2::its::ROframe frame){

  for(auto & tracklet : tracklets){
    auto Cluster = vertexer->getClusters()[layer][tracklet.firstClusterIndex] ;
    const auto MCLabel = frame.getClusterLabels(layer, Cluster.clusterId);
    if(LabelMap.find(MCLabel)!=LabelMap.end()){ //the cluster is in the LabelMap : it is real, we have to get the pT
      double Pt = LabelMap[MCLabel].Pt;
      RawPt->push_back(Pt);
      for ( std::map<double,int>::iterator it = PtBins->begin(); it != PtBins->end(); it++ ){
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