#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "iostream"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include <TNtuple.h>
#include <TObject.h>

#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPObject.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"
#include "ITStracking/Vertexer.h"


using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

//the aim of this macro is to compute tracklets using a tanLambdaCut value and then writes the number of real trackets found, of fakes tracklets 
//and of the cut value in a file 



void trackleter_phi_cut( 
                             const int inspEvt = -1,
                             const int numEvents = 1,
                             const std::string inputClustersITS = "o2clus_its.root",
                             const std::string inputGRP = "o2sim_grp.root",
                             const std::string simfilename = "o2sim.root",
                             const std::string paramfilename = "O2geometry.root",
                             const std::string path = "./")
{
  
  std::string outfile = "phi_cut_variation.root";

  std::cout<<outfile<<std::endl;
  
  const auto grp = o2::parameters::GRPObject::loadFrom(path + inputGRP);
  const bool isITS = grp->isDetReadOut(o2::detectors::DetID::ITS);
  const bool isContITS = grp->isDetContinuousReadOut(o2::detectors::DetID::ITS);
  std::cout << "ITS is in " << (isContITS ? "CONTINUOS" : "TRIGGERED") << " readout mode" << std::endl;
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



  TFile* outputfile = new TFile(outfile.data(), "update");
  
  TNtuple * Tracklets01 = nullptr;
  TNtuple * Tracklets12 = nullptr;
  
  outputfile->GetObject("Tracklets01", Tracklets01);
  outputfile->GetObject("Tracklets12", Tracklets12);


  if(Tracklets01==nullptr || Tracklets12==nullptr){
    std::cout<<"No objects to read from the file : they will be created \n";
    // here it is important to use dynamic allocation so the object does not disappear after the block
    Tracklets01 = new TNtuple("Tracklets01", "Tracklets01", "RecoMCvalid01:TotTracklets01:Cut:EntryNum01");
    Tracklets12 = new TNtuple("Tracklets12", "Tracklets12", "RecoMCvalid12:TotTracklets12:Cut:EntryNum12");
  }

 if(Tracklets01==nullptr || Tracklets12==nullptr){
    std::cout<<"Still no object \n";
    exit(0);
  }


std::vector<double> Cut{ 0.02, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.002, 1};
  //0.001,0.005}; //until 0.5

for (double cut : Cut){


  std::uint32_t roFrame = 0;
  const int stopAt = (inspEvt == -1) ? rofs->size() : inspEvt + numEvents;
  o2::its::ROframe frame(-123);
  o2::its::VertexerTraits* traits = nullptr;
  traits = o2::its::createVertexerTraits();
 
 
  o2::its::Vertexer vertexer(traits);
  int EntryNum = (inspEvt ==-1) ? 0 : inspEvt;

////////////////////// setting the cut value
  struct o2::its::VertexingParameters  par;
  par.phiCut=cut;
 

  itsClusters.GetEntry(0);
  mcHeaderTree.GetEntry(0);
  
  for (int iRof = (inspEvt == -1) ? 0 : inspEvt; iRof < stopAt; ++iRof) {

    int TGenerated01=0;
    int RecoMCvalidated01=0;
    int TGenerated12=0;
    int RecoMCvalidated12=0;

    auto rof = (*rofs)[iRof];

    std::cout << "Entry: " << EntryNum << std::endl;
    
   
    int nclUsed = o2::its::ioutils::loadROFrameData(rof, frame, clusters, labels);
    vertexer.initialiseVertexer(&frame);

    vertexer.setParameters(par); //to set Cut

    vertexer.findTracklets(true); //this is where the tracklets are created and then selected
    // the selected tracklets are stored in mTracklets

    //this is the mTracklets vector, which contains tracklets and the 2 clusters they contain
    // it is a vector of lines

    
    RecoMCvalidated01 += vertexer.getTracklets01().size(); //gets all the tracklets
    RecoMCvalidated12 += vertexer.getTracklets12().size(); 

    vertexer.initialiseVertexer(&frame);
   
   
    vertexer.findTrivialMCTracklets();
    TGenerated01 += vertexer.getTracklets01().size();
    TGenerated12 += vertexer.getTracklets12().size();

    float phiCut = vertexer.getVertParameters().phiCut;


    if(RecoMCvalidated01!=0 && TGenerated01!=0 && RecoMCvalidated12!=0 && TGenerated12!=0){ //do not write entries that are 0
      Tracklets01->Fill(RecoMCvalidated01,TGenerated01,phiCut,EntryNum);
      Tracklets12->Fill(RecoMCvalidated12,TGenerated12,phiCut,EntryNum);
      
    }
    std::cout<<"Tracklets reconstructed 01 :"<<RecoMCvalidated01<<"   Total real tracklets :"<<TGenerated01<<"  Cut :"<<phiCut<<std::endl;
    std::cout<<"Tracklets reconstructed 12 :"<<RecoMCvalidated12<<"   Total real tracklets :"<<TGenerated12<<"  Cut :"<<phiCut<<std::endl;
    ++EntryNum;
  }


  Tracklets01->Write(0,TObject::kOverwrite);
  Tracklets12->Write(0,TObject::kOverwrite);

}

  outputfile->Close();


  //return 0;
}


#endif



    