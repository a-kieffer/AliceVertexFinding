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



void basic_trackleter(       const int numPhiDivs = 16,
			     const int numZDivs = 0,
                             const int inspEvt = -1,
                             const int numEvents = 1,
                             const std::string inputClustersITS = "o2clus_its.root",
                             const std::string inputGRP = "o2sim_grp.root",
                             const std::string simfilename = "o2sim.root",
                             const std::string paramfilename = "O2geometry.root",
                             const std::string path = "./")
{


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


  std::uint32_t roFrame = 0;
  const int stopAt = (inspEvt == -1) ? rofs->size() : inspEvt + numEvents;

  o2::its::ROframe frame(-123);
  o2::its::VertexerTraits* traits = nullptr;

  traits = o2::its::createVertexerTraits();
  o2::its::Vertexer vertexer(traits);

  int EntryNum = (inspEvt ==-1) ? 0 : inspEvt;

  struct o2::its::VertexingParameters  par;
  // par.phiCut=3.f;
  // par.tanLambdaCut=0.025;
  for (int iRof = (inspEvt == -1) ? 0 : inspEvt; iRof < stopAt; ++iRof) {
    auto rof = (*rofs)[iRof];

    std::cout << "Entry: " << EntryNum << std::endl;
    EntryNum++;
    o2::its::ioutils::generateSimpleData(frame, numPhiDivs, numZDivs);
    vertexer.initialiseVertexer(&frame);
    // vertexer.setParameters(par);

    // vertexer.findTrivialMCTracklets();
    vertexer.findTracklets(false);
    int T01size = vertexer.getTracklets01().size();
    int T12size = vertexer.getTracklets12().size();
    int Linesize = vertexer.getLines().size();

    std::cout<<"Tracklets 01 : "<<T01size<<"  Tracklets 12: "<<T12size<<"  Lines:" <<Linesize<<std::endl;


    vertexer.processLines();
    vertexer.findVertices();
    vertexer.dumpTraits();
    auto vertices = vertexer.exportVertices();
    if ( vertices.size() > 0 ) std::cout<<"x : "<<vertices[0].getX()<<" y : "<<vertices[0].getY()<<" z "<<vertices[0].getX()<<std::endl;
    //std::vector<std::array<float, 4>> centroidsData = vertexer.getCentroids();
    //std::cout<<"x : "<<centroidsData[0][0]<<" y : "<<centroidsData[0][1]<<" z "<<centroidsData[0][2]<<std::endl;
  }

}
#endif
