#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "iostream"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include <TNtuple.h>
#include <TLeaf.h>
#include <TGraph.h>
#include <TH2F.h>
#include "TH1F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"


#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPObject.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"
#include "ITStracking/Vertexer.h"

#include "ITStracking/ClusterLines.h"
#include "ITStracking/Tracklet.h"
#include "ITStracking/Cluster.h"

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

//the aim of this macro is to read data from root files containing the reconstructed tracklets and compare them to the real tracklets to compute an efficiency and 
// plot this efficiency with respect to the delta tan lambda cut parameter


void plot_efficiency_lambda_cut(){

    TFile *f = new TFile("lambda_cut_variation_phi_0_05.root");

    TNtuple * data ;
    f->GetObject("Results", data);

    int nEntries=data->GetEntries();

    TLeaf *leafFake = data->GetLeaf("FakeTracklets");
    TLeaf *leafReco = data->GetLeaf("RecoMCvalid");
    TLeaf *leafTot = data->GetLeaf("TotTracklets");
    TLeaf* leafCut = data->GetLeaf("Cut");

    double ArrEffReco [nEntries];
    double ArrEffFake [nEntries];
    double ArrInvCut [nEntries];

	for (int i{0}; i < nEntries; ++i){

        data->GetEntry(i);
        double FakeTracklets = leafFake->GetValue();
        double RecoMCvalidated = leafReco->GetValue(); 
        double TGenerated = leafTot->GetValue(); 
        double lambdaCut = leafCut->GetValue(); 

        double EffReco=RecoMCvalidated/TGenerated;
        double EffFake= FakeTracklets/TGenerated;

        ArrEffReco[i]=EffReco;
        ArrEffFake[i]= EffFake;
        ArrInvCut[i]=1/lambdaCut;

        //std::cout<<FakeTracklets<<"\n";
        //std::cout<<EffReco<<"\n"
        

	}

    int NumDistinctCuts=0;
    double CurrentMeanEffReco=0;
    double CurrentInvCut;
    int CountCurrentInvCut=0;
    std::vector<double> DistinctCuts;
    std::vector<double> MeanEffRecoCut;
    double mean;

    for(int i=0; i<nEntries; i++){
        if( std::find(DistinctCuts.begin(), DistinctCuts.end(), ArrInvCut[i])==DistinctCuts.end() || i==0){ //we have not processed this value
            CurrentInvCut=ArrInvCut[i];
            DistinctCuts.push_back(CurrentInvCut);
            CurrentMeanEffReco=0;
            CountCurrentInvCut=0;
            NumDistinctCuts++;
            for(int j=i; j<nEntries; j++){
                if(ArrInvCut[j]==CurrentInvCut){
                    CurrentMeanEffReco+=ArrEffReco[j];
                    CountCurrentInvCut++;
                }
            }
            mean = CurrentMeanEffReco/CountCurrentInvCut;
            std::cout<<"Cut :"<<1/CurrentInvCut<<"  mean : "<<mean<<std::endl;
            MeanEffRecoCut.push_back(mean);
        }else{
            continue;
        }  
    }



    TCanvas * c1=new TCanvas("c1","A Zoomed Graph",1000,700);
  // TH2F * hpx = new TH2F("hpx","Zoomed Graph Example",10,0,0.05,10,1,1);
   //hpx->SetStats(kFALSE);   // no statistics
   //hpx->Draw();


    c1->SetLogx();

    TGraph * graph= new TGraph(nEntries, ArrInvCut,ArrEffReco);
    graph->SetMarkerStyle(7);
    graph->Draw("AP");

    TGraph * graphMean= new TGraph(NumDistinctCuts,DistinctCuts.data(),MeanEffRecoCut.data());
    graphMean->SetMarkerStyle(20);
    graphMean-> SetMarkerColor(2);
    graph->GetXaxis()->SetTitle("1/TanLambaCut");
    graph->GetYaxis()->SetTitle("RecoMCValidated/Generated");

    graphMean->Draw("P");

    //std::cout<<"to draw : "<<NumDistinctCuts<<" "<<*DistinctCuts.data()<<"   "<<*MeanEffRecoCut.data()<<std::endl;
}
#endif