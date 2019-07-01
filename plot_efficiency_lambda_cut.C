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

    TFile *f = new TFile("lambda_cut_variation_phi_0_005.root");

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
    double CurrentMeanEffFake=0;
    double CurrentInvCut;
    int CountCurrentInvCut=0;
    std::vector<double> DistinctCuts;
    std::vector<double> MeanEffRecoCut;
    std::vector<double> MeanEffFakeCut;
    double meanEffReco;
    double meanEffFake;

    for(int i=0; i<nEntries; i++){
        if( std::find(DistinctCuts.begin(), DistinctCuts.end(), ArrInvCut[i])==DistinctCuts.end() || i==0){ //we have not processed this value
            CurrentInvCut=ArrInvCut[i];
            DistinctCuts.push_back(CurrentInvCut);
            CurrentMeanEffReco=0;
            CurrentMeanEffFake=0;
            CountCurrentInvCut=0;
            NumDistinctCuts++;
            for(int j=i; j<nEntries; j++){
                if(ArrInvCut[j]==CurrentInvCut){
                    CurrentMeanEffReco+=ArrEffReco[j];
                    CurrentMeanEffFake+= ArrEffFake[j];
                    CountCurrentInvCut++;
                }
            }
            meanEffReco = CurrentMeanEffReco/CountCurrentInvCut;
            meanEffFake = CurrentMeanEffFake/CountCurrentInvCut;
            std::cout<<"Cut :"<<1/CurrentInvCut<<"  mean reco : "<<meanEffReco<<"  mean fake : "<<meanEffFake<<std::endl;
            MeanEffRecoCut.push_back(meanEffReco);
            MeanEffFakeCut.push_back(meanEffFake);
        }else{
            continue;
        }  
    }



    TCanvas * c1=new TCanvas("c1","Efficiency Reconstruction",1000,700);
  // TH2F * hpx = new TH2F("hpx","Zoomed Graph Example",10,0,0.05,10,1,1);
   //hpx->SetStats(kFALSE);   // no statistics
   //hpx->Draw();


    c1->SetLogx();

    TGraph * graphReco= new TGraph(nEntries, ArrInvCut,ArrEffReco);
    graphReco->SetMarkerStyle(7);
    graphReco->GetXaxis()->SetTitle("1/TanLambaCut");
    graphReco->GetYaxis()->SetTitle("RecoMCValidated/Generated");
    graphReco->SetTitle("Reconstructed Monte Carlo validated tracklets");
    graphReco->Draw("AP");

    TGraph * graphMeanReco= new TGraph(NumDistinctCuts,DistinctCuts.data(),MeanEffRecoCut.data());
    graphMeanReco->SetMarkerStyle(20);
    graphMeanReco-> SetMarkerColor(2);
    graphMeanReco->Draw("P");


    TCanvas * c2=new TCanvas("c2","Fake tracklets",1000,700);

    c2->SetLogx();
    c2->SetLogy();

    TGraph * graphFake= new TGraph(nEntries, ArrInvCut,ArrEffFake);
    graphFake->SetMarkerStyle(7);
    graphFake->GetXaxis()->SetTitle("1/TanLambaCut");
    graphFake->GetYaxis()->SetTitle("Fake Tracklets/Generated");
    graphFake->SetTitle("Fake tracklets");
    graphFake->Draw("AP");

    TGraph * graphMeanFake= new TGraph(NumDistinctCuts,DistinctCuts.data(),MeanEffFakeCut.data());
    graphMeanFake->SetMarkerStyle(20);
    graphMeanFake-> SetMarkerColor(4);
    graphMeanFake->Draw("P");





    //std::cout<<"to draw : "<<NumDistinctCuts<<" "<<*DistinctCuts.data()<<"   "<<*MeanEffRecoCut.data()<<std::endl;
}
#endif