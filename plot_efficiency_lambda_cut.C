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


void plot_efficiency_lambda_cut(float PhiAngle = 0.005f){

    std::string Phi = std::to_string(PhiAngle);
    std::size_t pos = Phi.find(".");
    if (pos == std::string::npos){

    }else{
        Phi.replace(pos, 1, "_");
    }
    std::string StdOutfile = "lambda_cut_variation";
    std::string outfile = StdOutfile + "_phi_" + Phi +".root";
    std::cout<<outfile<<std::endl;

    outfile = "lambda_cut_variation_phi_0_05.root";

    TFile *f = new TFile(outfile.data());

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

  //sorting the vectors to plot a curve later
    std::vector<std::size_t> index_vec;

    for (std::size_t i = 0; i < nEntries; ++i) { index_vec.push_back(i); }

    std::sort( index_vec.begin(), index_vec.end(), [&](std::size_t a, std::size_t b) { return ArrInvCut[a] < ArrInvCut[b]; });

    double SortedArrEffReco [nEntries];
    double SortedArrEffFake [nEntries];
    double SortedArrInvCut [nEntries];

    for (std::size_t i = 0; i <nEntries; ++i)
    {
        SortedArrEffReco[i]=ArrEffReco[index_vec[i]];
        SortedArrEffFake[i]=ArrEffFake[index_vec[i]];
        SortedArrInvCut[i]=ArrInvCut[index_vec[i]];
    }

    

    int NumDistinctCuts=0;
    double CurrentMeanEffReco=0;
    double CurrentMeanEffFake=0;
    double CurrentInvCut;

    std::vector <int> CountCurrentInvCut;
    std::vector<double> DistinctCuts;
    std::vector<double> MeanEffRecoCut;
    std::vector<double> MeanEffFakeCut;
    double meanEffReco;
    double meanEffFake;
 
    for(int i=0; i<nEntries; i++){
        if( std::find(DistinctCuts.begin(), DistinctCuts.end(), SortedArrInvCut[i])==DistinctCuts.end() || i==0){ //we have not processed this value
            //std:cout<<NumDistinctCuts<<std::endl;
        
            //CountCurrentInvCut[NumDistinctCuts]=0; //this line has a problem
            CountCurrentInvCut.push_back(0);
            
            CurrentInvCut=SortedArrInvCut[i];
            DistinctCuts.push_back(CurrentInvCut);
            CurrentMeanEffReco=0;
            CurrentMeanEffFake=0;
            for(int j=i; j<nEntries; j++){
                if(SortedArrInvCut[j]==CurrentInvCut){
                    //ok works
                    CurrentMeanEffReco+=SortedArrEffReco[j];
                    CurrentMeanEffFake+= SortedArrEffFake[j];
                    CountCurrentInvCut[NumDistinctCuts]++; 
                }
            }
            
            meanEffReco = CurrentMeanEffReco/CountCurrentInvCut[NumDistinctCuts];
            meanEffFake = CurrentMeanEffFake/CountCurrentInvCut[NumDistinctCuts];
            //std::cout<<"Cut :"<<1/CurrentInvCut<<"  mean reco : "<<meanEffReco<<"  mean fake : "<<meanEffFake<<std::endl;
            MeanEffRecoCut.push_back(meanEffReco);
            MeanEffFakeCut.push_back(meanEffFake);
            NumDistinctCuts++; 
        }else{
            continue;
        }  
    }

    
    for (int c : CountCurrentInvCut){
        std:cout<<c<<std::endl;
    }
    /*
    TCanvas * c1=new TCanvas("c1","Efficiency Reconstruction",1000,700);

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
    graphMeanReco->Draw("PC");


    TCanvas * c2=new TCanvas("c2","Fake tracklets",1000,700);

    c2->SetLogx();
    //c2->SetLogy();

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

    */

    //std::cout<<"to draw : "<<NumDistinctCuts<<" "<<*DistinctCuts.data()<<"   "<<*MeanEffRecoCut.data()<<std::endl;
}
#endif