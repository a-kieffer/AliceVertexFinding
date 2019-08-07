#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "iostream"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include <TNtuple.h>
#include <TLeaf.h>
#include <TGraph.h>
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

void plot_quantiles(int nEntries, double * ArrInvCut, double *ArrEff,int * ArrNum, std::vector <int> CountCurrentInvCut,std::vector<double> DistinctCuts);
void plot_efficiency_lambda_cut(float PhiAngle);
void connect_points(double *ArrEff, int * ArrNum, int Num, int NbCuts, std::vector <int>  CountCurrentInvCut, std::vector <double> DistinctCuts,
    int * FirstCutArrNum, int FirstCutEntries, int offset);


void plot_efficiency_lambda_cut(float PhiAngle = 0.005f){

    //setting the title of the file
    std::string Phi = std::to_string(PhiAngle);
    std::size_t pos = Phi.find(".");
    if (pos == std::string::npos){
    }else{
        Phi.replace(pos, 1, "_");
    }
    std::string Stdfile = "lambda_cut_variation";
    std::string file = Stdfile + "_phi_" + Phi +".root";
    std::cout<<"file to be read :"<<file<<std::endl;

    //getting the data from the file
    TFile *f = new TFile(file.data());
    TNtuple * data ;
    f->GetObject("Results", data);
    int nEntries=data->GetEntries();

    TLeaf *leafFake = data->GetLeaf("FakeTracklets");
    TLeaf *leafReco = data->GetLeaf("RecoMCvalid");
    TLeaf *leafTot = data->GetLeaf("TotTracklets");
    TLeaf *leafCut = data->GetLeaf("Cut");
    TLeaf *leafNum = data->GetLeaf("EntryNum");
    
    double ArrEffReco[nEntries];
    double ArrEffFake[nEntries];
    double ArrInvCut[nEntries];
    int ArrNum[nEntries];

    std::vector <int> alertReco;
    
	for (int i{0}; i < nEntries; ++i){
        data->GetEntry(i);
        double FakeTracklets = leafFake->GetValue();
        double RecoMCvalidated = leafReco->GetValue(); 
        double TGenerated = leafTot->GetValue(); 
        double lambdaCut = leafCut->GetValue(); 
        double EntryNumber = leafNum->GetValue();

        double EffReco=RecoMCvalidated/TGenerated;
        double EffFake= FakeTracklets/TGenerated;

        if(EffReco>1){
                alertReco.push_back(EntryNumber);
            }

        ArrEffReco[i]=EffReco;
        ArrEffFake[i]= EffFake;
        ArrInvCut[i]=1/lambdaCut;
        ArrNum[i]=EntryNumber;
	}


  //sorting the vectors to plot a curve later : sorting according to the Tan Lambda Cut value
    

    std::vector<std::size_t> index_vec;
    for (std::size_t i = 0; i < (unsigned int)nEntries; ++i) { index_vec.push_back(i); }
    std::stable_sort( index_vec.begin(), index_vec.end(), [&](std::size_t a, std::size_t b) { return ArrInvCut[a] < ArrInvCut[b]; });

    double SortedArrEffReco [nEntries];
    double SortedArrEffFake [nEntries];
    double SortedArrInvCut [nEntries];
    int SortedArrNum [nEntries];

    for (std::size_t i = 0; i < (unsigned int) nEntries; ++i)
    {
        SortedArrEffReco[i]=ArrEffReco[index_vec[i]];
        SortedArrEffFake[i]=ArrEffFake[index_vec[i]];
        SortedArrInvCut[i]=ArrInvCut[index_vec[i]];
        SortedArrNum[i]=ArrNum[index_vec[i]];
    }

    // Computing the mean for each Tan Lambda Cut value to plot it 

    int NumDistinctCuts=0;
    double CurrentMeanEffReco=0;
    double CurrentMeanEffFake=0;
    double CurrentInvCut;

    std::vector <int> CountCurrentInvCut; //contains the number of entries for each cut value
    std::vector<double> DistinctCuts; //contains the distinct cuts, not sorted
    std::vector<double> MeanEffRecoCut; //contains the mean EffReco for each cut
    std::vector<double> MeanEffFakeCut;
    double meanEffReco;
    double meanEffFake;
 
    for(int i=0; i<nEntries; i++){
        if( std::find(DistinctCuts.begin(), DistinctCuts.end(), SortedArrInvCut[i])==DistinctCuts.end() || i==0){ //we have not processed this value yet

            CountCurrentInvCut.push_back(0); //always use push back for the 1st element of the vector
            CurrentInvCut=SortedArrInvCut[i];
            DistinctCuts.push_back(CurrentInvCut);
            CurrentMeanEffReco=0;
            CurrentMeanEffFake=0;
            for(int j=i; j<nEntries; j++){ //all the entries after this one
                if(SortedArrInvCut[j]==CurrentInvCut){ //if the cut value is the same
                    CurrentMeanEffReco+=SortedArrEffReco[j];
                    CurrentMeanEffFake+= SortedArrEffFake[j];
                    CountCurrentInvCut[NumDistinctCuts]++; 
                }
            }
            
            meanEffReco = CurrentMeanEffReco/CountCurrentInvCut[NumDistinctCuts];
            meanEffFake = CurrentMeanEffFake/CountCurrentInvCut[NumDistinctCuts];
            MeanEffRecoCut.push_back(meanEffReco);
            MeanEffFakeCut.push_back(meanEffFake);
            NumDistinctCuts++; 
            std::cout<<"Cut :"<<1/CurrentInvCut<<"     Reco efficiency :"<<meanEffReco<<"      Fake tracklets :"<<meanEffFake<<std::endl;
        }else{ //we have already processed this cut
            continue;
        }  
    }

    
   
    for (int c : CountCurrentInvCut){
        std:cout<<"Number of Entries :"<<c<<std::endl;
    }
    
     
    TCanvas * c1=new TCanvas("c1","Efficiency Reconstruction",1000,700);

    c1->SetLogx();

    TGraph * graphReco= new TGraph(nEntries, ArrInvCut,ArrEffReco);
    graphReco->SetMarkerStyle(7);
    graphReco->GetXaxis()->SetTitle("1/TanLambaCut");
    graphReco->GetYaxis()->SetTitle("Validated tracklets/Total");
    graphReco->SetTitle("Selected Monte Carlo validated tracklets");
    graphReco->Draw("AP");

    TGraph * graphMeanReco= new TGraph(NumDistinctCuts,DistinctCuts.data(),MeanEffRecoCut.data());
    graphMeanReco->SetMarkerStyle(20);
    graphMeanReco-> SetMarkerColor(2);
    graphMeanReco->Draw("PL");

    plot_quantiles(nEntries, SortedArrInvCut, SortedArrEffReco,  SortedArrNum, CountCurrentInvCut,DistinctCuts);


    TCanvas * c2=new TCanvas("c2","Fake tracklets",1000,700);

    c2->SetLogx();
    
    TGraph * graphFake= new TGraph(nEntries, ArrInvCut,ArrEffFake);
    graphFake->SetMarkerStyle(7);
    graphFake->GetXaxis()->SetTitle("1/TanLambaCut");
    graphFake->GetYaxis()->SetTitle("Fake Tracklets/Total");
    graphFake->SetTitle("Fake tracklets");
    graphFake->Draw("AP");

    plot_quantiles(nEntries, SortedArrInvCut, SortedArrEffFake,  SortedArrNum, CountCurrentInvCut,DistinctCuts);

    TGraph * graphMeanFake= new TGraph(NumDistinctCuts,DistinctCuts.data(),MeanEffFakeCut.data());
    graphMeanFake->SetMarkerStyle(20);
    graphMeanFake-> SetMarkerColor(4);
    graphMeanFake->Draw("PL");

     std::cout<<" Number of the entries for Reco : "<<alertReco.size()<<std::endl;
    for (int entry : alertReco){
        std::cout<<entry<<std::endl;
    }


    //saving the plots into a file

    TFile * outputfile = new TFile("efficiency_lambda_cut_plots.root", "recreate");
    outputfile->WriteTObject(graphReco, "GraphReconstruction");
    outputfile->WriteTObject(graphMeanReco, "GraphMeanReconstruction");
    outputfile->WriteTObject(graphFake, "GraphFakeTracklets");
    outputfile->WriteTObject(graphMeanFake, "GraphMeanFakeTracklets");
    outputfile->Close();

}



    void plot_quantiles(int nEntries, double * ArrInvCut, double *ArrEff, int * ArrNum, std::vector <int> CountCurrentInvCut,
     std::vector<double> DistinctCuts){

        // the goal is to find the Number of all the quantiles for the 1st cut value and to then get all the values for this number
        // an easy way is to sort the efficiency array for the 1st cut value

        int FirstCutEntries = CountCurrentInvCut[0];
        int offset = FirstCutEntries/5;

        double FirstCutArrEff [FirstCutEntries ];
        int FirstCutArrNum [FirstCutEntries ];
        

        std::vector<std::size_t> index_vec;
        for (std::size_t i = 0; i < (unsigned int) FirstCutEntries ; ++i) { index_vec.push_back(i); }
        std::sort( index_vec.begin(), index_vec.end(), [&](std::size_t a, std::size_t b) { return ArrEff[a] < ArrEff[b]; });
        for (std::size_t i = 0; i < (unsigned int) FirstCutEntries ; ++i){
            FirstCutArrEff[i]=ArrEff[index_vec[i]];
            FirstCutArrNum[i]=ArrNum[index_vec[i]];
         }

        int Num1= FirstCutArrNum[0]; //Number of the entry with the lowset Efficiency with the first lambda cut value
        int Num5= FirstCutArrNum[FirstCutEntries-1];
        int Num2= FirstCutArrNum[offset] ;
        int Num3= FirstCutArrNum[2*offset] ;
        int Num4 = FirstCutArrNum[FirstCutEntries-1-offset];

        int NbCuts = DistinctCuts.size();
     
        std::sort(DistinctCuts.begin(), DistinctCuts.end()); //we sort it for the plots later

        connect_points(ArrEff, ArrNum, Num1,  NbCuts, CountCurrentInvCut,DistinctCuts, FirstCutArrNum, FirstCutEntries, offset);
        connect_points(ArrEff, ArrNum, Num2,  NbCuts, CountCurrentInvCut,DistinctCuts, FirstCutArrNum, FirstCutEntries, offset);
        connect_points(ArrEff, ArrNum, Num3,  NbCuts, CountCurrentInvCut,DistinctCuts, FirstCutArrNum, FirstCutEntries, offset);
        connect_points(ArrEff, ArrNum, Num4,  NbCuts, CountCurrentInvCut,DistinctCuts, FirstCutArrNum, FirstCutEntries, offset);
        connect_points(ArrEff, ArrNum, Num5,  NbCuts, CountCurrentInvCut,DistinctCuts, FirstCutArrNum, FirstCutEntries, offset);
    }
    

    void connect_points(double *ArrEff, int * ArrNum, int Num, int NbCuts, std::vector <int>  CountCurrentInvCut, std::vector <double> DistinctCuts,
    int * FirstCutArrNum, int FirstCutEntries, int offset){
        std::vector <double> Vector;
        int jOffset=0;
        std::vector <double> DistinctCutsPlot;
        
        //std::cout<<"vector size :"<<Vector.size()<<std::endl;

            for(int i=0; i<NbCuts; i++){
                for(int j=jOffset; j<jOffset+CountCurrentInvCut[i]; j++){ //we only look at one tan lambda cut value
                    if(ArrNum[j]==Num){ //this is the entry related to this number
                        Vector.push_back(ArrEff[j]);
                        DistinctCutsPlot.push_back(DistinctCuts[i]);
                        std::cout<<"To plot :   Cut :"<<DistinctCuts[i]<<"  Eff "<<ArrEff[j]<<"     Entry Num : "<<Num<<std::endl;
                        break;
                    } 
                }
                jOffset+=CountCurrentInvCut[i];
            }
        //std::cout<<"vector size after :"<<Vector.size()<<std::endl;
        TGraph * graph1= new TGraph(Vector.size(), DistinctCutsPlot.data(),Vector.data());
        graph1->Draw("PL");
    }
#endif