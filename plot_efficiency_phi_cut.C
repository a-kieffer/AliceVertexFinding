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
// plot this efficiency with respect to the phi cut parameter

void plot_efficiency_phi_cut();
void plot_quantiles(int nEntries, double * ArrInvCut, double *ArrEff,int * ArrNum, std::vector <int> CountCurrentInvCut,std::vector<double> DistinctCuts);
void connect_points(double *ArrEff, int * ArrNum, int Num, int NbCuts, std::vector <int>  CountCurrentInvCut, std::vector <double> DistinctCuts,
    int * FirstCutArrNum, int FirstCutEntries, int offset);




void plot_efficiency_phi_cut(){

    //name of the file where the TNtuples will be stored
    std::string file = "phi_cut_variation.root";
 
    //getting the data from the file
    TFile *f = new TFile(file.data());
    TNtuple * Tracklets01 ;
    TNtuple * Tracklets12 ;
    f->GetObject("Tracklets01", Tracklets01);
    f->GetObject("Tracklets12", Tracklets12);

    int nEntries=Tracklets01->GetEntries();

    TLeaf *leafReco01 = Tracklets01->GetLeaf("RecoMCvalid01");
    TLeaf *leafTot01 = Tracklets01->GetLeaf("TotTracklets01");
    TLeaf *leafReco12 = Tracklets12->GetLeaf("RecoMCvalid12");
    TLeaf *leafTot12 = Tracklets12->GetLeaf("TotTracklets12");
    TLeaf *leafCut = Tracklets01->GetLeaf("Cut");
    TLeaf *leafNum = Tracklets01->GetLeaf("EntryNum01"); 
    //this last leaf will be used in the plotting macro to connect point corresponding to the same ROFrame
    
    //arrays that contain the useful information to plot
    double ArrEffReco01[nEntries];
    double ArrEffReco12[nEntries];
    double ArrInvCut[nEntries];
    int ArrNum[nEntries];
    


	for (int i{0}; i < nEntries; ++i){
        Tracklets01->GetEntry(i);
        double RecoMCvalidated01 = leafReco01->GetValue(); 
        double TGenerated01 = leafTot01->GetValue(); 
        double Cut = leafCut->GetValue(); 
        double EntryNumber = leafNum->GetValue();

        Tracklets12->GetEntry(i);
        double RecoMCvalidated12 = leafReco12->GetValue(); 
        double TGenerated12 = leafTot12->GetValue(); 

        double EffReco01=RecoMCvalidated01/TGenerated01;
        double EffReco12=RecoMCvalidated12/TGenerated12;


        ArrEffReco01[i]=EffReco01;
        ArrEffReco12[i]=EffReco12;
        ArrInvCut[i]=1/Cut;
        ArrNum[i]=EntryNumber;
	}


  //sorting the vectors to plot a curve later : sorting according to the Cut value
    

    std::vector<std::size_t> index_vec;
    for (std::size_t i = 0; i < (unsigned int)nEntries; ++i) { index_vec.push_back(i); }
    std::stable_sort( index_vec.begin(), index_vec.end(), [&](std::size_t a, std::size_t b) { return ArrInvCut[a] < ArrInvCut[b]; });

    //the same arrays as before, but sorted according to the Cut value
    double SortedArrEffReco01 [nEntries];
    double SortedArrEffReco12 [nEntries];
    double SortedArrInvCut [nEntries];
    int SortedArrNum [nEntries];

    for (std::size_t i = 0; i < (unsigned int) nEntries; ++i)
    {
        SortedArrEffReco01[i]=ArrEffReco01[index_vec[i]];
        SortedArrEffReco12[i]=ArrEffReco12[index_vec[i]];
        SortedArrInvCut[i]=ArrInvCut[index_vec[i]];
        SortedArrNum[i]=ArrNum[index_vec[i]];
    }

    // Computing the mean for each  Cut value to plot it 

    int NumDistinctCuts=0; //number of distinct cut values
    double CurrentMeanEffReco01=0; //mean of the efficiency for the Cut we are currently working on
    double CurrentMeanEffReco12=0;
    double CurrentInvCut;

    std::vector <int> CountCurrentInvCut; //contains the number of entries for each cut value
    std::vector<double> DistinctCuts; //contains the distinct cuts, not sorted
    std::vector<double> MeanEffReco01Cut; //contains the mean Efficiency for each cut
    std::vector<double> MeanEffReco12Cut; 
    double meanEffReco01;
    double meanEffReco12;


    //this for loop will iterate over all the different cut values and then over all the entries corresponding to the cut
 
    for(int i=0; i<nEntries; i++){
        if( std::find(DistinctCuts.begin(), DistinctCuts.end(), SortedArrInvCut[i])==DistinctCuts.end() || i==0){ //we have not processed this value yet

            CountCurrentInvCut.push_back(0); //always use push back for the 1st element of the vector
            CurrentInvCut=SortedArrInvCut[i];
            DistinctCuts.push_back(CurrentInvCut);
            CurrentMeanEffReco01=0;
            CurrentMeanEffReco12=0;
            for(int j=i; j<nEntries; j++){ //all the entries after this one
                if(SortedArrInvCut[j]==CurrentInvCut){ //if the cut value is the same
                    CurrentMeanEffReco01+=SortedArrEffReco01[j];
                    CurrentMeanEffReco12+=SortedArrEffReco12[j];
                    CountCurrentInvCut[NumDistinctCuts]++; 
                }
            }
            meanEffReco01 = CurrentMeanEffReco01/CountCurrentInvCut[NumDistinctCuts];
            meanEffReco12 = CurrentMeanEffReco12/CountCurrentInvCut[NumDistinctCuts];
            MeanEffReco01Cut.push_back(meanEffReco01);
            MeanEffReco12Cut.push_back(meanEffReco12);
            NumDistinctCuts++; 
            std::cout<<"Cut :"<<1/CurrentInvCut<<"     Efficiency 01 :"<<meanEffReco01<<"      Efficiency 12 :"<<meanEffReco12<<std::endl;
        }else{ //we have already processed this cut
            continue;
        }  
    }

    
   
    for (int c : CountCurrentInvCut){
        std::cout<<"Number of Entries :"<<c<<std::endl;
    }
    
     
    TCanvas * c1=new TCanvas("c1","Efficiency Reconstruction",1000,700);

    c1->Divide(2,1);
    c1->cd(1)->SetLogx();

    TGraph * graphReco01= new TGraph(nEntries, ArrInvCut,ArrEffReco01);
    graphReco01->SetMarkerStyle(7);
    graphReco01->GetXaxis()->SetTitle("1/PhiCut");
    graphReco01->GetYaxis()->SetTitle("Validated tracklets 01/Total");
    graphReco01->SetTitle("Reconstructed Monte Carlo validated tracklets 01");
    graphReco01->Draw("AP");

    TGraph * graphMeanReco01= new TGraph(NumDistinctCuts,DistinctCuts.data(),MeanEffReco01Cut.data());
    graphMeanReco01->SetMarkerStyle(20);
    graphMeanReco01-> SetMarkerColor(2);
    graphMeanReco01->Draw("PL");

    plot_quantiles(nEntries, SortedArrInvCut, SortedArrEffReco01,  SortedArrNum, CountCurrentInvCut,DistinctCuts);

    c1->cd(2)->SetLogx();
    
    TGraph * graphReco12= new TGraph(nEntries, ArrInvCut,ArrEffReco12);
    graphReco12->SetMarkerStyle(7);
    graphReco12->GetXaxis()->SetTitle("1/PhiCut");
    graphReco12->GetYaxis()->SetTitle("Validated tracklets 12/Total");
    graphReco12->SetTitle("Reconstructed Monte Carlo validated tracklets 12");
    graphReco12->Draw("AP");

    TGraph * graphMeanReco12= new TGraph(NumDistinctCuts,DistinctCuts.data(),MeanEffReco12Cut.data());
    graphMeanReco12->SetMarkerStyle(20);
    graphMeanReco12-> SetMarkerColor(2);
    graphMeanReco12->Draw("PL");

    plot_quantiles(nEntries, SortedArrInvCut, SortedArrEffReco12,  SortedArrNum, CountCurrentInvCut,DistinctCuts);


    //saving the plots into a file

    TFile * outputfile = new TFile("efficiency_phi_cut_plots.root", "recreate");
    outputfile->WriteTObject(graphReco01, "GraphReconstruction01");
    outputfile->WriteTObject(graphMeanReco01, "GraphMeanReconstruction01");
    outputfile->WriteTObject(graphReco12, "GraphReconstruction12");
    outputfile->WriteTObject(graphMeanReco12, "GraphMeanReconstruction12");
    outputfile->Close();

}
 

    // this function selects 5 ROFrames and then connects the points corresponding to them

    void plot_quantiles(int nEntries, double * ArrInvCut, double *ArrEff, int * ArrNum, std::vector <int> CountCurrentInvCut,
     std::vector<double> DistinctCuts){

        //the selected ROFrames correspond to quantiles for the 1st cut value, this could be changed
        // the goal is to find the ROFrame Number of all the quantiles for the 1st cut value and to then get all the values for this number
        // an easy way is to sort the efficiency array for the 1st cut value

        int FirstCutEntries = CountCurrentInvCut[0];
        int offset = FirstCutEntries/5;

        double FirstCutArrEff [FirstCutEntries ]; //array containing the efficiency
        int FirstCutArrNum [FirstCutEntries ]; //array containing the number of the ROFrame associated with the efficiency
        
        // we sort the entries corresponding to the 1st cut by efficiency value
        std::vector<std::size_t> index_vec;
        for (std::size_t i = 0; i < (unsigned int) FirstCutEntries ; ++i) { index_vec.push_back(i); }
        std::sort( index_vec.begin(), index_vec.end(), [&](std::size_t a, std::size_t b) { return ArrEff[a] < ArrEff[b]; });
        for (std::size_t i = 0; i < (unsigned int) FirstCutEntries ; ++i){
            FirstCutArrEff[i]=ArrEff[index_vec[i]];
            FirstCutArrNum[i]=ArrNum[index_vec[i]];
         }

        int Num1= FirstCutArrNum[0]; //Number of the entry with the lowset Efficiency with the first lambda cut value
        int Num5= FirstCutArrNum[FirstCutEntries-3]; //close to the highest efficiency
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
        std::vector <double> Vector; //vector to store the points that will be connected
        int jOffset=0;
        std::vector <double> DistinctCutsPlot;

            for(int i=0; i<NbCuts; i++){
                for(int j=jOffset; j<jOffset+CountCurrentInvCut[i]; j++){ //we only look at one cut value
                    if(ArrNum[j]==Num){ //this is the entry related to this number
                        Vector.push_back(ArrEff[j]);
                        DistinctCutsPlot.push_back(DistinctCuts[i]);
                        std::cout<<"To plot :   Cut :"<<DistinctCuts[i]<<"  Eff "<<ArrEff[j]<<"     Entry Num : "<<Num<<std::endl;
                        break;
                    } 
                }
                jOffset+=CountCurrentInvCut[i];
            }
        TGraph * graph1= new TGraph(Vector.size(), DistinctCutsPlot.data(),Vector.data());
        graph1->Draw("PL");
    }
    
    
#endif