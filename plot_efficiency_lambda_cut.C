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

void plot_quantiles(int nEntries, double * ArrInvCut, double *ArrEffReco, double *ArrEffFake,int * ArrNum, std::vector <int> CountCurrentInvCut,std::vector<double> DistinctCuts);
void plot_efficiency_lambda_cut(float PhiAngle);
void connect_points(double *ArrEffReco, int * ArrNum, int Num, int NbCuts, std::vector <int>  CountCurrentInvCut, std::vector <double> DistinctCuts);




void plot_efficiency_lambda_cut(float PhiAngle = 0.005f){

    //setting the title of the file
    std::string Phi = std::to_string(PhiAngle);
    std::size_t pos = Phi.find(".");
    if (pos == std::string::npos){

    }else{
        Phi.replace(pos, 1, "_");
    }
    std::string StdOutfile = "lambda_cut_variation";
    std::string outfile = StdOutfile + "_phi_" + Phi +".root";
    std::cout<<outfile<<std::endl;

    //outfile = "lambda_cut_variation_phi_0_05.root";

    //getting the data from the file

    TFile *f = new TFile(outfile.data());

    TNtuple * data ;
    f->GetObject("Results", data);

    int nEntries=data->GetEntries();

    TLeaf *leafFake = data->GetLeaf("FakeTracklets");
    TLeaf *leafReco = data->GetLeaf("RecoMCvalid");
    TLeaf *leafTot = data->GetLeaf("TotTracklets");
    TLeaf *leafCut = data->GetLeaf("Cut");
    TLeaf *leafNum = data->GetLeaf("EntryNum");
  
    
    double ArrEffReco [nEntries];
    double ArrEffFake [nEntries];
    double ArrInvCut [nEntries];
    int ArrNum[nEntries];
    

	for (int i{0}; i < nEntries; ++i){

        data->GetEntry(i);
        double FakeTracklets = leafFake->GetValue();
        double RecoMCvalidated = leafReco->GetValue(); 
        double TGenerated = leafTot->GetValue(); 
        double lambdaCut = leafCut->GetValue(); 
        double EntryNumber = leafNum->GetValue();

        double EffReco=RecoMCvalidated/TGenerated;
        double EffFake= FakeTracklets/TGenerated;

        ArrEffReco[i]=EffReco;
        ArrEffFake[i]= EffFake;
        ArrInvCut[i]=1/lambdaCut;
        ArrNum[i]=EntryNumber;
        
        //std::cout<<FakeTracklets<<"\n";
        //std::cout<<EffReco<<"\n"
	}


/* 

   for(int i=0; i<60; i++){
            std::cout<<"Arrnum : "<<ArrNum[i]<<"    EffReco :"<<ArrEffReco[i]<<std::endl;
        }


*/

  //sorting the vectors to plot a curve later
    

    std::vector<std::size_t> index_vec;

    for (std::size_t i = 0; i < (unsigned int)nEntries; ++i) { index_vec.push_back(i); }
/* 
      for(int i=0; i<nEntries; i++){
            std::cout<<"index_vec : "<<index_vec[i]<<std::endl;
        }
*/
    std::stable_sort( index_vec.begin(), index_vec.end(), [&](std::size_t a, std::size_t b) { return ArrInvCut[a] < ArrInvCut[b]; });
/* 
     for(int i=0; i<nEntries; i++){
            std::cout<<"sorted index_vec : "<<index_vec[i]<<std::endl;
        }
*/
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


    

    /*

      for(int i=0; i<nEntries; i++){
            std::cout<<"Sorted Arrnum : "<<SortedArrNum[i]<<std::endl;
        }

     */

      for(int i=0; i<60; i++){
            std::cout<<"Sorted Arrnum : "<<SortedArrNum[i]<<"    EffReco :"<<SortedArrEffReco[i]<<std::endl;
        }


    //this worked

    // Computing the mean for each Tan Lambda Cut value

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
    
            CountCurrentInvCut.push_back(0); //always use push back for the 1st element of the vector
            
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
        std:cout<<"Number of Entries :"<<c<<std::endl;
    }
    
    

     
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

    plot_quantiles(nEntries, SortedArrInvCut, SortedArrEffReco, SortedArrEffFake, SortedArrNum, CountCurrentInvCut,DistinctCuts);


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
    
    

    //std::cout<<"to draw : "<<NumDistinctCuts<<" "<<*DistinctCuts.data()<<"   "<<*MeanEffRecoCut.data()<<std::endl;
}



    void plot_quantiles(int nEntries, double * ArrInvCut, double *ArrEffReco, double *ArrEffFake,int * ArrNum, std::vector <int> CountCurrentInvCut,
     std::vector<double> DistinctCuts){
        int Num1; //minimum
        int Num2;
        int Num3;
        int Num4;
        int Num5; //maximum

        // the goal is to find the Number of all the quantiles for the 1st cut value and to then get all the values for this number
        // an easy way is to sort the efficiency array for the 1st cut value

        int FirstCutEntries = CountCurrentInvCut[0];
        int offset = FirstCutEntries/5;

        double FirstCutArrEffReco [FirstCutEntries ];
        int FirstCutArrNum [FirstCutEntries ];
        

        std::vector<std::size_t> index_vec;
        for (std::size_t i = 0; i < (unsigned int) FirstCutEntries ; ++i) { index_vec.push_back(i); }
        std::sort( index_vec.begin(), index_vec.end(), [&](std::size_t a, std::size_t b) { return ArrEffReco[a] < ArrEffReco[b]; });
        for (std::size_t i = 0; i < (unsigned int) FirstCutEntries ; ++i){
            FirstCutArrEffReco[i]=ArrEffReco[index_vec[i]];
            FirstCutArrNum[i]=ArrNum[index_vec[i]];
         }

        for(int i=0; i<FirstCutEntries; i++){
            std::cout<<"Arrnum : "<<FirstCutArrNum[i]<<"    EffReco :"<<FirstCutArrEffReco[i]<<std::endl;
        }


        Num1= FirstCutArrNum[0]; //Number of the entry with the lowset Efficiency with the first lambda cut value
        Num5= FirstCutArrNum[FirstCutEntries-1];
        Num2= FirstCutArrNum[offset] ;
        Num3= FirstCutArrNum[2*offset] ;
        Num4 = FirstCutArrNum[FirstCutEntries-1-offset];

        std::cout<<"offset :"<<offset<<std::endl;
        std::cout<<"Num 1 :"<<Num1<<std::endl;
        std::cout<<"Num 2 :"<<Num2<<std::endl;
        std::cout<<"Num 3 :"<<Num3<<std::endl;
        std::cout<<"Num 4 :"<<Num4<<std::endl;
        std::cout<<"Num 5 :"<<Num5<<std::endl;

        int NbCuts = DistinctCuts.size();
     
        std::sort(DistinctCuts.begin(), DistinctCuts.end()); //we sort it for the plots later

        for(int i=0; i<NbCuts; i++){
            std::cout<<DistinctCuts[i]<<std::endl;
        }

        connect_points(ArrEffReco, ArrNum, Num1,  NbCuts, CountCurrentInvCut,DistinctCuts);
        connect_points(ArrEffReco, ArrNum, Num2,  NbCuts, CountCurrentInvCut,DistinctCuts);
        connect_points(ArrEffReco, ArrNum, Num3,  NbCuts, CountCurrentInvCut,DistinctCuts);
        connect_points(ArrEffReco, ArrNum, Num4,  NbCuts, CountCurrentInvCut,DistinctCuts);
        connect_points(ArrEffReco, ArrNum, Num5,  NbCuts, CountCurrentInvCut,DistinctCuts);

        /* 
        double * Arr1  = fill_array(ArrEffReco, ArrNum, Num1,  NbCuts, CountCurrentInvCut);
        double * Arr2  = fill_array(ArrEffReco, ArrNum, Num2,  NbCuts, CountCurrentInvCut);
        double * Arr3  = fill_array(ArrEffReco, ArrNum, Num3,  NbCuts, CountCurrentInvCut);
        double * Arr4  = fill_array(ArrEffReco, ArrNum, Num4,  NbCuts, CountCurrentInvCut);
        double * Arr5  = fill_array(ArrEffReco, ArrNum, Num5,  NbCuts, CountCurrentInvCut);

        for(int i=0; i<NbCuts; i++){
            std::cout<<Arr2[i]<<std::endl;
        }

        TGraph * graph1= new TGraph(NbCuts, DistinctCuts.data(),Arr1);
        graph1->SetMarkerStyle(20);
        graph1-> SetMarkerColor(4);
        graph1->Draw("PC");

        TGraph * graph2= new TGraph(NbCuts, DistinctCuts.data(),Arr2);
        graph2->SetMarkerStyle(20);
        graph2-> SetMarkerColor(4);
        graph2->Draw("PC");

        delete Arr1;
        delete Arr2;
        delete Arr3;
        delete Arr4;
        delete Arr5;
*/
        //filling the arrays with values
        /* 
        int jOffset=0;
        for(int i=0; i<NbCuts; i++){ 
            for(int j=jOffset; j<jOffset+CounctCurrentInvCut[i]; j++){ //we only look at one tan lambda cut value
                switch(ArrNum[j]){ //we fill the arrays with the data corresponding to the same entry
                    case Num1 : Arr1[i]=ArrEffReco[j]; break;
                    case Num2 : Arr2[i]=ArrEffReco[j]; break;
                    case Num3 : Arr3[i]=ArrEffReco[j]; break;
                    case Num4 : Arr4[i]=ArrEffReco[j]; break;
                    case Num5 : Arr5[i]=ArrEffReco[j]; break;
                }
            }
            jOffset+=CounctCurrentInvCut[i];
      
        }*/
    }
    

    void connect_points(double *ArrEffReco, int * ArrNum, int Num, int NbCuts, std::vector <int>  CountCurrentInvCut, std::vector <double> DistinctCuts){
        double * Array = new double [NbCuts];
        int jOffset=0;
        
        for(int i=0; i<NbCuts; i++){
             //std::cout<<"J offset :"<<jOffset<<std::endl;
            for(int j=jOffset; j<jOffset+CountCurrentInvCut[i]; j++){ //we only look at one tan lambda cut value
                if(ArrNum[j]==Num){ //this is the entry related to this number
                    Array[i]=ArrEffReco[j];
                    break;
                }  
            }
            jOffset+=CountCurrentInvCut[i];
        }
        TGraph * graph1= new TGraph(NbCuts, DistinctCuts.data(),Array);
        graph1->SetMarkerStyle(20);
        graph1-> SetMarkerColor(4);
        graph1->Draw("PC");
        delete [] Array;
    }
    
    
#endif