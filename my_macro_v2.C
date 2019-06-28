#include <iostream>
#include <TFile.h>
#include "TSystem.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"

//the aim of this macro is to read data from a tree and to display the data as an histogram

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

void my_macro_v2()
{
	//opening the file
	TFile *f = new TFile("vertexer_serial_data.root");
	//f->ShowStreamerInfo(); //seems to work

	//creating and initialazing the tree
	TTree *T;
	f->GetObject("o2sim", T);

	//creating a pointer to a vector of the vertices
	// every element of the vector is 0,1 or serveral vertices
	std::vector<Vertex> *verticesITS_vector = nullptr;

	//check if the Tree works
	if (!T->GetBranch("ITSVertices"))
	{
		std::cout << "data not found!\n";
	}

	//creating a pointer to a TBranch which points to the branch we need
	TBranch *vertices = T->GetBranch("ITSVertices");
	//we set the adress of the branch to the adress of our vector to store the data from the branch there
	vertices->SetAddress(&verticesITS_vector);

	// T->Branch("ITSVertices", &verticesITS_vector); this syntax is used when writing a Tree

	std::cout << "This Tree has " << T->GetEntries() << " entries\n";
	//int counter{0};

	int nEntries=T->GetEntries();

    //creating histograms
    TH1F *hX= new TH1F("hX","Histogram for X coordinates; X coordinates; Count",100,-0.2,0.5);
    TH1F *hY= new TH1F("hY","Histogram for Y coordinates; Y coordinates; Count",100,-0.2,0.4);
    TH1F *hZ= new TH1F("hZ","Histogram for Z coordinates; Z coordinates; Count",100,-30,30);

	for (int i{0}; i < nEntries; ++i)
	{
		// puts 0,1 or several vertices into verticesITS_vector
		T->GetEntry(i);
		// std::cout << "vertices vector size: " << verticesITS_vector->size() << std::endl;

        // for each vertex in the vector pointed by verticesITS_vector
        // "vertex" is a reference to the vertices in the vector
		for (auto &vertex : *verticesITS_vector)
		{
			//std::cout << counter <<" vertex -> X: " << vertex.getX() << std::endl;
			//++counter;
            hX->Fill(vertex.getX());
            hY->Fill(vertex.getY());
            hZ->Fill(vertex.getZ());
		}
	}


	TCanvas *c1= new TCanvas ("c", "Histograms", 1600, 600); // doesn't work without a pointer
    c1->Divide(3,1);

    c1->cd(1);

    hX->SetLineColor(kRed);
    hX->SetFillColor(kRed-7);
	hX->Draw();

    c1->cd(2);

    hY->SetLineColor(kCyan+2);
    hY->SetFillColor(kCyan);
    hY->Draw();

    c1->cd(3);

    hZ->SetLineColor(kGreen+3);
    hZ->SetFillColor(kGreen);
    hZ->Draw();
}
