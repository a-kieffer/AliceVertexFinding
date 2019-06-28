#include <iostream>
#include <TFile.h>
#include "TSystem.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "TH1F.h"

//the aim of this macro is to read data from a tree and to display the data as an histogram

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

void my_macro()
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
	int counter{0};
	for (unsigned int i{0}; i < T->GetEntries(); ++i)
	{
		T->GetEntry(i);
		// std::cout << "vertices vector size: " << verticesITS_vector->size() << std::endl;
		for (auto &vertex : *verticesITS_vector)
		{
			std::cout << counter <<" vertex -> X: " << vertex.getX() << std::endl;
			++counter;
		}
	}

	// //T->StartViewer(); // it works
	// //T->Print();
	// TH1F *hX= new TH1F("hX","Histogram for X coordinates",100,-0.2,0.5);
	//
	// TTreeReader Treader("T", f);
	//
	// TTreeReaderValue<Float_t> Px(Treader, "ITSVertices.mPos.fCoordinates.fX");
	//
	//
	// cout <<"state:"<<Treader.Next()<<endl; //false
	//
	// while(Treader.Next()){
	// 	hX->Fill(*Px);
	// }
	//
	//
	//
	// hX->Draw();
}
