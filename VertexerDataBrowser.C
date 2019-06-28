#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <iostream>
#include <string>
#include <cstring>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphPolar.h"
#include <TH1F.h>
#include "TH2F.h"
#include "TLeaf.h"
#include "TNtuple.h"
#include "TPaveStats.h"
#include "TFrame.h"

void Drawing(double DefaultMin, bool compare, TH1F * hist, TH1F * hist_ref);


void VertexerDataBrowser(
    const std::string fileNameResults, /*, const std::string fileNameData*/
    const std::string fileNameCompareTo = "") {

  TFile vertexerData(fileNameResults.data());
  // Comparison block
  char compare{!fileNameCompareTo.empty()};
  TFile *vertexerComparisonData{NULL};

  if (compare) {
    vertexerComparisonData = TFile::Open(fileNameCompareTo.data());
  }


  //Creation of the Canvases

  auto canvasClusters = new TCanvas("Clusters", "Clusters data", 1000, 1000);
  auto canvasTracklets = new TCanvas("Tracklets", "Tracklets data", 1000, 1000);
  auto canvasLines = new TCanvas("Lines", "Lines data", 1000, 1000);
  auto canvasPolar = new TCanvas("Polar", "Polar", 1000, 1000);
  auto canvasDebugTanL = new TCanvas("dbTanL", "dbTanL", 1000, 1000);


  //Definition of the TNtuples

  TNtuple *lines = (TNtuple *)vertexerData.Get("Tracklets");
  TNtuple *comb01 = (TNtuple *)vertexerData.Get("comb01");
  TNtuple *comb12 = (TNtuple *)vertexerData.Get("comb12");
  TNtuple *phi01 = (TNtuple *)vertexerData.Get("clus_phi01");
  TNtuple *phi12 = (TNtuple *)vertexerData.Get("clus_phi12");
  TNtuple *deltaTanLambda = (TNtuple *)vertexerData.Get("dtl");
  TNtuple *centroids = (TNtuple *)vertexerData.Get("centroids");
  TNtuple *linesData = (TNtuple *)vertexerData.Get("ld");

  //Reference Tntuples

  TNtuple *linesRef{NULL};
  TNtuple *comb01Ref{NULL};
  TNtuple *comb12Ref{NULL};
  TNtuple *phi01Ref{NULL};
  TNtuple *phi12Ref{NULL};
  TNtuple *deltaTanLambdaRef{NULL};
  TNtuple *centroidsRef{NULL};
  TNtuple *linesDataRef{NULL};

  

  // Comparison data
  if (compare) {
    linesRef = (TNtuple *)vertexerComparisonData->Get("Tracklets");
    comb01Ref = (TNtuple *)vertexerComparisonData->Get("comb01");
    comb12Ref = (TNtuple *)vertexerComparisonData->Get("comb12");
    phi01Ref = (TNtuple *)vertexerComparisonData->Get("clus_phi01");
    phi12Ref = (TNtuple *)vertexerComparisonData->Get("clus_phi12");
    deltaTanLambdaRef = (TNtuple *)vertexerComparisonData->Get("dtl");
    centroidsRef = (TNtuple *)vertexerComparisonData->Get("centroids");
    linesDataRef = (TNtuple *)vertexerComparisonData->Get("ld");
  }

  //Definition of the histograms

  TH1F *histComb01Phi =
      new TH1F("hComb01phi", "Comb01 Phi angle; Phi (rad) ; Number of Tracklets", 200, -3.15f, 3.15f);
 
  TH1F *histComb12Phi =
      new TH1F("hComb12phi", "Comb12 Phi angle; Phi (rad) ; Number of Tracklets", 200, -3.15f, 3.15f);
  TH1F *histClus0Phi =
      new TH1F("hCl0phi", "Clusters Layer 0 Phi;Phi (rad); Number of Clusters", 200, 0.f, 6.3f);   
  TH1F *histClus1Phi =
      new TH1F("hCl1phi", "Clusters Layer 1 Phi;Phi (rad); Number of Clusters", 200, 0.f, 6.3f);
  TH1F *histClus2Phi =
      new TH1F("hCl2phi", "Clusters Layer 2 Phi;Phi (rad); Number of Clusters ", 200, 0.f, 6.3f);
  TH1F *histComb01TanLambda =
      new TH1F("hComb01tlambda", "Comb01 tanLambda; TanLambda; Number of Tracklets", 100, -65.f, 65.f);
  TH1F *histComb12TanLambda =
      new TH1F("hComb12tlambda", "Comb12 tanLambda; TanLambda ; Number of Tracklets", 100, -65.f, 65.f);
  TH1F *histDeltaTanLambda =
      new TH1F("hdtlambda", "DeltaTanLambda; Difference in TanLambda; Number of Tracklets", 100, 0, 0.03f);    
  TH2F *histDelta01Phi1 = new TH2F("hDelta01Phi1", "DeltaPhi 01 vs Phi1", 100,
                                   0.f, 6.29f, 200, -1.f, 1.f);
  TH2F *histDelta12Phi1 = new TH2F("hDelta12Phi1", "DeltaPhi 12 vs Phi1", 100,
                                   0.f, 6.29f, 200, -1.f, 1.f);
  TH2F *histPhi0VsPhi1 =
      new TH2F("hPhi0Phi1", "phi0 vs phi1", 200, 0.f, 6.29f, 200, 0.f, 6.29f);
  TH2F *histPhi1VsPhi2 =
      new TH2F("hPhi1Phi2", "phi1 vs phi2", 200, 0.f, 6.29f, 200, 0.f, 6.29f);
  TH1F *histPhi01 = new TH1F("hPhi01", "Phi0 - Phi1; Phi0 - Phi1 (rad); Number of Clusters", 200, -1.f, 1.f);
  TH1F *histPhi12 = new TH1F("hPhi12", "Phi1 - Phi2; Phi1 - Phi2 (rad); Number of Clusters", 200, -1.f, 1.f);
  TH1F *histOrigX = new TH1F("hOrigX", "Origin X; X (cm); Number of Lines ", 200, -3, 3);
  TH1F *histOrigY = new TH1F("hOrigY", "Origin Y; Y (cm); Number of Lines", 200, -3, 3);
  TH1F *histOrigZ = new TH1F("hOrigZ", "Origin Z; Z (cm); Number of Lines", 200, -25, 25);
  TH1F *histCosDir1 =
      new TH1F("hCosDir1", "Cosine director 1; Cosine ; Number of Lines", 100, -1.2f, 1.2f);
  TH1F *histCosDir2 =
      new TH1F("hCosDir2", "Cosine director 2; Cosine ; Number of Lines", 100, -1.2f, 1.2f);
  TH1F *histCosDir3 =
      new TH1F("hCosDir3", "Cosine director 3; Cosine ; Number of Lines", 100, -1.2f, 1.2f);
  TH2F *histClusTanLambda = new TH2F("hClusTanLambda", "Clusters on 3 Layers",
                                     300, -16.33f, 16.33f, 300, 0.f, 5.f);
  
  //Setting and defining minima

  histComb01Phi->SetMinimum(0);
  histComb12Phi->SetMinimum(0);
  histClus0Phi->SetMinimum(0);
  histClus1Phi->SetMinimum(0);
  histClus2Phi->SetMinimum(0);
  double MinLogLines = (linesData->GetEntries() < 1000) ? 0.1 : 10;
  double MinLogTracklets = (lines->GetEntries() < 1000) ? 0.1 : 1;

  TH1F *histComb01Phi_ref{NULL};
  TH1F *histComb12Phi_ref{NULL};
  TH1F *histComb01Phi_diff{NULL};
  TH1F *histComb12Phi_diff{NULL};

  TH1F *histClus0Phi_ref{NULL};
  TH1F *histClus1Phi_ref{NULL};
  TH1F *histClus2Phi_ref{NULL};

  TH1F *histClus0Phi_diff{NULL};
  TH1F *histClus1Phi_diff{NULL};
  TH1F *histClus2Phi_diff{NULL};

  TH1F *histComb01TanLambda_ref{NULL};
  TH1F *histComb12TanLambda_ref{NULL};
  TH1F *histDeltaTanLambda_ref{NULL};
  TH1F *histComb01TanLambda_diff{NULL};
  TH1F *histComb12TanLambda_diff{NULL};
  TH1F *histDeltaTanLambda_diff{NULL};

  TH1F *histPhi01_ref{NULL};
  TH1F *histPhi12_ref{NULL};
  TH1F *histPhi01_diff{NULL};
  TH1F *histPhi12_diff{NULL};

  TH2F *histDelta01Phi1_ref{NULL};
  TH2F *histDelta12Phi1_ref{NULL};

  TH2F *histPhi0VsPhi1_ref{NULL};
  TH2F *histPhi1VsPhi2_ref{NULL};

  TH1F *histOrigX_ref{NULL};
  TH1F *histOrigY_ref{NULL};
  TH1F *histOrigZ_ref{NULL};
  TH1F *histCosDir1_ref{NULL};
  TH1F *histCosDir2_ref{NULL};
  TH1F *histCosDir3_ref{NULL};

  TH1F *histOrigX_diff{NULL};
  TH1F *histOrigY_diff{NULL};
  TH1F *histOrigZ_diff{NULL};
  TH1F *histCosDir1_diff{NULL};
  TH1F *histCosDir2_diff{NULL};
  TH1F *histCosDir3_diff{NULL};

  TH2F *histClusTanLambda_ref{NULL};

 

  // Comparison data
  if (compare) {
    histComb01Phi_ref =
        new TH1F("hComb01phiRef", "Comb01 Phi angle", 200, -3.15f, 3.15f);
    histComb12Phi_ref =
        new TH1F("hComb12phiRef", "Comb12 Phi angle", 200, -3.15f, 3.15f);
    histClus0Phi_ref =
        new TH1F("hCl0phiRef", "Clusters Layer 0 Phi; Phi (rad); Count", 200, 0.f, 6.3f);
    histClus1Phi_ref =
        new TH1F("hCl1phiRef", "Clusters Layer 1 Phi;X;Y", 200, 0.f, 6.3f);
    histClus2Phi_ref =
        new TH1F("hCl2phiRef", "Clusters Layer 2 Phi;X;Y", 200, 0.f, 6.3f);
    histComb01TanLambda_ref =
        new TH1F("hComb01tlambdaRef", "Comb01 tanLambda;X;Y", 100, -65.f, 65.f);
    histComb12TanLambda_ref =
        new TH1F("hComb12tlambdaRef", "Comb12 tanLambda;X;Y", 100, -65.f, 65.f);
    histDeltaTanLambda_ref =
      new TH1F("hdtlambdaRef", "DeltaTanLambda;X;Y", 100, -0.03f, 0.03f);  
    histDelta01Phi1_ref = new TH2F("hDelta01Phi1Ref", "DeltaPhi 01 vs Phi1",
                                   100, 0.f, 6.29f, 200, -1.f, 1.f);
    histDelta12Phi1_ref = new TH2F("hDelta12Phi1Ref", "DeltaPhi 12 vs Phi1",
                                   100, 0.f, 6.29f, 200, -1.f, 1.f);
    histPhi0VsPhi1_ref = new TH2F("hPhi0Phi1Ref", "phi0 vs phi1", 200, 0.f,
                                  6.29f, 200, 0.f, 6.29f);
    histPhi1VsPhi2_ref = new TH2F("hPhi1Phi2Ref", "phi1 vs phi2", 200, 0.f,
                                  6.29f, 200, 0.f, 6.29f);
    histPhi01_ref = new TH1F("hPhi01Ref", "Phi0 - Phi1", 200, -1.f, 1.f);
    histPhi12_ref = new TH1F("hPhi12Ref", "Phi1 - Phi2", 200, -1.f, 1.f);
    histOrigX_ref = new TH1F("hOrigXRef", "Origin X", 200, -3, 3);
    histOrigY_ref = new TH1F("hOrigYRef", "Origin Y", 200, -3, 3);
    histOrigZ_ref = new TH1F("hOrigZRef", "Origin Z", 200, -25, 25);
    histCosDir1_ref =
        new TH1F("hCosDir1Ref", "Cosine director 1", 100, -1.2f, 1.2f);
    histCosDir2_ref =
        new TH1F("hCosDir2Ref", "Cosine director 2", 100, -1.2f, 1.2f);
    histCosDir3_ref =
        new TH1F("hCosDir3Ref", "Cosine director 3", 100, -1.2f, 1.2f);
    histClusTanLambda_ref =
        new TH2F("hClusTanLambdaRef", "Clusters on 3 Layers", 300, -16.33f,
                 16.33f, 300, 0.f, 5.f);
  }

  
  //Dividing the canvases

  canvasTracklets->Divide(3, 2);
  canvasClusters->Divide(2, 3);
  canvasLines->Divide(3, 2);
  canvasPolar->Divide(2, 2);

  //Projection of the TNtuples on the histograms

  comb01->Project("hComb01phi", "phi"); 
  comb12->Project("hComb12phi", "phi");
  comb01->Project("hComb01tlambda", "tanLambda");
  comb12->Project("hComb12tlambda", "tanLambda");
  deltaTanLambda->Project("hdtLambda", "deltatanlambda");
  
  phi01->Project("hDelta01Phi1", "phi0-phi1:phi1");
  phi12->Project("hDelta12Phi1", "phi1-phi2:phi1");
  phi01->Project("hPhi01", "phi0-phi1");
  phi12->Project("hPhi12", "phi1-phi2");
  phi01->Project("hPhi0Phi1", "phi0:phi1");
  phi12->Project("hPhi1Phi2", "phi1:phi2");

  phi01->Project("hCl0phi", "phi0");
  phi01->Project("hCl1phi", "phi1");
  phi12->Project("hCl2phi", "phi2");

  lines->Project("hOrigX", "oX");
  lines->Project("hOrigY", "oY");
  lines->Project("hOrigZ", "oZ");
  lines->Project("hCosDir1", "c1");
  lines->Project("hCosDir2", "c2");
  lines->Project("hCosDir3", "c3");



  // Comparison data
  if (compare) {
    comb01Ref->Project("hComb01phiRef", "phi");
    comb12Ref->Project("hComb12phiRef", "phi");
    comb01Ref->Project("hComb01tlambdaRef", "tanLambda");
    comb12Ref->Project("hComb12tlambdaRef", "tanLambda");
    deltaTanLambdaRef->Project("hdtLambdaRef", "deltatanlambda");

    phi01Ref->Project("hDelta01Phi1Ref", "phi0-phi1:phi1");
    phi12Ref->Project("hDelta12Phi1Ref", "phi1-phi2:phi1");
    phi01Ref->Project("hPhi01Ref", "phi0-phi1");
    phi12Ref->Project("hPhi12Ref", "phi1-phi2");
    phi01Ref->Project("hPhi0Phi1Ref", "phi0:phi1");
    phi12Ref->Project("hPhi1Phi2Ref", "phi1:phi2");

    phi01Ref->Project("hCl0phiRef", "phi0");
    phi01Ref->Project("hCl1phiRef", "phi1");
    phi12Ref->Project("hCl2phiRef", "phi2");

    linesRef->Project("hOrigXRef", "oX");
    linesRef->Project("hOrigYRef", "oY");
    linesRef->Project("hOrigZRef", "oZ");
    linesRef->Project("hCosDir1Ref", "c1");
    linesRef->Project("hCosDir2Ref", "c2");
    linesRef->Project("hCosDir3Ref", "c3");
  }

  
  //Filling the special histograms

  for (auto i{0}; i < deltaTanLambda->GetEntries(); ++i) {
    deltaTanLambda->GetEntry(i);
    // std::cout<<deltaTanLambda->GetLeaf("deltatanlambda")->GetValue()<<std::endl;
    if (deltaTanLambda->GetLeaf("deltatanlambda")->GetValue() > 0.0f) {
      histClusTanLambda->Fill(deltaTanLambda->GetLeaf("c0z")->GetValue(),
                              deltaTanLambda->GetLeaf("c0r")->GetValue());
      histClusTanLambda->Fill(deltaTanLambda->GetLeaf("c1z")->GetValue(),
                              deltaTanLambda->GetLeaf("c1r")->GetValue());
      histClusTanLambda->Fill(deltaTanLambda->GetLeaf("c2z")->GetValue(),
                              deltaTanLambda->GetLeaf("c2r")->GetValue());
      // break;
    }
    histDeltaTanLambda->Fill(deltaTanLambda->GetLeaf("deltatanlambda")->GetValue());
  }

  // Comparison data
  if (compare) {
    for (auto i{0}; i < deltaTanLambdaRef->GetEntries(); ++i) {
      deltaTanLambdaRef->GetEntry(i);
      // std::cout<<deltaTanLambda->GetLeaf("deltatanlambda")->GetValue()<<std::endl;
      if (deltaTanLambdaRef->GetLeaf("deltatanlambda")->GetValue() > 0.0f) {
        histClusTanLambda_ref->Fill(
            deltaTanLambdaRef->GetLeaf("c0z")->GetValue(),
            deltaTanLambdaRef->GetLeaf("c0r")->GetValue());
        histClusTanLambda_ref->Fill(
            deltaTanLambdaRef->GetLeaf("c1z")->GetValue(),
            deltaTanLambdaRef->GetLeaf("c1r")->GetValue());
        histClusTanLambda_ref->Fill(
            deltaTanLambdaRef->GetLeaf("c2z")->GetValue(),
            deltaTanLambdaRef->GetLeaf("c2r")->GetValue());
        // break;
      }
      histDeltaTanLambda_ref->Fill(deltaTanLambdaRef->GetLeaf("deltatanlambda")->GetValue());
    }
  }

  

  //Polar Canvas

  //Definition

  const int npoints = 200;

  double r01[npoints];
  double phiArr01[npoints];

  double r12[npoints];
  double phiArr12[npoints];

  for (int ipt = 1; ipt < npoints + 1; ipt++) {
    r01[ipt - 1] = histComb01Phi->GetXaxis()->GetBinCenter(ipt);
    phiArr01[ipt - 1] = histComb01Phi->GetBinContent(ipt);
    r12[ipt - 1] = histComb12Phi->GetXaxis()->GetBinCenter(ipt);
    phiArr12[ipt - 1] = histComb12Phi->GetBinContent(ipt);
  }

  // Comparison data

  double r01Ref[npoints];
  double phiArr01Ref[npoints];

  double r12Ref[npoints];
  double phiArr12Ref[npoints];

  if (compare) {
    for (int ipt = 1; ipt < npoints + 1; ipt++) {
      r01Ref[ipt - 1] = histComb01Phi_ref->GetXaxis()->GetBinCenter(ipt);
      phiArr01Ref[ipt - 1] = histComb01Phi_ref->GetBinContent(ipt);
      r12Ref[ipt - 1] = histComb12Phi_ref->GetXaxis()->GetBinCenter(ipt);
      phiArr12Ref[ipt - 1] = histComb12Phi_ref->GetBinContent(ipt);
    }
  }

  TGraphPolar *grP01 = new TGraphPolar(npoints, r01, phiArr01);
  TGraphPolar *grP12 = new TGraphPolar(npoints, r12, phiArr12);
  grP01->SetLineWidth(2);
  grP01->SetLineColor(kBlue + 2);
  grP01->SetMinRadial(0);
  grP12->SetLineWidth(2);
  grP12->SetLineColor(kRed - 3);
  grP12->SetMinRadial(0);
  canvasPolar->cd(4);
  grP01->Draw("LEP");
  grP12->Draw("LEP");
  canvasPolar->Update();
  grP01->GetPolargram()->SetRadialLabelSize(0.025);
  grP12->GetPolargram()->SetRadialLabelSize(0.025);

  histClus0Phi->SetDirectory(0);
  histClus0Phi->SetLineColor(kOrange + 8);
  histClus1Phi->SetDirectory(0);
  histClus1Phi->SetLineColor(kOrange + 8);
  histClus2Phi->SetDirectory(0);
  histClus2Phi->SetLineColor(kOrange + 8);

  // Comparison data
  if (compare) {
    TGraphPolar *grP01Ref = new TGraphPolar(npoints, r01Ref, phiArr01Ref);
    TGraphPolar *grP12Ref = new TGraphPolar(npoints, r12Ref, phiArr12Ref);
    grP01Ref->SetLineWidth(2);
    grP01Ref->SetLineColor(kAzure + 6);
    grP01Ref->SetMinRadial(0);
    grP12Ref->SetLineWidth(2);
    grP12Ref->SetLineColor(kMagenta + 1);
    grP12Ref->SetMinRadial(0);
    canvasPolar->cd(4);
    grP01Ref->Draw("LEP");
    grP12Ref->Draw("LEP");
    canvasPolar->Update();
    grP01Ref->GetPolargram()->SetRadialLabelSize(0.025);
    grP12Ref->GetPolargram()->SetRadialLabelSize(0.025);
    
  }


  //Drawing 

  canvasPolar->cd(1);
  Drawing(0, compare, histClus0Phi, histClus0Phi_ref);
  canvasPolar->cd(2);
  Drawing(0, compare, histClus1Phi, histClus1Phi_ref);
  canvasPolar->cd(3);
  Drawing(0, compare, histClus2Phi, histClus2Phi_ref);


  // Tracklets

  canvasTracklets->cd(1);
  Drawing(0, compare, histComb01Phi, histComb01Phi_ref);
  canvasTracklets->cd(2);
  Drawing(0, compare, histComb12Phi, histComb12Phi_ref);
  canvasTracklets->cd(4)->SetLogy();
  Drawing(MinLogTracklets, compare, histComb01TanLambda, histComb01TanLambda_ref);
  canvasTracklets->cd(5)->SetLogy();
  Drawing(MinLogTracklets, compare, histComb12TanLambda, histComb12TanLambda_ref);
  canvasTracklets->cd(6)->SetLogy();
  Drawing(MinLogTracklets, compare, histDeltaTanLambda, histDeltaTanLambda_ref);


  // Phi Clusters block


  histDelta01Phi1->SetDirectory(0);
  canvasClusters->cd(1);
  histDelta01Phi1->Draw("colz");
  histDelta12Phi1->SetDirectory(0);
  canvasClusters->cd(2);
  histDelta12Phi1->Draw("colz");
  histPhi0VsPhi1->SetDirectory(0);
  canvasClusters->cd(3);
  histPhi0VsPhi1->Draw("colz");
  histPhi1VsPhi2->SetDirectory(0);
  canvasClusters->cd(4);
  histPhi1VsPhi2->Draw("colz");

  canvasClusters->cd(5)->SetLogy();
  Drawing(100, compare, histPhi01, histPhi01_ref);
  canvasClusters->cd(6)->SetLogy();
  Drawing(100, compare, histPhi12, histPhi12_ref);

  // Lines block

  canvasLines->cd(1)->SetLogy();
  Drawing(MinLogLines, compare, histOrigX, histOrigX_ref);
  canvasLines->cd(2)->SetLogy();
  Drawing(MinLogLines, compare, histOrigY, histOrigY_ref);
  canvasLines->cd(3)->SetLogy();
  Drawing(MinLogLines, compare, histOrigZ, histOrigZ_ref);
  canvasLines->cd(4)->SetLogy();
  Drawing(MinLogLines, compare, histCosDir1, histCosDir1_ref);
  canvasLines->cd(5)->SetLogy();
  Drawing(MinLogLines, compare, histCosDir2, histCosDir2_ref);
  canvasLines->cd(6)->SetLogy();
  Drawing(MinLogLines, compare, histCosDir3, histCosDir3_ref);

  // Debug clusters block
  histClusTanLambda->SetDirectory(0);
  canvasDebugTanL->cd();
  histClusTanLambda->Draw("colz");

}


void Drawing(double DefaultMin, bool compare, TH1F * hist, TH1F * hist_ref ){

TH1F * hist_diff;

  if(compare){
    string HistDiffName = (string)(hist->GetName());
    char diff []= "Diff";
    HistDiffName = HistDiffName + diff;
    hist_diff = (TH1F *)hist_ref->Clone(HistDiffName.c_str());
    hist_diff->Add(hist, -1);

    double Max= 1.1* std::max(hist_ref->GetMaximum(), hist->GetMaximum());
    double Min= DefaultMin;
    hist->GetYaxis()->SetRangeUser(Min, Max);
    hist_ref->GetYaxis()->SetRangeUser(Min, Max);
    hist_diff->GetYaxis()->SetRangeUser(Min, Max);

    hist_ref->SetDirectory(0);
    hist_ref->SetLineColor(kBlue + 2);

    hist_diff->SetDirectory(0);
    hist_diff->SetFillStyle(3001);
    hist_diff->SetFillColor(kBlue);

  }else{
    hist->SetMinimum(DefaultMin);
  }

  hist->SetLineColor(kOrange + 8);
  hist->SetDirectory(0);
  hist->Draw();
  
  if(compare){
    hist_ref->Draw("sames");
    hist_diff->Draw("sames");
  }
}



#endif
