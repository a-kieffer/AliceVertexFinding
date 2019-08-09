# AliceVertexFinding
This repository is used as a backup to store the code I am developing while working on algorithms for Vertex finding in the ALICE experiment at CERN


The important macros are :

  * trackleter_lambda_cut.C
  * trackleter_phi_cut.C
  * plot_efficiency_lambda_cut.C
  * plot_efficiency_phi_cut.C
  * plot_efficiency_Pt.C
  
The two "trackleter" macros run the vertexer to reconstruct and/or select tracklets and store the information in a TNTuple in a file. The associated "plot" macros read the information from the file and produce plots.
The last macro creates and plots the data without storing it in a file. 
