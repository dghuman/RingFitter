// Just a macro that uses the PMT hit positions from the fiTQun modification until I can get it to work without fitqun to fit rings

#include <iostream>

void fitter(char *filename=NULL) {
  TFile* file;
  TNtuple* Ntuple;
  TH2D* ThetaPhi;   // The 2-Dim histogram that will keep track of the possible muon tracks that could produce the PMT hit observed at for a fixed position of interest. Will be remade for each position.

  Double_t pmtxpos,pmtypos,pmtzpos,pmttime;
  Double_t* chrnkv[3];
  int entries;

  const Double_t posgran = 25.; // The spatial cube dimension that will be used as the starting position of each event track check
  const Double_t thetaerror = 2.; // Error in the track angle reconstruction. This should be more dynamic based off of the position of interest and the pmt being hit. Const for now.
  const Double_t phierror = 2.;
  const Double_t pi = 3.14159265535897;
  const Double_t chrnkvangle = pi/4.; // Angle to track at which we assume that the cherenkov light is being emmitted.
  const Double_t xmax = 380.;
  const Double_t ymax = 380.; 
  const Double_t zmax = 525.;

  // Open the file
  if (filename != NULL) {
    file = new TFile(filename);
  } else {
    std::cout << "File could not be found. Exiting ..." << std::endl;
    return -1;
  }
  // Check that the file is open
  if (!f->IsOpen()) {
    std::cout << "Could not open " << filename << ". Exiting ... " << std::endl;
    return -1;
  }
  
  // Read in the Ntuple from the file
  Ntuple = (TNtuple*)TheFile.Get("hitpmts");
  Ntuple->SetBranchAddress("xpos",&pmtxpos);
  Ntuple->SetBranchAddress("ypos",&pmtypos);
  Ntuple->SetBranchAddress("zpos",&pmtzpos);
  Ntuple->SetBranchAddress("time",&pmttime);

  entries = (int) Ntuple->GetEntries();
  // Loop through the positions inside the detector and check the tracks that could have led to each pmt hit
  for (float xpos = -xmax; xpos < xmax; xpos+=posgran) {
    for (float ypos = -ymax; ypos < ymax; ypos+=posgran) {
      for (float zpos = -zmax; zpos < zmax; zpos+=posgran) {
	for (int i = 0; i < entries; i++) {
	  // Read in the Ntuple
	  Ntuple->GetEntry(i);   
	  // Compute the Cherenkov vector
	  chrnkv[0] = pmtxpos - xpos;
	  chrnkv[1] = pmtypos - ypos;
	  chrnkv[2] = pmtzpos - zpos;
	  
      }
    }
  }
  
}
