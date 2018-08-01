// Just a macro that uses the PMT hit positions from the fiTQun modification until I can get it to work without fitqun to fit rings

#include <iostream>
#include <math.h>

void fitter(char *filename=NULL) {
  TFile* file;
  TNtuple* Ntuple;
  TH2D* ThetaPhi;   // The 2-Dim histogram that will keep track of the possible muon tracks that could produce the PMT hit observed at for a fixed position of interest. Will be remade for each position.

  Double_t pmtxpos,pmtypos,pmtzpos,pmttime,length;
  Double_t* chrnkv[3];
  Double_t* polarchrknv[2]; // (theta,phi) where theta is the azimuthal (0 - 2pi)

  // Make the rotation matrix
  Double_t** rotmat = new Double_t*[3];
  rotmat[0] = new Double_t[9];

  // Actually initialize the other coloumns of the matrix
  for (int i = 1; i < 3; i++) {
    rotmat[i] = rotmat[0] + i*3;
  }

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
	  // Normalize chrnkv vector
	  length = sqrt(chrnkv[0]*chrnkv[0] + chrnkv[1]*chrnkv[1] + chrnkv[2]*chrnkv[2]);
	  chrnkv[0] = chrnkv[0]/length;
	  chrnkv[1] = chrnkv[1]/length;
	  chrnkv[2] = chrnkv[2]/length;
	  // Find the Spherical coordinates of this unit vector
	  polarchrnkv[0] = atan(chrnkv[1]/chrnkv[0]);
	  polarchrnkv[1] = acos(chrnkv[2]);
	  // Now to find the set of unit track vectors with a parameterization of a circle chrnkvangle off of the muon track. To do this, we take the usual circle parameteriztion in R^2 and push it in the z direction by 1 unit, then apply a rotation to the parameterization so that the central vector (0,0,1) points along our cherenkov path. The rotation matrices are uniquely defined by the theta and phi found earlier.
	  // First we make our rotation matrix, which was defined by a rotation about the y-axis through the angle phi and then a rotation about the z-axis by an angle theta.
	  rotmat[0][0] = cos(polarchrnkv[0])*cos(polarchrnkv[1]);
	  rotmat[0][1] = -sin(polarchrnkv[0]);
	  rotmat[0][2] = cos(polarchrnkv[0])*sin(polarchrnkv[1]);
	  rotmat[1][0] = sin(polarchrnkv[0])*cos(polarchrnkv[1]);
	  rotmat[1][1] = cos(polarchrnkv[0]);
	  rotmat[1][2] = sin(polarchrnkv[0])*sin(polarchrnkv[1]);
	  rotmat[2][0] = -sin(polarchrnkv[1]);
	  rotmat[2][1] = 0;
	  rotmat[2][2] = cos(polarchrnkv[1]);
	  // The idea here is to have a vector on our circle, say v, and find what it looks like in this rotated frame so that it is a predicted muon track. The circle is defined by using the predicted cherenkov path
	  
	}
      }
    }
  }
  // Clean up
  delete [] rotmat[0];
  delete [] rotmat;
}
