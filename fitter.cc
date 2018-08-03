// Just a macro that uses the PMT hit positions from the fiTQun modification until I can get it to work without fitqun to fit rings

#include <iostream>
#include <math.h>

void fitter(char *filename=NULL, char *outfile=NULL) {
  TFile* file;
  TNtuple* Ntuple;
  TFile* NewFile;
  TNtuple* Track_Pos = new TNtuple("Track_Pos", "Predicted Muon Track and start Position", "xpos:ypos:zpos:xvec:yvec:zvec:nhits"); // An Ntuple to store the max hit vec, position and the number of corresponding hits.
  TH2D* ThetaPhi;   // The 2-Dim histogram that will keep track of the possible muon tracks that could produce the PMT hit observed at for a fixed position of interest. Will be remade for each position.

  Int_t counter = 0;
  Int_t maxbin,btheta,bphi,z,maxeve;
  Float_t pmtxpos,pmtypos,pmtzpos,pmttime,length;
  Double_t xfin = 0;
  Double_t yfin = 0;
  Double_t zfin = 0;
  Double_t vxmax = 0;
  Double_t vymax = 0;
  Double_t vzmax = 0;
  Double_t nmax = 0;; // For storing the predicted position and vector
  Double_t chrnkv[3];
  Double_t polarchrnkv[2]; // (theta,phi) where theta is the azimuthal (0 - 2pi)
  Double_t tempvec[3]; 
  Double_t vec[3]; 

  // Make the rotation matrix
  Double_t** rotmat = new Double_t*[3];
  rotmat[0] = new Double_t[9];

  // Actually initialize the other coloumns of the matrix
  for (int i = 1; i < 3; i++) {
    rotmat[i] = rotmat[0] + i*3;
  }

  int entries;

  const Double_t posgran = 20.; // The spatial cube dimension that will be used as the starting position of each event track check
  const Double_t pi = 3.14159265535897;
  const Double_t thetaerror = 10.*pi/180.; // Error in the track angle reconstruction. This should be more dynamic based off of the position of interest and the pmt being hit. Const for now.
  const Double_t phierror = 10.*pi/180.;
  const Double_t chrnkvangle = 43.*pi/180..; // Angle to track at which we assume that the cherenkov light is being emmitted.
  const Double_t radius = cos(chrnkvangle); // The corresponding radius and height of the ring that would be produced on the unitsphere given a cherenkov path along the z-axis
  const Double_t height = sin(chrnkvangle);
  const Double_t xmax = 380.;
  const Double_t ymax = 380.; 
  const Double_t zmax = 525.;
  const Int_t num = 100; // Number of points to take along the cone of muon tracks

  // Make a file to write to
  if (outfile != NULL) {
    NewFile = new TFile(outfile, "RECREATE");
    std::cout << "<OUTPUT> Will be writing the output to " << outfile << "." << std::endl;
  } else {
    NewFile = new TFile("output_fitter.root", "RECREATE");
    std::cout << "<OUTPUT> Will be writing the output to output_fitter.root." << std::endl;
  }

  // Check that the new file is open
  if (!NewFile->IsOpen()) { 
    std::cout << "<ERROR> Could not open the output file. Exiting." << std::endl;
    return -1;
  }

  // Open the file
  if (filename != NULL) {
    file = new TFile(filename);
  } else {
    std::cout << "<ERROR> File could not be found. Exiting ..." << std::endl;
    return -2;
  }
  // Check that the file is open
  if (!file->IsOpen()) {
    std::cout << "<ERROR> Could not open " << filename << ". Exiting ... " << std::endl;
    return -3;
  }
  
  // Read in the Ntuple from the file
  Ntuple = (TNtuple*)file->Get("hitpmts");
  Ntuple->SetBranchAddress("xpos",&pmtxpos);
  Ntuple->SetBranchAddress("ypos",&pmtypos);
  Ntuple->SetBranchAddress("zpos",&pmtzpos);
  Ntuple->SetBranchAddress("time",&pmttime);

  entries = (int) Ntuple->GetEntries();
  // Loop through the positions inside the detector and check the tracks that could have led to each pmt hit
  for (float xpos = -xmax; xpos < xmax; xpos+=posgran) {
    std::cout << "<STATUS> " << ((xpos + xmax)*100)/((xmax*2)) << " % complete" << std::endl; // Just a status check
    for (float ypos = -ymax; ypos < ymax; ypos+=posgran) {
      for (float zpos = -zmax; zpos < zmax; zpos+=posgran) {
	// Make a histogram for this position
	ThetaPhi = new TH2D(Form("thetaphi_%d",counter), "Theta Phi Parameter Space", ((int) 2*pi/thetaerror) + 1, 0.0, 2*pi, ((int) pi/phierror), 0.0, pi);
	file->cd();
	for (int entry = 0; entry < entries; entry++) {
	  // Read in the Ntuple for the first PMT position
	  Ntuple->GetEntry(entry);   
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
	  for (double t = 0; t < 2*pi; t+=(2*pi/((double)num))) {
	    tempvec[0] = radius*cos(t);
	    tempvec[1] = radius*sin(t);
	    tempvec[2] = height;
	    // Now we rotate the vector by multiplying it through the rotation matrix Mv = v_new, which is just matrix vector multiplication
	    for (int i = 0; i < 3; i++) {
	      for (int j = 0; j < 3; j++) {
		vec[i] = rotmat[i][j]*tempvec[j];
	      }
	    }
	    // Convert the vector to a spherical coordinate and feed it into the histogram
	    ThetaPhi->Fill(atan(vec[1]/vec[0]),acos(vec[2]));
	  }
	}
	counter++;
	maxbin = ThetaPhi->GetMaximumBin(); // This should be the most number of PMTs hit
	ThetaPhi->GetBinXYZ(maxbin,btheta,bphi,z);
	//std::cout << "Estimated theta: " << btheta*thetaerror << " phi: " << bphi*phierror << std::endl;
	Track_Pos->Fill((Float_t)xpos,(Float_t)ypos,(Float_t)zpos,(Float_t)cos(btheta*thetaerror),(Float_t)sin(btheta*thetaerror),(Float_t)cos(bphi*phierror), (Float_t)maxbin);
	// Keep track of the position and vector with the most matched PMTs
	if (maxbin > nmax) {
	  xfin = xpos;
	  yfin = ypos;
	  zfin = zpos;
	  vxmax = cos(btheta*thetaerror);
	  vymax = sin(btheta*thetaerror);
	  vzmax = cos(bphi*phierror);
	  nmax = maxbin;
	  maxeve = counter;
	}
	// This is if I feel it necessary to look at all of the plots produced. Will be VERY memory intensive.
	if (xpos == 0 && ypos == 0 && zpos == 0) {ThetaPhi->Write();}
	//NewFile->cd();
	//ThetaPhi->Write();  
	ThetaPhi->Reset();
      }
    }
  }
  Track_Pos->Write();
  std::cout << "<RESULT> The predicted pos of the muon is [x,y,z] = [" << xfin << "," << yfin << "," << zfin << "] and the predicted track is along [vx,vy,vz] = [" << vxmax << "," << vymax << "," << vzmax << "] with " << nmax << " number of hit pmts (probably). This is event number " << maxeve << ". Done!" << std::endl;
  // Clean up
  delete [] rotmat[0];
  delete [] rotmat;
}
