// Just a macro that uses the PMT hit positions from the fiTQun modification until I can get it to work without fitqun to fit rings

#include <iostream>
#include <math.h>
#include <stdlib>
#include <vector>

void bruteforce(char *filename=NULL, char *outfile=NULL,bool debug=0) {
  TFile* file;
  TNtuple* Ntuple;
  TFile* NewFile;
  TNtuple* Track_Pos = new TNtuple("Track_Pos", "Predicted Muon Track and start Position", "xpos:ypos:zpos:xvec:yvec:zvec:nhits"); // An Ntuple to store the max hit vec, position and the number of corresponding hits.
  TH2D* ThetaPhi;   // The 2-Dim histogram that will keep track of the possible muon tracks that could produce the PMT hit observed at for a fixed position of interest. Will be remade for each position.

  Int_t counter = 0;
  Int_t maxbin,btheta,bphi,z,maxeve;
  Float_t pmtxpos,pmtypos,pmtzpos,pmttime;
  Double_t length,innerp;
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

  int entries;
  int bincounter;

  const Double_t posgran = 95.; // The spatial cube dimension that will be used as the starting position of each event track check
  const Double_t pi = 3.14159265535897;
  const Double_t thetaerror = 2.*pi/180.; // Error in the track angle reconstruction. This should be more dynamic based off of the position of interest and the pmt being hit. Const for now.
  const Double_t phierror = 2.*pi/180.;
  const Double_t chrnkvangle = 43.*pi/180..; // Angle to track at which we assume that the cherenkov light is being emmitted.
  const Double_t radius = sin(chrnkvangle); // The corresponding radius and height of the ring that would be produced on the unitsphere given a cherenkov path along the z-axis
  const Double_t height = cos(chrnkvangle);
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

  // Not sure how fast it is to read from the ntuple, so will be reading the ntuple into a vector to make things faster
  vector< vector<Float_t> > pmt;
  vector<Float_t> temppmt;


  for (int entry = 0; entry < entries; entry++) {
    temppmt.clear();
    Ntuple->GetEntry(entry);
    temppmt.push_back(pmtxpos);
    temppmt.push_back(pmtypos);
    temppmt.push_back(pmtzpos);
    temppmt.push_back(pmttime);
    pmt.push_back(temppmt);
  }

  // Switch to the new file
  NewFile->cd();
  // Loop through the positions inside the detector and check the tracks that could have led to each pmt hit
  for (Double_t xpos = -xmax; xpos < xmax; xpos+=posgran) {
    std::cout << "<STATUS> " << ((xpos + xmax)*100)/((xmax*2)) << " % complete" << std::endl; // Just a status check
    for (Double_t ypos = -ymax; ypos < ymax; ypos+=posgran) {
      for (Double_t zpos = -zmax; zpos < zmax; zpos+=posgran) {
	// Make a histogram for this position
	ThetaPhi = new TH2D(Form("thetaphi_%d",counter), "Theta Phi Parameter Space", ((int) 2*pi/thetaerror) + 1, 0.0 - thetaerror/2., 2*pi + thetaerror/2., ((int) pi/phierror) + 1, 0.0 - phierror/2., pi + phierror/2.);
	// ---------------------------- DEBUG --------------------
	if (debug == 1) {
	  std::cout << "<DEBUG> Entering DEBUG mode. Will check pos [x,y,z] = [0,0,0] only." << std::endl;
	  xpos = 0;
	  ypos = 0;
	  zpos = 0;
	}
	// Loop through all of the bins in this histogram and check the number of PMTs that match this possible track unit vector
	for (int tbin = 0; tbin < (((int) 2*pi/thetaerror) + 1); tbin++) {
	  for (int pbin = 0; pbin < (((int) pi/phierror) + 1); pbin++) {
	    // Convert to R3 vector to take inner product and check if it matches 
	    tempvec[0] = sin(pbin*phierror)*cos(tbin*thetaerror);
	    tempvec[1] = sin(pbin*phierror)*sin(tbin*thetaerror);
	    tempvec[2] = cos(pbin*phierror);
	    if (fabs(sqrt(tempvec[0]*tempvec[0] + tempvec[0]*tempvec[0] + tempvec[0]*tempvec[0]) - 1.) > 0.05) { 
	      std::cout << "<WARNING> Check computation for [theta,phi] = [" << tbin*thetaerror << "," << pbin*phierror << "]." << std::endl;
	    }
	    bincounter = 0;
	    for (int entry = 0; entry < entries; entry++) {
	      innerp = 0;
	      chrnkv[0] = pmt[entry][0] - xpos;
	      chrnkv[1] = pmt[entry][1] - ypos;
	      chrnkv[2] = pmt[entry][2] - zpos;	
	      // normalize for a simpler innerproduct
	      length = sqrt(chrnkv[0]*chrnkv[0] + chrnkv[1]*chrnkv[1] + chrnkv[2]*chrnkv[2]);
	      chrnkv[0] = chrnkv[0]/length;
	      chrnkv[1] = chrnkv[1]/length;
	      chrnkv[2] = chrnkv[2]/length;
	      // Check if the inner product matches what we would expect for chrnkv
	      for (int l = 0; l < 3; l++) {
		innerp += chrnkv[l]*tempvec[l];
	      }
	      if (fabs(innerp - height) < cos(phierror)) { bincounter++; } // use height for the expected innerproduct value since for unit vectors we would expect cos(43) to be what we want, which is exactly the height when working with a unit vector.
	    }
	    ThetaPhi->SetBinContent(tbin,pbin,(Double_t)bincounter);
	    // Keep track of the position and vector with the most matched PMTs
	    if (bincounter > nmax) {
	      xfin = xpos;
	      yfin = ypos;
	      zfin = zpos;
	      vxmax = tempvec[0];
	      vymax = tempvec[1];
	      vzmax = tempvec[2];
	      nmax = bincounter;
	      maxeve = counter;
	    }
	  }
	}
	counter++;
	maxbin = ThetaPhi->GetMaximumBin(); // This should be the most number of PMTs hit
	ThetaPhi->GetBinXYZ(maxbin,btheta,bphi,z);
	//std::cout << "Estimated theta: " << btheta*thetaerror << " phi: " << bphi*phierror << std::endl;
	Track_Pos->Fill((Float_t)xpos,(Float_t)ypos,(Float_t)zpos,(Float_t)sin(bphi*phierror)*(Float_t)cos(btheta*thetaerror),(Float_t)sin(bphi*phierror)*(Float_t)sin(btheta*thetaerror),(Float_t)cos(bphi*phierror),(Float_t)maxbin);
	// This is if I feel it necessary to look at all of the plots produced. Will be VERY memory intensive.
	ThetaPhi->Write();  
	ThetaPhi->Reset();
	// ---------------------------------- DEBUG ---------------------------------
	if (debug == 1) { break;}
	// --------------------------------------------------------------------------
      }
      // ---------------------------------- DEBUG ---------------------------------
      if (debug == 1) { break;}
      // --------------------------------------------------------------------------
    }
    // ---------------------------------- DEBUG ---------------------------------
    if (debug == 1) { break;}
    // --------------------------------------------------------------------------
  }
  Track_Pos->Write();
  std::cout << "<RESULT> The predicted pos of the muon is [x,y,z] = [" << xfin << "," << yfin << "," << zfin << "] and the predicted track is along [vx,vy,vz] = [" << vxmax << "," << vymax << "," << vzmax << "] with " << nmax << " number of hit pmts (probably). This is event number " << maxeve << ". Done!" << std::endl;
  if (debug == 1) {
    std::cout << "<DEBUG> DONE" << std::endl;
  }
}
