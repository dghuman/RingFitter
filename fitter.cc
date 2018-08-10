
// Just a macro that uses the PMT hit positions from the fiTQun modification until I can get it to work without fitqun to fit rings

#include <iostream>
#include <math.h>
#include <stdlib>
#include <vector>

void fitter(char *filename=NULL, char *outfile=NULL,bool debug=0) {
  TFile* file;
  TNtuple* Ntuple;
  TFile* NewFile;
  TNtuple* Track_Pos = new TNtuple("Track_Pos", "Predicted Muon Track and start Position", "xpos:ypos:zpos:xvec:yvec:zvec:nhits"); // An Ntuple to store the max hit vec, position and the number of corresponding hits.
  TH2D* ThetaPhi;   // The 2-Dim histogram that will keep track of the possible muon tracks that could produce the PMT hit observed at for a fixed position of interest. Will be remade for each position.

  Int_t counter = 0;
  Int_t btheta,bphi,z,maxeve;
  Float_t pmtxpos,pmtypos,pmtzpos,pmttime;
  Double_t length,maxbin,binval,t,deltat,deltaphi,deltatheta,ratio,theta,phi;
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
  Double_t oldvec[2];

  // Make the rotation matrix
  Double_t rotmat[3][3];

  int entries;
  bool fix;

  const Double_t posgran = 95.; // The spatial cube dimension that will be used as the starting position of each event track check
  const Double_t pi = 3.14159265535897;
  const Double_t thetaerror = 1.*pi/180.; // Error in the track angle reconstruction. This should be more dynamic based off of the position of interest and the pmt being hit. Const for now.
  const Double_t phierror = 1.*pi/180.;
  const Double_t chrnkvangle = 43.*pi/180..; // Angle to track at which we assume that the cherenkov light is being emmitted.
  const Double_t radius = sin(chrnkvangle); // The corresponding radius and height of the ring that would be produced on the unitsphere given a cherenkov path along the z-axis
  const Double_t height = cos(chrnkvangle);
  const Double_t xmax = 380.;
  const Double_t ymax = 380.; 
  const Double_t zmax = 525.;
  const Int_t num = 1000; // Number of points to take along the cone of muon tracks

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
  vector< vector<Double_t> > pmt;
  vector<Double_t> temppmt;


  for (int entry = 0; entry < entries; entry++) {
    Ntuple->GetEntry(entry);
    temppmt.clear();
    temppmt.push_back((Double_t)pmtxpos);
    temppmt.push_back((Double_t)pmtypos);
    temppmt.push_back((Double_t)pmtzpos);
    temppmt.push_back((Double_t)pmttime);
    pmt.push_back(temppmt);
  }
  // Switch to the new file
  NewFile->cd();
  TCanvas c1;

  std::cout << "ENTRIES == " << entries << std::endl;

  // Loop through the positions inside the detector and check the tracks that could have led to each pmt hit
  for (Double_t xpos = -xmax; xpos < xmax; xpos+=posgran) {
    std::cout << "<STATUS> " << ((xpos + xmax)*100)/((xmax*2)) << " % complete" << std::endl; // Just a status check
    for (Double_t ypos = -ymax; ypos < ymax; ypos+=posgran) {
      for (Double_t zpos = -zmax; zpos < zmax; zpos+=posgran) {
	// Make a histogram for this position
	ThetaPhi = new TH2D(Form("thetaphi_%d",counter), "Theta Phi Parameter Space", ((int) 2*pi/thetaerror) + 1, 0.0 - thetaerror/2., 2*pi + thetaerror/2., ((int) pi/phierror) + 1, 0.0 - phierror/2., pi - phierror/2.);
	// ---------------------------- DEBUG --------------------
	if (debug == 1) {
	  std::cout << "<DEBUG> Entering DEBUG mode. Will check pos [x,y,z] = [0,0,0] only." << std::endl;
	  xpos = 0.;
	  ypos = 0.;
	  zpos = 0.;
	}
	// -------------------------------------------------------
	for (int entry = 0; entry < entries; entry++) {
	  // Compute the Cherenkov vector
	  chrnkv[0] = pmt[entry][0] - xpos;
	  chrnkv[1] = pmt[entry][1] - ypos;
	  chrnkv[2] = pmt[entry][2] - zpos;
	  
	  // Normalize chrnkv vector
	  length = sqrt(chrnkv[0]*chrnkv[0] + chrnkv[1]*chrnkv[1] + chrnkv[2]*chrnkv[2]);
	  chrnkv[0] = chrnkv[0]/length;
	  chrnkv[1] = chrnkv[1]/length;
	  chrnkv[2] = chrnkv[2]/length;
	  // Find the Spherical coordinates of this unit vector
	  polarchrnkv[0] = atan2(chrnkv[1],chrnkv[0]);
	  polarchrnkv[1] = acos(chrnkv[2]);
	  // Check to see if the polar conversion happened correctly
	  if (fabs(sin(polarchrnkv[1])*cos(polarchrnkv[0]) - chrnkv[0]) > 0.0001 || fabs(sin(polarchrnkv[1])*sin(polarchrnkv[0]) - chrnkv[1]) > 0.0001 || fabs(cos(polarchrnkv[1]) - chrnkv[2]) > 0.0001) {
	    std::cout << "<WARNING> For Cherenkov [" << chrnkv[0] << "," << chrnkv[1] << "," << chrnkv[2] << "], the polar conversion failed." << std::endl;
	    std::cout << "<WARNING> Instead we got [" << sin(polarchrnkv[1])*cos(polarchrnkv[0]) << "," << sin(polarchrnkv[1])*sin(polarchrnkv[0]) << "," << cos(polarchrnkv[1]) << "]." << std::endl;
	  }
	  // ----------------------------------------------------------------------
	  //if (debug == 1) {
	  //std::cout << "<DEBUG> Cherenkov is [vx,vy,vz] = [" << chrnkv[0] << "," << chrnkv[1] << "," << chrnkv[2] << "]" << std::endl;
	  //std::cout << "<DEBUG> Polar cherenkov coordinates are [theta,phi] = [" << polarchrnkv[0] << "," << polarchrnkv[1] << "]" << std::endl;
	  //}
	  // ----------------------------------------------------------------------
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
	  /*rotmat[0][0] = cos(polarchrnkv[1])*cos(polarchrnkv[0]);
	  rotmat[0][1] = -sin(polarchrnkv[1]);
	  rotmat[0][2] = cos(polarchrnkv[1])*sin(polarchrnkv[0]);
	  rotmat[1][0] = sin(polarchrnkv[1])*cos(polarchrnkv[0]);
	  rotmat[1][1] = cos(polarchrnkv[1]);
	  rotmat[1][2] = sin(polarchrnkv[1])*sin(polarchrnkv[0]);
	  rotmat[2][0] = -sin(polarchrnkv[0]);
	  rotmat[2][1] = 0;
	  rotmat[2][2] = cos(polarchrnkv[0]);*/
	  
	  
	  // The idea here is to have a vector on our circle, say v, and find what it looks like in this rotated frame so that it is a predicted muon track. The circle is defined by using the predicted cherenkov path
	  // ------------------------------------- DEBUG -------------------------------
	  /*if (debug == 1) {
	    std::cout << "<DEBUG> Checking pmt POS [x,y,z] = [" << pmt[entry][0] << "," << pmt[entry][1] << "," << pmt[entry][2] << "]" << std::endl;
	    std::cout << "<DEBUG> Cherenkov is [vx,vy,vz] = [" << chrnkv[0] << "," << chrnkv[1] << "," << chrnkv[2] << "]" << std::endl;
	    std::cout << "<DEBUG> Polar cherenkov coordinates are [theta,phi] = [" << polarchrnkv[0] << "," << polarchrnkv[1] << "]" << std::endl;
	    }*/
	  /*if (debug == 1) {
	    vec[0] = 0;
	    vec[1] = 0;
	    vec[2] = 0;
	    tempvec[0] = 1.;
	    tempvec[1] = 0.;
	    tempvec[2] = 0.;
	    for (int i = 0; i < 3; i++) {
	      for (int j = 0; j < 3; j++) {
	      vec[i] += rotmat[i][j]*tempvec[j];
	      }
	      }
	      //if (debug == 1) {
	      //std::cout << "<DEBUG> !!Testing Rotation matrix. Input vector of x-axis [1,0,0]." << std::endl;
	      //std::cout << "<DEBUG> !!Should go to [" << chrnkv[0] << "," << chrnkv[1] << "," << chrnkv[2] << "]." << std::endl;
	      //std::cout << "<DEBUG> !!After rotation we get [" << vec[0] << "," << vec[1] << "," << vec[2] << "]." << std::endl;
	      //}
	      }*/
	  // --------------------------------------------------------------------------
	  t = 0.;
	  deltat = (2*pi)/((Double_t)num);
	  oldvec[0] = NULL;
	  oldvec[1] = NULL;
	  fix = 0;
	  do {
	    //for (double t = 0; t < 2*pi; t+=(2*pi/((double)num))) {
	    tempvec[0] = radius*cos(t);
	    tempvec[1] = radius*sin(t);
	    tempvec[2] = height;
	    // Now we rotate the vector by multiplying it through the rotation matrix Mv = v_new, which is just matrix vector multiplication
	    vec[0] = 0;
	    vec[1] = 0;
	    vec[2] = 0;
	    for (int i = 0; i < 3; i++) {
	      for (int j = 0; j < 3; j++) {
		vec[i] += rotmat[i][j]*tempvec[j];
	      }
	    }
	    theta = atan2(vec[1],vec[0]);
	    phi = acos(vec[2]);
	    if (oldvec[0] == NULL || oldvec[1] == NULL) {
	      oldvec[0] = theta;
	      oldvec[1] = phi;
	    }
	    // Check if the stepsize was appropriate and change it accordingly to fill gaps or make up for overcompensation
	    deltaphi = fabs(oldvec[1] - phi);
	    deltatheta = fabs(oldvec[0] - theta);
	    //std::cout << "t0 = " << t << " at entry " << entry << std::endl;
	    //std::cout << "deltat = " << deltat << " deltatheta = " << deltatheta << " deltaphi = " << deltaphi << " error = " << phierror << " theta = " << theta << " phi = " << phi << std::endl;
	    //std::cout << "oldvec[0] = " << oldvec[0] << " oldvec[1] = " << oldvec[1] << std::endl;
	    // Cant be dividing by 0, so ... 
	    if (!(deltaphi == 0 || deltatheta == 0) && fix == 0) {
	      // First check if the step size is too small
	      if (deltatheta < thetaerror/2. && deltaphi < phierror/2.) {
		fix = 1;
		if (deltatheta > deltaphi) {
		  ratio = thetaerror/(2.*deltatheta);
		  t += deltat*(ratio - 1.);		
		  //std::cout << "t1 = " << t << std::endl;
		  deltat = deltat*(ratio);
		  continue;
		} else {
		  ratio = phierror/(2.*deltaphi);
		  t += deltat*(ratio - 1.);
		  //std::cout << "t1 = " << t << std::endl;
		  deltat = deltat*(ratio);
		  continue;
		}
	      }
	      // Now if the stepsize is too large
	      if (deltatheta > 2.*thetaerror && deltaphi > 2.*phierror) {
		fix = 1;
		if (deltatheta > deltaphi) {
		  ratio = phierror/deltaphi;
		  t -= deltat*(1. - ratio);
		  //std::cout << "t1 = " << t << std::endl;
		  deltat = deltat*(ratio);
		  continue;
		} else {
		  ratio = thetaerror/deltatheta;
		  t -= deltat*(1. - ratio);
		  //std::cout << "t1 = " << t << std::endl;
		  deltat = deltat*(ratio);
		  continue;
		}
	      }
	    } else {fix = 0;}
	    //std::cout << "Made it past the checks! Pushing t forward!" << std::endl;
	    // Convert the vector to spherical coordinate and feed it into the histogram
	    if (vec[1] < 0) {
	      ThetaPhi->Fill(2*pi + atan2(vec[1],vec[0]),acos(vec[2]));
	    } else {
	      ThetaPhi->Fill(atan2(vec[1],vec[0]),acos(vec[2]));
	    }
		/* ---------------------------------------------------------------------------
		   if (debug == 1) {
		   std::cout << "<DEBUG> Possible track [vx,vy,vz] = [" << vec[0] << "," << vec[1] << "," << vec[2] << "]" << std::endl;
		   }
		   ------------------------------------------------------------------------------*/
	    if (fabs(sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]) - 1) > 0.01) {
	      std::cout << "<WARNING> Track is of length != 1. This should never happen. Current position is [x,y,z] = [" << xpos << "," << ypos << "," << zpos << "]." << std::endl;
	    }
	    t += deltat;
	    oldvec[0] = theta;
	    oldvec[1] = phi;
	  } while ( t < 2*pi);
	  if (debug == 1) {
	    ThetaPhi->Draw("colz");
	    c1->Print(Form("./plots/testingring1/ring_%d.png",entry));
	  }
	}
	counter++;
	maxbin = ThetaPhi->GetMaximumBin(); // This should be the bin with the most number of PMTs hit
	binval = ThetaPhi->GetBinContent(maxbin);
	ThetaPhi->GetBinXYZ(maxbin,btheta,bphi,z);
	//std::cout << "Estimated theta: " << btheta*thetaerror << " phi: " << bphi*phierror << std::endl;
	Track_Pos->Fill((Float_t)xpos,(Float_t)ypos,(Float_t)zpos,(Float_t)sin(bphi*phierror)*(Float_t)cos(btheta*thetaerror),(Float_t)sin(bphi*phierror)*(Float_t)sin(btheta*thetaerror),(Float_t)cos(bphi*phierror),(Float_t)binval);
	// Keep track of the position and vector with the most matched PMTs
	if (binval > nmax) {
	  xfin = xpos;
	  yfin = ypos;
	  zfin = zpos;
	  vxmax = sin(bphi*phierror)*cos(btheta*thetaerror);
	  vymax = sin(bphi*phierror)*sin(btheta*thetaerror);
	  vzmax = cos(bphi*phierror);
	  nmax = binval;
	  maxeve = counter;
	}
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
