#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "interpol_golovach.h"
#include "interpol_gol2.h"
#include "get_xsect_ripani.h"
using namespace std;

//-------
Short_t getWbin_GOL (Float_t W) {
Short_t bin;
if ((W>=1.6125)&&(W<=2.1375))  bin = int((W-1.6125)/0.025);
if ((W>2.1375)&&(W<=2.5375))  bin = 21+int((W-2.1375)/0.05);

/*
if ((W>=1.6125)&&(W<=1.6375)) bin = 0;
if ((W>=1.6375)&&(W<=1.6625)) bin = 1;
if ((W>=1.6625)&&(W<=1.6875)) bin = 2;
if ((W>=1.6875)&&(W<=1.7125)) bin = 3;

if ((W>=1.7125)&&(W<=1.7375)) bin = 4;
if ((W>=1.7375)&&(W<=1.7625)) bin = 5;
if ((W>=1.7625)&&(W<=1.7875)) bin = 6;
if ((W>=1.7875)&&(W<=1.8125)) bin = 7;

if ((W>=1.8125)&&(W<=1.8375)) bin = 8;
if ((W>=1.8375)&&(W<=1.8625)) bin = 9;
if ((W>=1.8625)&&(W<=1.8875)) bin = 10;
if ((W>=1.8875)&&(W<=1.9125)) bin = 11;

if ((W>=1.9125)&&(W<=1.9375)) bin = 12;
if ((W>=1.9375)&&(W<=1.9625)) bin = 13;
if ((W>=1.9625)&&(W<=1.9875)) bin = 14;
if ((W>=1.9875)&&(W<=2.0125)) bin = 15;

if ((W>=2.0125)&&(W<=2.0375)) bin = 16;
if ((W>=2.0375)&&(W<=2.0625)) bin = 17;
if ((W>=2.0625)&&(W<=2.0875)) bin = 18;
if ((W>=2.0875)&&(W<=2.1125)) bin = 19;

if ((W>=2.1125)&&(W<=2.1375)) bin = 20;

if ((W>=2.1375)&&(W<=2.1875)) bin = 21;
if ((W>=2.1875)&&(W<=2.2375)) bin = 22;
if ((W>=2.2375)&&(W<=2.2875)) bin = 23;
if ((W>=2.2875)&&(W<=2.3375)) bin = 24;
if ((W>=2.3375)&&(W<=2.3875)) bin = 25;
if ((W>=2.3875)&&(W<=2.4375)) bin = 26;
if ((W>=2.4375)&&(W<=2.4875)) bin = 27;
if ((W>=2.4875)&&(W<=2.5375)) bin = 28;
*/
if ((W<1.6125)||(W>2.5375)) {
cout << "Error, wrong W range, Golovach" << "\n";
bin = -100;
};
return bin;
};

//-------
Short_t getWbin_GOL2 (Float_t W) {

Short_t bin;
if ((W>2.5875)&&(W<=3.0375))  bin = int((W-2.5875)/0.05);

if ((W<2.5875)||(W>3.0375)) {
cout << "Error, wrong W range, WGOL2" << "\n";
bin = -100;
}
return bin;
};

//---------
Short_t getsbin_GOL (Float_t sgen, Float_t Smax, Float_t Smin) {

Short_t bin;
if ((sgen>=Smin)&&(sgen<=Smax)) bin = int((sgen-Smin)/((Smax - Smin)/15.));
if (sgen<Smin) bin = 0;
if (sgen>Smax) bin = 14;

return bin;
};
//--------------------
Short_t getanglebin_GOL (Float_t anglegen, Float_t anglemax) {

Short_t bin;
if ((anglegen < 0.01)) bin = 0;
if ((anglegen > anglemax - 0.01)) bin = 12;
if ((anglegen >= 0.01) && (anglegen <= anglemax - 0.01)) bin = int((anglegen - 0.01)/((anglemax  - 0.02)/13.));

return bin;
};

//This file contains two similar subroutine for geting photoprofuction cross sections.
//get_xsect_gol_datamod(...) - gets Golovach cross sections for W from 1.6125 to 2.5375 GeV
//get_xsect_gol_model(...) - takes cross sections in the region W from 2.5875 to 3.0375 GeV and Q2 = 0 GeV^2 (pure-model based)

//This subroutine gets Golovach cross sections on the tabuted grid and interpolates them to a desired point within W from 1.6125 to 2.5375 GeV
//This is the grid and xsect array for Golovach cross sections:
//W_ARR_GOL[30];
//S12_ARR_GOL[16][30];
//S23_ARR_GOL[16][30];
//THETA_ARR_GOL[14]; 
//ALPHA_ARR_GOL[14];
//SIGMA_ARR_GOL[30][16][16][14][14];
 

void get_xsect_gol_datamod(Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Short_t &Wleft_bin, Float_t &sigma_wright_gol,Float_t &sigma_wleft_gol){ 

//using auxiliary functions (getWbin_GOL, getsbin_GOL, getanglebin_GOL) we identify the left and right point for each variable generated  (according to Golovach binning)
Wleft_bin = getWbin_GOL(Wgen);
Short_t Wright_bin = Wleft_bin+1;

//cout << Wgen << " "<< Wleft_bin<< " "<< W_ARR_GOL[Wleft_bin] <<" "<< W_ARR_GOL[Wright_bin]<< "\n";

Short_t s12left_wleft_bin = getsbin_GOL(s12gen, S12_ARR_GOL[15][Wleft_bin], S12_ARR_GOL[0][Wleft_bin]);
Short_t s12right_wleft_bin = s12left_wleft_bin +1;

//cout <<s12left_wleft_bin <<" "<<s12gen<<" "<<S12_ARR_GOL[s12left_wleft_bin][Wleft_bin]<< " "<< S12_ARR_GOL[s12right_wleft_bin][Wleft_bin]<<"\n";
Short_t s12left_wright_bin = getsbin_GOL(s12gen, S12_ARR_GOL[15][Wright_bin], S12_ARR_GOL[0][Wright_bin]);
Short_t s12right_wright_bin = s12left_wright_bin +1;

Short_t s23left_wleft_bin = getsbin_GOL(s23gen, S23_ARR_GOL[15][Wleft_bin], S23_ARR_GOL[0][Wleft_bin]);
Short_t s23right_wleft_bin = s23left_wleft_bin +1;

Short_t s23left_wright_bin = getsbin_GOL(s23gen, S23_ARR_GOL[15][Wright_bin], S23_ARR_GOL[0][Wright_bin]);
Short_t s23right_wright_bin = s23left_wright_bin +1;

//cout << s23left_wleft_bin <<" "<< Wgen << " "<<Wleft_bin<< " " << s23gen << " "<< S23_ARR_GOL[s23left_wleft_bin][Wleft_bin] << " "<< S23_ARR_GOL[s23right_wleft_bin][Wleft_bin] << "\n";

Short_t thetaleft_bin = getanglebin_GOL(thetagen,THETA_ARR_GOL[13]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin_GOL(alphagen,ALPHA_ARR_GOL[13]);
Short_t alpharight_bin = alphaleft_bin+1;

//cout << thetagen <<" "<< thetaleft_bin<<" "<< thetaright_bin<<" \n";

//then we are doing 4d-interpolation for each Wleft_bin and Wright_bin and obtain sigma_wright_gol and sigma_wleft_gol
interpol_golovach(Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_gol);

interpol_golovach(Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_gol);
};


//This subroutine takes cross sections in the region W from 2.5875 to 3.0375 GeV and Q2 = 0 GeV^2 (pure-model based)
//This is the grid and xsect array for this cross sections:
//W_ARR_RIP3[10];
//S12_ARR_RIP3[16][10];
//S23_ARR_RIP3[16][10];
//THETA_ARR[6]; 
//ALPHA_ARR[6];
//SIGMA_ARR_GOL2[10][16][16][6][6]; (the same as SIGMA_ARR_RIP3[6][10][16][16][6][6], but with scaling factors, which are W, s12/s23, and angles dependent)

void get_xsect_gol_model(Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Short_t &Wleft_bin, Float_t &sigma_wright_gol,Float_t &sigma_wleft_gol){ 

Wleft_bin = getWbin_GOL2(Wgen);
Short_t Wright_bin = Wleft_bin+1;
//cout << Wgen << " "<< Wleft_bin<< "\n";

Short_t s12left_wleft_bin = getsbin_GOL(s12gen, S12_ARR_RIP3[15][Wleft_bin], S12_ARR_RIP3[0][Wleft_bin]);
Short_t  s12right_wleft_bin = s12left_wleft_bin +1;

Short_t  s12left_wright_bin = getsbin_GOL(s12gen, S12_ARR_RIP3[15][Wright_bin], S12_ARR_RIP3[0][Wright_bin]);
Short_t  s12right_wright_bin = s12left_wright_bin +1;

Short_t  s23left_wleft_bin = getsbin_GOL(s23gen, S23_ARR_RIP3[15][Wleft_bin], S23_ARR_RIP3[0][Wleft_bin]);
Short_t  s23right_wleft_bin = s23left_wleft_bin +1;

Short_t  s23left_wright_bin = getsbin_GOL(s23gen, S23_ARR_RIP3[15][Wright_bin], S23_ARR_RIP3[0][Wright_bin]);
Short_t  s23right_wright_bin = s23left_wright_bin +1;

//cout << s23left_wleft_bin <<" "<< Wgen << " "<<Wleft_bin<< " " << s23gen << " "<< S23_ARR_RIP3[s23left_wleft_bin][Wleft_bin] << " "<< S23_ARR_RIP3[s23right_wleft_bin][Wleft_bin] << "\n";

Short_t thetaleft_bin = getanglebin(thetagen,THETA_ARR[5]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin(alphagen,ALPHA_ARR[5]);
Short_t alpharight_bin = alphaleft_bin+1;

//then we are doing 4d-interpolation for each Wleft_bin and Wright_bin and obtain sigma_wright_gol and sigma_wleft_gol
interpol_gol2(Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_gol);

interpol_gol2(Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_gol);
};

