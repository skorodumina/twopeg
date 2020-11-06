#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "interpol.h"


using namespace std;


//-------
Short_t getWbin (Float_t W) {
//return int(W*10000. - 1.4125*10000.)/250;
//return int((W-1.4125)/0.025);
Short_t bin;

if ((W>=1.4125)&&(W<=1.4375)) bin = 0;
if ((W>=1.4375)&&(W<=1.4625)) bin = 1;
if ((W>=1.4625)&&(W<=1.4875)) bin = 2;
if ((W>=1.4875)&&(W<=1.5125)) bin = 3;

if ((W>=1.5125)&&(W<=1.5375)) bin = 4;
if ((W>=1.5375)&&(W<=1.5625)) bin = 5;
if ((W>=1.5625)&&(W<=1.5875)) bin = 6;
if ((W>=1.5875)&&(W<=1.6125)) bin = 7;

if ((W>=1.6125)&&(W<=1.6375)) bin = 8;
if ((W>=1.6375)&&(W<=1.6625)) bin = 9;
if ((W>=1.6625)&&(W<=1.6875)) bin = 10;
if ((W>=1.6875)&&(W<=1.7125)) bin = 11;

if ((W>=1.7125)&&(W<=1.7375)) bin = 12;
if ((W>=1.7375)&&(W<=1.7625)) bin = 13;
if ((W>=1.7625)&&(W<=1.7875)) bin = 14;
if ((W>=1.7875)&&(W<=1.8125)) bin = 15;

if ((W<1.4125)||(W>1.8125)) {
cout << "Error, wrong W range, Ripani" << "\n";
bin = -100;
};
return bin;
};

//-------
Short_t getQ2bin (Float_t Q2) {

Short_t bin;

if ((Q2 >= 0.65) && (Q2 < 0.95)) bin = 0;
if ((Q2 >= 0.95) && (Q2 <= 1.3)) bin = 1;
if ((Q2 < 0.65) || (Q2 > 1.3)) {
cout << "Error, wrong Q2 range, Ripani" << "\n";
bin = -100;
};
return bin;
};

//---------
Short_t getsbin (Short_t Wbin, Float_t sgen, Float_t Smax, Float_t Smin) {

Short_t bin;
if ((sgen>=Smin)&&(sgen<=Smax)) bin = int((sgen-Smin)/((Smax - Smin)/11.));
if (sgen<Smin) bin = 0;
if (sgen>Smax) bin = 10;

return bin;
};
//--------------------
Short_t getanglebin (Float_t anglegen, Float_t anglemax) {

Short_t bin;
if ((anglegen < 0.01)) bin = 0;
if ((anglegen > anglemax - 0.01)) bin = 4;
if ((anglegen >= 0.01) && (anglegen <= anglemax - 0.01)) bin = int((anglegen - 0.01)/((anglemax  - 0.02)/5.));

return bin;
};

//This subroutine gets Ripani cross sections on the tabuted grid and interpolates them to a desired point within (W>=1.4125)&&(W<=1.8125)&&(Q2>=0.65)&&(Q2<=1.3)
//This is the grid and xsect arrays for Ripani cross sections:
// W_ARR[17];
// Q2_ARR[3];
// S12_ARR[12][17];
// S23_ARR[12][17];
// THETA_ARR[6]; 
// ALPHA_ARR[6];
// SIGMA_ARR[6][3][17][12][12][6][6];

//The subrouting is doing the following:
//1 - using auxiliary functions (getWbin, getQ2bin, getsbin, getanglebin) we identify the number of left and right point 
//2 - then we are doing 4d-interpolation for each (Wleft_bin, Q2_left_bin),  (Wright_bin, Q2_left_bin), (Wright_bin, Q2_right_bin) and (Wleft_bin, Q2_right_bin) and obtain cross-secton in that points (4 GLOBAL 6dim arrays)
//3 - then we are doing 2d-interpolation for the points written above and obtain sigma_final[6]


void get_xsect_ripani(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final){


//using auxiliary functions (getWbin, getQ2bin, getsbin, getanglebin) we identify the number of left and right point 
Short_t Wleft_bin = getWbin(Wgen);
Short_t Wright_bin = Wleft_bin+1;

//cout << Wgen <<" "<< Wleft_bin << " " <<Wright_bin <<"\n";

Short_t Q2left_bin = getQ2bin(Q2gen);
Short_t Q2right_bin = Q2left_bin+1;
//cout << Q2gen <<" "<< Q2left_bin << " " <<Q2right_bin <<"\n";

Short_t s12left_wleft_bin = getsbin(Wleft_bin, s12gen, S12_ARR[11][Wleft_bin], S12_ARR[0][Wleft_bin]);
Short_t s12right_wleft_bin = s12left_wleft_bin +1;

Short_t s12left_wright_bin = getsbin(Wright_bin, s12gen, S12_ARR[11][Wright_bin], S12_ARR[0][Wright_bin]);
Short_t s12right_wright_bin = s12left_wright_bin +1;

Short_t s23left_wleft_bin = getsbin(Wleft_bin, s23gen, S23_ARR[11][Wleft_bin], S23_ARR[0][Wleft_bin]);
Short_t s23right_wleft_bin = s23left_wleft_bin +1;

Short_t s23left_wright_bin = getsbin(Wright_bin, s23gen, S23_ARR[11][Wright_bin], S23_ARR[0][Wright_bin]);
Short_t s23right_wright_bin = s23left_wright_bin +1;

//cout << s23left_wleft_bin <<" "<< Wgen << " "<<Wleft_bin<< " " << s23gen << " "<< S23_ARR[s23left_wleft_bin][Wleft_bin] << " "<< S23_ARR[s23right_wleft_bin][Wleft_bin] << "\n";


Short_t thetaleft_bin = getanglebin(thetagen,THETA_ARR[5]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin(alphagen,ALPHA_ARR[5]);
Short_t alpharight_bin = alphaleft_bin+1;

//cout << alphaleft_bin << " "<< alphagen<< " "<< ALPHA_ARR[alphaleft_bin] << " "<< ALPHA_ARR[ alpharight_bin] << "\n";

Float_t sigma_final[6];

//then we are doing 4d-interpolation for each (Wleft_bin, Q2_left_bin),  (Wright_bin, Q2_left_bin), (Wright_bin, Q2_right_bin) and (Wleft_bin, Q2_right_bin) and obtain cross-secton in that points (4 GLOBAL 6dim arrays)
//0 - sigma_t, 1 - sigma_l, 2 - sigma_c2f, 3 - sigma_s2f, 4 - sigma_cf, 5 - sigma_sf
for (Short_t i=0;i<6;i++){
interpol(4,Q2left_bin,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_q2left[i],i);

interpol(4,Q2right_bin,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_q2right[i],i);

interpol(4,Q2right_bin,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_q2right[i],i);

interpol(4,Q2left_bin,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_q2left[i],i);
};


//cout <<sigma_wleft_q2left[0] <<"  qqgett\n";
//then we are doing 2d-interpolation for the points written above and obtain sigma_final[6]
//0 - sigma_t, 1 - sigma_l, 2 - sigma_c2f, 3 - sigma_s2f, 4 - sigma_cf, 5 - sigma_sf


for (Short_t i=0;i<6;i++){
interpol(2,0,0, Wleft_bin, Wright_bin, Q2left_bin, Q2right_bin, 0,  0, 0,  0, Wgen, Q2gen, 0., 0.,sigma_final[i], i);
};
//cout << sigma_wright_q2left[0] << "  "<< sigma_wright_q2right[0] << "  "<< sigma_wleft_q2right[0]<< " "<< sigma_wleft_q2left[0] << " \n";
//cout << Wgen<<" "<< sigma_final[0] << " qqq \n";

//We get explicitly different sigmas from the array
 sigma_t_final = sigma_final[0];
 sigma_l_final = sigma_final[1];
 sigma_c2f_final = sigma_final[2];
 sigma_s2f_final = sigma_final[3];
 sigma_cf_final = sigma_final[4];
 sigma_sf_final = sigma_final[5];


};


