#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"

#include "interpol_gol2.h"
#include "get_xsect_gol2.h"
#include "get_xsect_golovach.h"

#include "get_xsect_ripani.h"
using namespace std;




//-------
Short_t getWbin_GOL2 (Float_t W) {
//return int(W*10000. - 1.4125*10000.)/250;
//if ((W>=1.6125)&&(W<=2.1375))  return int((W-1.6125)/0.025);
if ((W>2.5875)&&(W<=3.0375))  return int((W-2.5875)/0.05);

/*

if ((W>=1.6125)&&(W<=1.6375)) return 0;
if ((W>=1.6375)&&(W<=1.6625)) return 1;
if ((W>=1.6625)&&(W<=1.6875)) return 2;
if ((W>=1.6875)&&(W<=1.7125)) return 3;

if ((W>=1.7125)&&(W<=1.7375)) return 4;
if ((W>=1.7375)&&(W<=1.7625)) return 5;
if ((W>=1.7625)&&(W<=1.7875)) return 6;
if ((W>=1.7875)&&(W<=1.8125)) return 7;

if ((W>=1.8125)&&(W<=1.8375)) return 8;
if ((W>=1.8375)&&(W<=1.8625)) return 9;
if ((W>=1.8625)&&(W<=1.8875)) return 10;
if ((W>=1.8875)&&(W<=1.9125)) return 11;

if ((W>=1.9125)&&(W<=1.9375)) return 12;
if ((W>=1.9375)&&(W<=1.9625)) return 13;
if ((W>=1.9625)&&(W<=1.9875)) return 14;
if ((W>=1.9875)&&(W<=2.0125)) return 15;

if ((W>=2.0125)&&(W<=2.0375)) return 16;
if ((W>=2.0375)&&(W<=2.0625)) return 17;
if ((W>=2.0625)&&(W<=2.0875)) return 18;
if ((W>=2.0875)&&(W<=2.1125)) return 19;

if ((W>=2.1125)&&(W<=2.1375)) return 20;

if ((W>=2.1375)&&(W<=2.1875)) return 21;
if ((W>=2.1875)&&(W<=2.2375)) return 22;
if ((W>=2.2375)&&(W<=2.2875)) return 23;
if ((W>=2.2875)&&(W<=2.3375)) return 24;
if ((W>=2.3375)&&(W<=2.3875)) return 25;
if ((W>=2.3875)&&(W<=2.4375)) return 26;
if ((W>=2.4375)&&(W<=2.4875)) return 27;
if ((W>=2.4875)&&(W<=2.5375)) return 28;

*/
if ((W<2.5875)||(W>3.0375)) {
cout << "Error, wrong W range" << "\n";
return -100;
}

};

/*
//---------
Short_t getsbin_GOL (Float_t sgen, Float_t Smax, Float_t Smin) {
if ((sgen>=Smin)&&(sgen<=Smax)) return int((sgen-Smin)/((Smax - Smin)/15.));
if (sgen<Smin) return 0;
if (sgen>Smax) return 14;


};
//--------------------
Short_t getanglebin_GOL (Float_t anglegen, Float_t anglemax) {

if ((anglegen < 0.01)) return 0;
if ((anglegen > anglemax - 0.01)) return 12;
if ((anglegen >= 0.01) && (anglegen <= anglemax - 0.01)) return int((anglegen - 0.01)/((anglemax  - 0.02)/13.));


};
*/

//This subrouting is doing the following:
//1 - using auxiliary functions (getWbin_GOL, getsbin_GOL, getanglebin_GOL) we identify the left and right point for each variable generated  (according to Golovach binning)
//2 - then we are doing 4d-interpolation for each Wleft_bin and Wright_bin and obtain sigma_wright_gol and sigma_wleft_gol


void get_xsect_gol2(Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Short_t &Wleft_bin, Float_t &sigma_wright_gol,Float_t &sigma_wleft_gol){ 


//using auxiliary functions (getWbin_GOL, getsbin_GOL, getanglebin_GOL) we identify the left and right point for each variable generated  (according to Golovach binning)
Wleft_bin = getWbin_GOL2(Wgen);
Short_t Wright_bin = Wleft_bin+1;

//cout << Wgen << " "<< Wleft_bin<< " "<< W_ARR_GOL[Wleft_bin] <<" "<< W_ARR_GOL[Wright_bin]<< "\n";
Short_t s12left_wleft_bin = getsbin_GOL(s12gen, S12_ARR_RIP3[15][Wleft_bin], S12_ARR_RIP3[0][Wleft_bin]);
Short_t  s12right_wleft_bin = s12left_wleft_bin +1;

Short_t  s12left_wright_bin = getsbin_GOL(s12gen, S12_ARR_RIP3[15][Wright_bin], S12_ARR_RIP3[0][Wright_bin]);
Short_t  s12right_wright_bin = s12left_wright_bin +1;
 


Short_t  s23left_wleft_bin = getsbin_GOL(s23gen, S23_ARR_RIP3[15][Wleft_bin], S23_ARR_RIP3[0][Wleft_bin]);
Short_t  s23right_wleft_bin = s23left_wleft_bin +1;

Short_t  s23left_wright_bin = getsbin_GOL(s23gen, S23_ARR_RIP3[15][Wright_bin], S23_ARR_RIP3[0][Wright_bin]);
Short_t  s23right_wright_bin = s23left_wright_bin +1;


Short_t thetaleft_bin = getanglebin(thetagen,THETA_ARR[5]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin(alphagen,ALPHA_ARR[5]);
Short_t alpharight_bin = alphaleft_bin+1;

//cout << s23left_wleft_bin <<" "<< Wgen << " "<<Wleft_bin<< " " << s23gen << " "<< S23_ARR_GOL[s23left_wleft_bin][Wleft_bin] << " "<< S23_ARR_GOL[s23right_wleft_bin][Wleft_bin] << "\n";

//then we are doing 4d-interpolation for each Wleft_bin and Wright_bin and obtain sigma_wright_gol and sigma_wleft_gol
interpol_gol2(Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_gol);

interpol_gol2(Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_gol);

//cout <<Wgen<<" "<< sigma_wright_gol << " "<<sigma_wleft_gol<<"\n";

//cout << Wgen << " "<< s12left_wright_bin << " "<<s12right_wright_bin<<" "<<S12_ARR_GOL[s12right_wright_bin][Wright_bin] << "  "<<s23left_wright_bin<< " "<< s23right_wright_bin<< " "<< sigma_wright_gol<<"\n";
//cout << Wgen << " "<< s12left_wleft_bin << " "<<s12right_wleft_bin<<" "<<S12_ARR_GOL[s12right_wleft_bin][Wleft_bin]<< "  "<<s23left_wleft_bin<< " "<< s23right_wleft_bin<< " "<< sigma_wleft_gol<<"\n";

//for (Short_t i=0;i<=13;i++){

//cout<< i<<" "<<S12_ARR[i][2]<<"\n";
//};


//sigma_t_gol = 1./fabs(W_ARR_GOL[Wright_bin]-W_ARR_GOL[Wleft_bin]);
//sigma_t_gol = sigma_t_gol*(sigma_wright_gol*fabs(W_ARR_GOL[Wleft_bin]-Wgen)+sigma_wleft_gol*fabs(W_ARR_GOL[Wright_bin]-Wgen));
//cout << sigma_t_gol << " "<< sigma_wright_gol << " "<< sigma_wleft_gol<<"\n";
return;
};
