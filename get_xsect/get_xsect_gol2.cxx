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

Short_t bin;
if ((W>2.5875)&&(W<=3.0375))  bin = int((W-2.5875)/0.05);

if ((W<2.5875)||(W>3.0375)) {
cout << "Error, wrong W range, WGOL2" << "\n";
bin = -100;
}

return bin;
};

//Note that SIGMA_ARR_GOL2 is the same as SIGMA_ARR_RIP3, but with scaling coefficients (W, s12, 23, and angles dependent). Both SIGMA sets are taken from the model without relying on the data.


void get_xsect_gol2(Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Short_t &Wleft_bin, Float_t &sigma_wright_gol,Float_t &sigma_wleft_gol){ 


//using auxiliary functions (getWbin_GOL, getsbin_GOL, getanglebin_GOL) we identify the left and right point for each variable generated  (according to Golovach binning)
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


Short_t thetaleft_bin = getanglebin(thetagen,THETA_ARR[5]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin(alphagen,ALPHA_ARR[5]);
Short_t alpharight_bin = alphaleft_bin+1;

//cout << s23left_wleft_bin <<" "<< Wgen << " "<<Wleft_bin<< " " << s23gen << " "<< S23_ARR_RIP3[s23left_wleft_bin][Wleft_bin] << " "<< S23_ARR_RIP3[s23right_wleft_bin][Wleft_bin] << "\n";

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
