#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"

#include "get_xsect_ripani.h"
#include "get_xsect_golovach.h"
#include "get_xsect_scale_gol_18_25.h"
#include "get_xsect_14_18_lowq2_fit.h"
using namespace std;
void get_xsect_scale_gol_18_25(Float_t E_beam, Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_total){

Float_t  sigma_t_wright_gol,sigma_t_wleft_gol;
Float_t Ratio_l, Ratio_c2f, Ratio_cf;
Short_t w_left_bin;
Short_t w_right_bin;


 get_xsect_golovach(Wgen, s12gen,s23gen, thetagen, alphagen, w_left_bin, sigma_t_wright_gol,sigma_t_wleft_gol );
 w_right_bin = w_left_bin+1;
 
// cout << Wgen << " "<< w_left_bin<< "\n";
// cout << func_sigma_t (Q2gen,16) <<"\n";
 
sigma_total = 1./fabs(W_ARR_GOL[w_right_bin]-W_ARR_GOL[w_left_bin]);
sigma_total = sigma_total*(sigma_t_wright_gol*fabs(W_ARR_GOL[w_left_bin]-Wgen)+sigma_t_wleft_gol*fabs(W_ARR_GOL[w_right_bin]-Wgen));


Ratio_l = getEpsL(E_beam,W_ARR[16],Q2gen)*pol2(Q2gen,16,1)/getEpsL(2.445,W_ARR[16],Q2gen)/func_sigma_t(Q2gen,16);

Ratio_c2f = getEpsT(E_beam,W_ARR[16],Q2gen)*(pol2(Q2gen,16,2)-1.)/getEpsT(2.445,W_ARR[16],Q2gen);

Ratio_cf = sqrt(2*getEpsL(E_beam,W_ARR[16],Q2gen)*(getEpsT(E_beam,W_ARR[16],Q2gen)+1.))*(pol2(Q2gen,16,4)-1.)/ sqrt(2*getEpsL(2.445,W_ARR[16],Q2gen)*(getEpsT(2.445,W_ARR[16],Q2gen)+1.));




sigma_total = sigma_total*func_sigma_t(Q2gen,16)/func_sigma_t(0,16)*(1. + Ratio_l + Ratio_c2f*cos(2.*phigen) + Ratio_cf*cos(phigen));
//sigma_total = sigma_total*func_sigma_t(Q2gen,16)/func_sigma_t(0,16);

 //cout <<Wgen<<" "<< sigma_t_wright_gol << " "<<sigma_t_wleft_gol<<"\n";
}; 
