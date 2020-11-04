#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "interpol.h"
#include "get_xsect_ripani.h"
#include "get_xsect_golovach.h"
#include "interpol_int.h"

using namespace std;

Float_t getEpsL (Float_t E_beam, Float_t W, Float_t Q2) {

Float_t nu = (W*W + Q2 - MP*MP)/2./MP;
Float_t theta_el = acos(1.- Q2/E_beam/(E_beam - nu)/2.);
Float_t eps_t = 1./(1.+ 2.*(1. + nu*nu/Q2)*tan(theta_el/2.)*tan(theta_el/2.));
Float_t eps_l = Q2*eps_t/nu/nu;

//cout << E_beam<<" "<<W<<" "<<Q2<<" "<<nu<<" "<<theta_el<<" " <<Q2/E_beam/(E_beam - nu)/2 <<" "<<eps_t<<" "<< eps_l << "\n";
return eps_l;
};

Float_t getEpsT (Float_t E_beam, Float_t W, Float_t Q2) {

Float_t nu = (W*W + Q2 - MP*MP)/2./MP;
Float_t theta_el = acos(1.- Q2/E_beam/(E_beam - nu)/2.);
Float_t eps_t = 1./(1.+ 2.*(1. + nu*nu/Q2)*tan(theta_el/2.)*tan(theta_el/2.));
Float_t eps_l = Q2*eps_t/nu/nu;
return eps_t;
};


//-------------------------------------------------------------------------
Float_t func_sigma_t (Float_t  x, Short_t wbin) {
Float_t func = pow((x+FIT_PARAM_SIGMA_T[0][wbin]), FIT_PARAM_SIGMA_T[1][wbin])*FIT_PARAM_SIGMA_T[2][wbin]+FIT_PARAM_SIGMA_T[3][wbin];
return func;
};

Float_t pol2 (Float_t  x, Short_t wbin, Short_t i) {
 Float_t func;
if (i==1) func =  FIT_PARAM_SIGMA_L[2][wbin]*x*x + FIT_PARAM_SIGMA_L[1][wbin]*x + FIT_PARAM_SIGMA_L[0][wbin]; 
if (i==2) func =  FIT_PARAM_SIGMA_C2F[2][wbin]*x*x + FIT_PARAM_SIGMA_C2F[1][wbin]*x + FIT_PARAM_SIGMA_C2F[0][wbin];  
if (i==4) func =  FIT_PARAM_SIGMA_CF[2][wbin]*x*x + FIT_PARAM_SIGMA_CF[1][wbin]*x + FIT_PARAM_SIGMA_CF[0][wbin];  
return func;
};

//This subroutine estimates cross section in the region (W>=1.6125)&&(W<=1.8125)&&(Q2>around 0)&&(Q2<0.65)
//First, it takes Ripani xsect for Q2 = 0.65 GeV^{2} (Q2_bin = 0) for Wleft_bin and Wright_bin
//Then it scales all sigma[0 to 6] to Q2 = zero with particular scaling functions and performs their W-interpolation
//After that it takes Golovach photoproduction cross sections for Wright_bin and Wleft_bin and scales them with Q2 (only sigma_t). Their W interpolation is also performed.
//Then mixing of scaled Ripani sigma_t and scaled Golovach sigma_t takes place

void get_xsect_w16_18_lowq2_fit(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen,Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final){


//using auxiliary functions (getWbin, getQ2bin, getsbin, getanglebin) we identify the number of left and right point 
Short_t Wleft_bin = getWbin(Wgen);
Short_t Wright_bin = Wleft_bin+1;

Short_t w_left_bin_gol;
Float_t A_tmp[9];
Float_t eps_t_tmp,eps_l_tmp,nu_tmp,theta_el_tmp;
Float_t diff_xsect_gol, sigma_t_wright_gol,sigma_t_wleft_gol;
//cout << Wgen <<" "<< Wleft_bin << " " <<Wright_bin <<"\n";

Short_t Q2_bin = 0;


Short_t s12left_wleft_bin = getsbin(Wleft_bin, s12gen, S12_ARR[11][Wleft_bin], S12_ARR[0][Wleft_bin]);
Short_t s12right_wleft_bin = s12left_wleft_bin +1;

Short_t s12left_wright_bin = getsbin(Wright_bin, s12gen, S12_ARR[11][Wright_bin], S12_ARR[0][Wright_bin]);
Short_t s12right_wright_bin = s12left_wright_bin +1;

Short_t s23left_wleft_bin = getsbin(Wleft_bin, s23gen, S23_ARR[11][Wleft_bin], S23_ARR[0][Wleft_bin]);
Short_t s23right_wleft_bin = s23left_wleft_bin +1;

Short_t s23left_wright_bin = getsbin(Wright_bin, s23gen, S23_ARR[11][Wright_bin], S23_ARR[0][Wright_bin]);
Short_t s23right_wright_bin = s23left_wright_bin +1;

Short_t thetaleft_bin = getanglebin(thetagen,THETA_ARR[5]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin(alphagen,ALPHA_ARR[5]);
Short_t alpharight_bin = alphaleft_bin+1;


Float_t sigma_final[6];
Float_t sigma_wright[6],sigma_wleft[6];
//then we are doing 4d-interpolation for Wleft_bin and Wright_bin and obtain cross-secton in that points (sigma_wright[i] and sigma_wleft[i])
//0 - sigma_t, 1 - sigma_l, 2 - sigma_c2f, 3 - sigma_s2f, 4 - sigma_cf, 5 - sigma_sf
for (Short_t i=0;i<6;i++){
interpol(4,Q2_bin,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright[i],i);

interpol(4,Q2_bin,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft[i],i);
};

//Here we Q2-scale sigma_t, sigma_t, sigma_c2f and sigma_cf with the fit_functions (func_sigma_t and corresponding pol2) for Wright_bin and Wleft_bin

//sigma_t
sigma_wright[0] = sigma_wright[0]*func_sigma_t(Q2gen,Wright_bin)/func_sigma_t(0.65,Wright_bin);
sigma_wleft[0] = sigma_wleft[0]*func_sigma_t(Q2gen,Wleft_bin)/func_sigma_t(0.65,Wleft_bin);


//sigma_l
sigma_wright[1] = sigma_wright[1]*pol2(Q2gen,Wright_bin,1)/pol2(0.65,Wright_bin,1);
sigma_wright[1] = sigma_wright[1]*getEpsL(2.445,W_ARR[Wright_bin],0.65)/getEpsL(2.445,W_ARR[Wright_bin],Q2gen);

//sigma_c2f
sigma_wright[2] = sigma_wright[2]*Func_q2_dep(Q2gen)/Func_q2_dep(0.65);
//sigma_wright[2] = sigma_wright[2]*(pol2(Q2gen,Wright_bin,2) -1.)/(pol2(0.65,Wright_bin,2) -1.);
//sigma_wright[2] = sigma_wright[2]*getEpsT(2.445,W_ARR[Wright_bin],0.65)/getEpsT(2.445,W_ARR[Wright_bin],Q2gen);
//sigma_wright[2] = sigma_wright[2]* func_sigma_t(Q2gen,Wright_bin)/SIGMA_T_INT_RIPANI[Wright_bin];

//sigma_cf
sigma_wright[4] = sigma_wright[4]*Func_q2_dep(Q2gen)/Func_q2_dep(0.65);
//sigma_wright[4] = sigma_wright[4]*(pol2(Q2gen,Wright_bin,4) -1.)/(pol2(0.65,Wright_bin,4) -1.);
//sigma_wright[4] = //sigma_wright[4]*sqrt(2*getEpsL(2.445,W_ARR[Wright_bin],0.65)*(getEpsT(2.445,W_ARR[Wright_bin],0.65)+1))/sqrt(2*getEpsL(2.445,W_ARR[Wright_bin],Q2gen)*(getEpsT(2.445,W_ARR[Wright_bin],Q2gen)+1));
//sigma_wright[4] = sigma_wright[4]*func_sigma_t(Q2gen,Wright_bin)/SIGMA_T_INT_RIPANI[Wright_bin];

//sigma_s2f
sigma_wright[3] = sigma_wright[3]*Func_q2_dep(Q2gen)/Func_q2_dep(0.65);

//sigma_sf
sigma_wright[5] = sigma_wright[5]*Func_q2_dep(Q2gen)/Func_q2_dep(0.65);

////sigma_l
sigma_wleft[1] = sigma_wleft[1]*pol2(Q2gen,Wleft_bin,1)/pol2(0.65,Wleft_bin,1);
sigma_wleft[1] = sigma_wleft[1]*getEpsL(2.445,W_ARR[Wleft_bin],0.65)/getEpsL(2.445,W_ARR[Wleft_bin],Q2gen);
//sigma_c2f
sigma_wleft[2] = sigma_wleft[2]*Func_q2_dep(Q2gen)/Func_q2_dep(0.65);
//sigma_wleft[2] = sigma_wleft[2]*(pol2(Q2gen,Wleft_bin,2) -1.)/(pol2(0.65,Wleft_bin,2) -1.);
//sigma_wleft[2] = sigma_wleft[2]*getEpsT(2.445,W_ARR[Wleft_bin],0.65)/getEpsT(2.445,W_ARR[Wleft_bin],Q2gen);
//sigma_wleft[2] = sigma_wleft[2]* func_sigma_t(Q2gen,Wleft_bin)/SIGMA_T_INT_RIPANI[Wleft_bin];

//sigma_cf
sigma_wleft[4] = sigma_wleft[4]*Func_q2_dep(Q2gen)/Func_q2_dep(0.65);
//sigma_wleft[4] = sigma_wleft[4]*(pol2(Q2gen,Wleft_bin,4) -1.)/(pol2(0.65,Wleft_bin,4) -1.);
//sigma_wleft[4] = sigma_wleft[4]*sqrt(2*getEpsL(2.445,W_ARR[Wleft_bin],0.65)*(getEpsT(2.445,W_ARR[Wleft_bin],0.65)+1))/sqrt(2*getEpsL(2.445,W_ARR[Wleft_bin],Q2gen)*(getEpsT(2.445,W_ARR[Wleft_bin],Q2gen)+1));
//sigma_wleft[4] = sigma_wleft[4]* func_sigma_t(Q2gen,Wleft_bin)/SIGMA_T_INT_RIPANI[Wleft_bin];

//sigma_s2f
sigma_wleft[3] = sigma_wleft[3]*Func_q2_dep(Q2gen)/Func_q2_dep(0.65);

//sigma_sf
sigma_wleft[5] = sigma_wleft[5]*Func_q2_dep(Q2gen)/Func_q2_dep(0.65);


//We are doing 1dim linear W-interpolation
for (Short_t i=0;i<6;i++){
sigma_final[i] = 1./fabs(W_ARR[Wright_bin]-W_ARR[Wleft_bin]);
sigma_final[i] = sigma_final[i]*(sigma_wright[i]*fabs(W_ARR[Wleft_bin]-Wgen)+sigma_wleft[i]*fabs(W_ARR[Wright_bin]-Wgen));
};


//We get explicitly different sigmas from the array
sigma_t_final = sigma_final[0];
sigma_l_final = sigma_final[1];
sigma_c2f_final = sigma_final[2];
sigma_s2f_final = sigma_final[3];
sigma_cf_final = sigma_final[4];
sigma_sf_final = sigma_final[5];
Float_t sigma_t_gol;


//Here we take Golovach photoproduction cross sections for Wright and Wleft and scale them with Q2  (only sigma_t)
if ((Wgen>=1.6125)&&(Wgen<=1.8125)) {
 
get_xsect_gol_datamod(Wgen, s12gen,s23gen, thetagen, alphagen, w_left_bin_gol, sigma_t_wright_gol,sigma_t_wleft_gol );
 
A_tmp[0] = 75.002141;
A_tmp[1] = 73.512912;
A_tmp[2] = 72.50197;
A_tmp[3] = 72.790675;
A_tmp[4] = 70.36;
A_tmp[5] = 65.2;
A_tmp[6] = 60.77;
A_tmp[7] = 58.02;
A_tmp[8] = 57.32;

//cout <<Wgen<<" "<< sigma_t_wright_gol << " "<<sigma_t_wleft_gol<<"\n";

//if ((Wgen>=1.7125)&&(Wgen<=1.8125)) {  
//sigma_t_wright_gol = sigma_t_wright_gol*func_sigma_t(Q2gen,Wright_bin)/func_sigma_t(0.005,Wright_bin);
//sigma_t_wleft_gol = sigma_t_wleft_gol*func_sigma_t(Q2gen,Wleft_bin)/func_sigma_t(0.005,Wleft_bin);
//};
//if ((Wgen>1.6125)&&(Wgen<=1.7125)) {
sigma_t_wright_gol = sigma_t_wright_gol*func_sigma_t(Q2gen,Wright_bin)/A_tmp[Wright_bin-8];
sigma_t_wleft_gol = sigma_t_wleft_gol*func_sigma_t(Q2gen,Wleft_bin)/A_tmp[Wleft_bin-8];
//};

//1dim W-interpolation of scalecd Golovach sigma_t
sigma_t_gol = 1./fabs(W_ARR[Wright_bin]-W_ARR[Wleft_bin]);
sigma_t_gol = sigma_t_gol*(sigma_t_wright_gol*fabs(W_ARR[Wleft_bin]-Wgen)+sigma_t_wleft_gol*fabs(W_ARR[Wright_bin]-Wgen));

//Mixing Ripani and Golovach cross sections
sigma_t_final = (Q2gen-0.0003)/0.65*sigma_t_final + (0.65-Q2gen)/0.65*sigma_t_gol;

};


};


