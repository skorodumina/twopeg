#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "get_xsect_fedotov.h"
#include "get_xsect_w16_18_lowq2_fit.h"
#include "get_xsect_ripani.h"
#include "interpol.h"
#include "interpol_fedotov.h"
#include "interpol_int.h"
using namespace std;

Float_t pol1_fed_0575_065 (Float_t  x, Short_t wbin, Short_t i) {
Float_t func;
 
//if (i==0) func =  -16.0533*x + 20.7038; 
//if (i==1) func =  34.7661*x + -12.3857; 
//if (i==2) func =  -2.03151*x + 0.273179;  
//if (i==4) func =  13.8965*x + -5.00247;  

if (i==0) func =  -28.2913*x +28.6584 ; 
if (i==1) func =  26.6543*x -7.113 ; 
if (i==2) func =  -1.07691*x -0.34731; 
if (i==4) func =  10.7093*x -2.93076;  

return func;
};

Float_t pol1_fed_065_095 (Float_t  x, Short_t wbin, Short_t i) {
Float_t func;
 
if (i==0) func =  -18.3052*x +22.1675; 
if (i==1) func = -2.26447*x + 11.6842; 
if (i==2) func = 1.8495*x -2.24948;  
if (i==4) func =  -6.70793*x +8.39044;  

return func;
};


Float_t pol1_fed_095_130 (Float_t  x, Short_t wbin, Short_t i) {
Float_t func;
 
if (i==0) func =  -8.76446*x +13.1038; 
if (i==1) func = -15.716*x + 24.4631; 
if (i==2) func = 0.846243*x -1.29638;  
if (i==4) func =  -4.27051*x +6.07489;  

return func;
};


//-------------------------------------------------------------------------
Float_t func_sigma_t_fed (Float_t  x, Short_t wbin) {
Float_t func = FIT_PARAM_SIGMA_T_FED[1][wbin]*x + FIT_PARAM_SIGMA_T_FED[0][wbin];
//Float_t func = pow((x+FIT_PARAM_SIGMA_T_FED[0][wbin]), FIT_PARAM_SIGMA_T_FED[1][wbin])*FIT_PARAM_SIGMA_T_FED[2][wbin]+FIT_PARAM_SIGMA_T_FED[3][wbin];
return func;
};

Float_t pol2_fed (Float_t  x, Short_t wbin, Short_t i) {
Float_t func;
if (i==1) func =  FIT_PARAM_SIGMA_L_FED[2][wbin]*x*x + FIT_PARAM_SIGMA_L_FED[1][wbin]*x + FIT_PARAM_SIGMA_L_FED[0][wbin]; 

return func;
};

//This subroutione performs joining between Ripani and Fedotov cross sections in the regions:
//REGION 1. W from 1.4125 to 1.5875 GeV in between the coverages
//REGION 2. (Wgen>=1.3125)&&(Wgen<1.4125)&&(Q2gen>=0.575)&&(Q2gen<=1.3) -- to the left from Ripani's coverage
//REGION 3. W from 1.315 to 1.5875 GeV and Q2 from around 0.25 to zero -- below Fedotov's coverage


void get_xsect_rip_fed_join(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final){

//REGION 1. W from 1.4125 to 1.5875 GeV in between the coverages
Float_t sigma_t_final_2, sigma_l_final_2, sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2;
Float_t sigma_wr_fed[6];
Float_t sigma_wl_fed[6];
Float_t sigma_wr_rip[6];
Float_t sigma_wl_rip[6];

Float_t sigma_wr_both[6];
Float_t sigma_wl_both[6];

Float_t sigma_wr_fed1[6];
Float_t sigma_wr_fed2[6];

Short_t Wleft_bin_fed,Wright_bin_fed, Wleft_bin_rip,Wright_bin_rip,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s12left_wright_bin_fed,s12right_wright_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,s23left_wright_bin_fed,s23right_wright_bin_fed,thetaleft_bin_fed ,thetaright_bin_fed, alphaleft_bin_fed,alpharight_bin_fed ;

Short_t s12left_wleft_bin_rip,s12right_wleft_bin_rip,s12left_wright_bin_rip,s12right_wright_bin_rip,s23left_wleft_bin_rip,s23right_wleft_bin_rip,s23left_wright_bin_rip,s23right_wright_bin_rip,thetaleft_bin_rip,thetaright_bin_rip,alphaleft_bin_rip,alpharight_bin_rip;

Short_t Q2_bin_fed_l,Q2_bin_fed_r ,Q2_bin_rip;

Float_t Q2_edge_r,Q2_edge_l, Q2_up,Q2_down ;

Float_t sigma_final[6];

Wleft_bin_fed = getWbin_fed(Wgen);
Wright_bin_fed = Wleft_bin_fed+1;


s12left_wleft_bin_fed = getsbin_fed(Wleft_bin_fed, s12gen, S12_ARR_FED[9][Wleft_bin_fed], S12_ARR_FED[0][Wleft_bin_fed]);
s12right_wleft_bin_fed = s12left_wleft_bin_fed +1;

s12left_wright_bin_fed = getsbin_fed(Wright_bin_fed, s12gen, S12_ARR_FED[9][Wright_bin_fed], S12_ARR_FED[0][Wright_bin_fed]);
s12right_wright_bin_fed = s12left_wright_bin_fed +1;

s23left_wleft_bin_fed = getsbin_fed(Wleft_bin_fed, s23gen, S23_ARR_FED[9][Wleft_bin_fed], S23_ARR_FED[0][Wleft_bin_fed]);
s23right_wleft_bin_fed = s23left_wleft_bin_fed +1;

s23left_wright_bin_fed = getsbin_fed(Wright_bin_fed, s23gen, S23_ARR_FED[9][Wright_bin_fed], S23_ARR_FED[0][Wright_bin_fed]);
s23right_wright_bin_fed = s23left_wright_bin_fed +1;

thetaleft_bin_fed = getanglebin_fed(thetagen,THETA_ARR_FED[7]);
thetaright_bin_fed = thetaleft_bin_fed+1;

alphaleft_bin_fed = getanglebin_fed(alphagen,ALPHA_ARR_FED[7]);
alpharight_bin_fed = alphaleft_bin_fed+1;

//-----------------------------------------


if (Wgen >= 1.4125){

Wleft_bin_rip = getWbin(Wgen);
Wright_bin_rip = Wleft_bin_rip+1;

s12left_wleft_bin_rip = getsbin(Wleft_bin_rip, s12gen, S12_ARR[11][Wleft_bin_rip], S12_ARR[0][Wleft_bin_rip]);
s12right_wleft_bin_rip = s12left_wleft_bin_rip +1;

s12left_wright_bin_rip = getsbin(Wright_bin_rip, s12gen, S12_ARR[11][Wright_bin_rip], S12_ARR[0][Wright_bin_rip]);
s12right_wright_bin_rip = s12left_wright_bin_rip +1;


s23left_wleft_bin_rip = getsbin(Wleft_bin_rip, s23gen, S23_ARR[11][Wleft_bin_rip], S23_ARR[0][Wleft_bin_rip]);
s23right_wleft_bin_rip = s23left_wleft_bin_rip +1;

s23left_wright_bin_rip = getsbin(Wright_bin_rip, s23gen, S23_ARR[11][Wright_bin_rip], S23_ARR[0][Wright_bin_rip]);
s23right_wright_bin_rip = s23left_wright_bin_rip +1;


thetaleft_bin_rip = getanglebin(thetagen,THETA_ARR[5]);
thetaright_bin_rip = thetaleft_bin_rip+1;

alphaleft_bin_rip = getanglebin(alphagen,ALPHA_ARR[5]);
alpharight_bin_rip = alphaleft_bin_rip+1;
Q2_bin_rip = 0;


if ((Wgen>=1.4125)&&(Wgen<=1.4375)&&(Q2gen>=0.575)&&(Q2gen<=0.65)){
Q2_bin_fed_l =6;
Q2_bin_fed_r =6;
};

if ((Wgen>=1.4375)&&(Wgen<=1.4625)&&(Q2gen>=0.525)&&(Q2gen<=0.65)){
Q2_bin_fed_l =6;
Q2_bin_fed_r =5;
};

if ((Wgen>=1.4625)&&(Wgen<=1.4875)&&(Q2gen>=0.525)&&(Q2gen<=0.65)){
Q2_bin_fed_l =5;
Q2_bin_fed_r =5;
};

if ((Wgen>=1.4875)&&(Wgen<=1.5125)&&(Q2gen>=0.425)&&(Q2gen<=0.65)){
Q2_bin_fed_l =5;
Q2_bin_fed_r =3;
};

if ((Wgen>=1.5125)&&(Wgen<=1.5375)&&(Q2gen>=0.325)&&(Q2gen<=0.65)){
Q2_bin_fed_l =3;
Q2_bin_fed_r =2;
};

if ((Wgen>=1.5375)&&(Wgen<=1.5625)&&(Q2gen>=0.275)&&(Q2gen<=0.65)){
Q2_bin_fed_l =2;
Q2_bin_fed_r =1;
};

if ((Wgen>=1.5625)&&(Wgen<=1.5875)&&(Q2gen>=0.225)&&(Q2gen<=0.65)){
Q2_bin_fed_l =1;
Q2_bin_fed_r =0;
};

for (Short_t i=0;i<6;i++){

interpol_fedotov(4,Q2_bin_fed_r,Wright_bin_fed,s12left_wright_bin_fed,s12right_wright_bin_fed,s23left_wright_bin_fed,s23right_wright_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed[i],i);

interpol_fedotov(4,Q2_bin_fed_l,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);

interpol(4,Q2_bin_rip,Wright_bin_rip,s12left_wright_bin_rip,s12right_wright_bin_rip,s23left_wright_bin_rip,s23right_wright_bin_rip,thetaleft_bin_rip,thetaright_bin_rip,alphaleft_bin_rip,alpharight_bin_rip,s12gen,s23gen,thetagen,alphagen,sigma_wr_rip[i],i);

interpol(4,Q2_bin_rip,Wleft_bin_rip,s12left_wleft_bin_rip,s12right_wleft_bin_rip,s23left_wleft_bin_rip,s23right_wleft_bin_rip,thetaleft_bin_rip,thetaright_bin_rip,alphaleft_bin_rip,alpharight_bin_rip,s12gen,s23gen,thetagen,alphagen,sigma_wl_rip[i],i);

Q2_up = 0.65;
Q2_down = Q2_ARR_FED[Q2_bin_fed_l];

if (((Wgen>=1.4375)&&(Wgen<=1.4625)&&(Q2gen>=0.525)&&(Q2gen<=0.575))||((Wgen>=1.5125)&&(Wgen<=1.5375)&&(Q2gen>=0.325)&&(Q2gen<=0.425))||((Wgen>=1.5375)&&(Wgen<=1.5625)&&(Q2gen>=0.275)&&(Q2gen<=0.325))||((Wgen>=1.5625)&&(Wgen<=1.5875)&&(Q2gen>=0.225)&&(Q2gen<=0.275))){

interpol_fedotov(4,Q2_bin_fed_l,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_rip[i],i);

interpol_fedotov(4,Q2_bin_fed_r,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);
Q2_up = Q2_ARR_FED[Q2_bin_fed_l];
Q2_down = Q2_ARR_FED[Q2_bin_fed_r];
};

if (((Wgen>=1.4875)&&(Wgen<=1.5125)&&(Q2gen>=0.425)&&(Q2gen<=0.475))){

interpol_fedotov(4,4,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_rip[i],i);

interpol_fedotov(4,3,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);
Q2_up = Q2_ARR_FED[4];
Q2_down = Q2_ARR_FED[3];
};

if (((Wgen>=1.4875)&&(Wgen<=1.5125)&&(Q2gen>=0.475)&&(Q2gen<=0.525))){

interpol_fedotov(4,5,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_rip[i],i);

interpol_fedotov(4,4,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);
Q2_up = Q2_ARR_FED[5];
Q2_down = Q2_ARR_FED[4];
};


Float_t factor;

sigma_wl_both[i] = 1./fabs(Q2_up - Q2_down);
sigma_wl_both[i] = sigma_wl_both[i]*(sigma_wl_rip[i]*fabs(Q2_down-Q2gen) + sigma_wl_fed[i]*fabs(Q2_up-Q2gen));

sigma_wr_both[i] = 1./fabs(0.65 - Q2_ARR_FED[Q2_bin_fed_r]);
sigma_wr_both[i] =sigma_wr_both[i]*(sigma_wr_rip[i]*fabs(Q2_ARR_FED[Q2_bin_fed_r]-Q2gen) + sigma_wr_fed[i]*fabs(0.65-Q2gen));

sigma_final[i] = 1./0.025;
sigma_final[i] = sigma_final[i]*(sigma_wr_both[i]*fabs(W_ARR[Wleft_bin_rip]-Wgen)+sigma_wl_both[i]*fabs(W_ARR[Wright_bin_rip]-Wgen));

};
};//end if (Wgen>=1.4125)


//----------------------------
//REGION 2. (Wgen>=1.3125)&&(Wgen<1.4125)&&(Q2gen>=0.575)&&(Q2gen<=1.3) -- to the left from Ripani's coverage

if ((Wgen>=1.3125)&&(Wgen<1.4125)&&(Q2gen>=0.575)&&(Q2gen<=1.3)){

for (Short_t i=0;i<6;i++){
interpol_fedotov(4,6,Wright_bin_fed,s12left_wright_bin_fed,s12right_wright_bin_fed,s23left_wright_bin_fed,s23right_wright_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed[i],i);

interpol_fedotov(4,6,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);
};

//Here we Q2-scale sigma_t, sigma_t, sigma_c2f and sigma_cf with the fit_functions (corresponding pol1) for Wright_bin_fed and Wleft_bin_fed
if ((Q2gen>=0.575)&&(Q2gen<=0.65)){
//sigma_t
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_0575_065(Q2gen,Wright_bin_fed,0)/pol1_fed_0575_065(0.575,Wright_bin_fed,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_0575_065(Q2gen,Wleft_bin_fed,0)/pol1_fed_0575_065(0.575,Wleft_bin_fed,0);
////sigma_l
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_0575_065(Q2gen,Wright_bin_fed,1)/pol1_fed_0575_065(0.575,Wright_bin_fed,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_0575_065(Q2gen,Wleft_bin_fed,1)/pol1_fed_0575_065(0.575,Wleft_bin_fed,1);
//sigma_c2f
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_0575_065(Q2gen,Wright_bin_fed,2)/pol1_fed_0575_065(0.575,Wright_bin_fed,2);
sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_0575_065(Q2gen,Wleft_bin_fed,2)/pol1_fed_0575_065(0.575,Wleft_bin_fed,2);
//sigma_cf
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_0575_065(Q2gen,Wright_bin_fed,4)/pol1_fed_0575_065(0.575,Wright_bin_fed,4);
sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_0575_065(Q2gen,Wleft_bin_fed,4)/pol1_fed_0575_065(0.575,Wleft_bin_fed,4);

};

if ((Q2gen>=0.65)&&(Q2gen<=0.95)){
//sigma_t
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_0575_065(0.65,Wright_bin_fed,0)/pol1_fed_0575_065(0.575,Wright_bin_fed,0);
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_065_095(Q2gen,Wleft_bin_fed,0)/pol1_fed_065_095(0.65,Wleft_bin_fed,0);

sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_0575_065(0.65,Wleft_bin_fed,0)/pol1_fed_0575_065(0.575,Wleft_bin_fed,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_065_095(Q2gen,Wleft_bin_fed,0)/pol1_fed_065_095(0.65,Wleft_bin_fed,0);

////sigma_l
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_0575_065(0.65,Wright_bin_fed,1)/pol1_fed_0575_065(0.575,Wright_bin_fed,1);
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_065_095(Q2gen,Wright_bin_fed,1)/pol1_fed_065_095(0.65,Wright_bin_fed,1);

sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_0575_065(0.65,Wleft_bin_fed,1)/pol1_fed_0575_065(0.575,Wleft_bin_fed,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_065_095(Q2gen,Wleft_bin_fed,1)/pol1_fed_065_095(0.65,Wleft_bin_fed,1);

//sigma_c2f
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_0575_065(0.65,Wright_bin_fed,2)/pol1_fed_0575_065(0.575,Wright_bin_fed,2);
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_065_095(Q2gen,Wright_bin_fed,2)/pol1_fed_065_095(0.65,Wright_bin_fed,2);

sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_0575_065(0.65,Wleft_bin_fed,2)/pol1_fed_0575_065(0.575,Wleft_bin_fed,2);
sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_065_095(Q2gen,Wleft_bin_fed,2)/pol1_fed_065_095(0.65,Wleft_bin_fed,2);

//sigma_cf
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_0575_065(0.65,Wright_bin_fed,4)/pol1_fed_0575_065(0.575,Wright_bin_fed,4);
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_065_095(Q2gen,Wright_bin_fed,4)/pol1_fed_065_095(0.65,Wright_bin_fed,4);

sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_0575_065(0.65,Wleft_bin_fed,4)/pol1_fed_0575_065(0.575,Wleft_bin_fed,4);
sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_065_095(Q2gen,Wleft_bin_fed,4)/pol1_fed_065_095(0.65,Wleft_bin_fed,4);

};


if ((Q2gen>=0.95)&&(Q2gen<=1.3)){

//sigma_t
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_0575_065(0.65,Wright_bin_fed,0)/pol1_fed_0575_065(0.575,Wright_bin_fed,0);
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_065_095(0.95,Wleft_bin_fed,0)/pol1_fed_065_095(0.65,Wleft_bin_fed,0);
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_095_130(Q2gen,Wright_bin_fed,0)/pol1_fed_095_130(0.95,Wright_bin_fed,0);

sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_0575_065(0.65,Wleft_bin_fed,0)/pol1_fed_0575_065(0.575,Wleft_bin_fed,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_065_095(0.95,Wleft_bin_fed,0)/pol1_fed_065_095(0.65,Wleft_bin_fed,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_095_130(Q2gen,Wleft_bin_fed,0)/pol1_fed_095_130(0.95,Wleft_bin_fed,0);

////sigma_l
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_0575_065(0.65,Wright_bin_fed,1)/pol1_fed_0575_065(0.575,Wright_bin_fed,1);
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_065_095(0.95,Wright_bin_fed,1)/pol1_fed_065_095(0.65,Wright_bin_fed,1);
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_095_130(Q2gen,Wright_bin_fed,1)/pol1_fed_095_130(0.95,Wright_bin_fed,1);

sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_0575_065(0.65,Wleft_bin_fed,1)/pol1_fed_0575_065(0.575,Wleft_bin_fed,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_065_095(0.95,Wleft_bin_fed,1)/pol1_fed_065_095(0.65,Wleft_bin_fed,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_095_130(Q2gen,Wleft_bin_fed,1)/pol1_fed_095_130(0.95,Wleft_bin_fed,1);
//sigma_c2f
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_0575_065(0.65,Wright_bin_fed,2)/pol1_fed_0575_065(0.575,Wright_bin_fed,2);
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_065_095(0.95,Wright_bin_fed,2)/pol1_fed_065_095(0.65,Wright_bin_fed,2);
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_095_130(Q2gen,Wright_bin_fed,2)/pol1_fed_095_130(0.95,Wright_bin_fed,2);

sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_0575_065(0.65,Wleft_bin_fed,2)/pol1_fed_0575_065(0.575,Wleft_bin_fed,2);
sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_065_095(0.95,Wleft_bin_fed,2)/pol1_fed_065_095(0.65,Wleft_bin_fed,2);
sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_095_130(Q2gen,Wleft_bin_fed,2)/pol1_fed_095_130(0.95,Wleft_bin_fed,2);
//sigma_cf
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_0575_065(0.65,Wright_bin_fed,4)/pol1_fed_0575_065(0.575,Wright_bin_fed,4);
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_065_095(0.95,Wright_bin_fed,4)/pol1_fed_065_095(0.65,Wright_bin_fed,4);
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_095_130(Q2gen,Wright_bin_fed,4)/pol1_fed_095_130(0.95,Wright_bin_fed,4);

sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_0575_065(0.65,Wleft_bin_fed,4)/pol1_fed_0575_065(0.575,Wleft_bin_fed,4);
sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_065_095(0.95,Wleft_bin_fed,4)/pol1_fed_065_095(0.65,Wleft_bin_fed,4);
sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_095_130(Q2gen,Wleft_bin_fed,4)/pol1_fed_095_130(0.95,Wleft_bin_fed,4);

};


//sigma_s2f
sigma_wr_fed[3] = sigma_wr_fed[3]*Func_q2_dep(Q2gen)/Func_q2_dep(0.575);
sigma_wl_fed[3] = sigma_wl_fed[3]*Func_q2_dep(Q2gen)/Func_q2_dep(0.575);

//sigma_sf
sigma_wr_fed[5] = sigma_wr_fed[5]*Func_q2_dep(Q2gen)/Func_q2_dep(0.575);
sigma_wl_fed[5] = sigma_wl_fed[5]*Func_q2_dep(Q2gen)/Func_q2_dep(0.575);


//We are doing 1dim linear W-interpolation
for (Short_t i=0;i<6;i++){
sigma_final[i] = 1./fabs(W_ARR_FED[Wright_bin_fed]-W_ARR_FED[Wleft_bin_fed]);
sigma_final[i] = sigma_final[i]*(sigma_wr_fed[i]*fabs(W_ARR_FED[Wleft_bin_fed]-Wgen)+sigma_wl_fed[i]*fabs(W_ARR_FED[Wright_bin_fed]-Wgen));
};

//sigma_final[3] = 0.;
//sigma_final[5] = 0.;

};//end if ((Wgen>=1.3125)&&(Wgen<1.4125)&&(Q2gen>=0.575)&&(Q2gen<=1.3))


//--------------------------------------------------
//REGION 3. W from 1.315 to 1.5875 GeV and Q2 from around 0.25 to zero -- below Fedotov's coverage

if (((Wgen>=1.3125)&&(Wgen<1.5125)&&(Q2gen>=0.00002)&&(Q2gen<=0.275))||((Wgen>=1.5125)&&(Wgen<=1.5875)&&(Q2gen>=0.00002)&&(Q2gen<=0.225))){

for (Short_t i=0;i<6;i++){

if ((Wgen>=1.3125)&&(Wgen<=1.5125)&&(Q2gen>=0.00002)&&(Q2gen<=0.275)){
Q2_edge_r = 0.275;
Q2_edge_l = 0.275;
interpol_fedotov(4,1,Wright_bin_fed,s12left_wright_bin_fed,s12right_wright_bin_fed,s23left_wright_bin_fed,s23right_wright_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed[i],i);

interpol_fedotov(4,1,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);
};


if ((Wgen>=1.5125)&&(Wgen<=1.5875)&&(Q2gen>=0.00002)&&(Q2gen<=0.225)){
Q2_edge_l = 0.225;
Q2_edge_r = 0.225;
interpol_fedotov(4,0,Wright_bin_fed,s12left_wright_bin_fed,s12right_wright_bin_fed,s23left_wright_bin_fed,s23right_wright_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed[i],i);

interpol_fedotov(4,0,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);
};


if ((Wgen>=1.4875)&&(Wgen<1.5125)&&(Q2gen>=0.00002)&&(Q2gen<=0.275)){
Q2_edge_l = 0.275;
Q2_edge_r = 0.225;

interpol_fedotov(4,1,Wleft_bin_fed,s12left_wleft_bin_fed,s12right_wleft_bin_fed,s23left_wleft_bin_fed,s23right_wleft_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);

interpol_fedotov(4,0,Wright_bin_fed,s12left_wright_bin_fed,s12right_wright_bin_fed,s23left_wright_bin_fed,s23right_wright_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed[i],i);
//cout<< Wleft_bin_fed << "\n";


if ((Wgen>=1.4875)&&(Wgen<1.5125)&&(Q2gen>=0.225)&&(Q2gen<=0.275)){

interpol_fedotov(4,0,Wright_bin_fed,s12left_wright_bin_fed,s12right_wright_bin_fed,s23left_wright_bin_fed,s23right_wright_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed1[i],i);
//cout<< Wleft_bin_fed << "\n";

interpol_fedotov(4,1,Wright_bin_fed,s12left_wright_bin_fed,s12right_wright_bin_fed,s23left_wright_bin_fed,s23right_wright_bin_fed,thetaleft_bin_fed,thetaright_bin_fed,alphaleft_bin_fed,alpharight_bin_fed,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed2[i],i);

sigma_wr_fed[i] = 1./0.05;
sigma_wr_fed[i] = sigma_wr_fed[i]*(sigma_wr_fed1[i]*fabs(0.275-Q2gen) + sigma_wr_fed2[i]*fabs(0.225-Q2gen));
};

};
};//end for i

//Here we Q2-scale sigma_t, sigma_t, sigma_c2f and sigma_cf with the fit_functions (func_sigma_t and corresponding pol2) for Wright_bin and Wleft_bin

if ((Wgen>=1.3125)&&(Wgen<1.3625)&&(Q2gen>=0.000001)&&(Q2gen<=0.275)){
//sigma_t
sigma_wr_fed[0] = sigma_wr_fed[0]*func_sigma_t_fed(Q2gen,0)/func_sigma_t_fed(Q2_edge_r,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*func_sigma_t_fed(Q2gen,0)/func_sigma_t_fed(Q2_edge_l,0);

//sigma_l
sigma_wr_fed[1] = sigma_wr_fed[1]*pol2_fed(Q2gen,0,1)/pol2_fed(Q2_edge_r,0,1);
sigma_wr_fed[1] = sigma_wr_fed[1]*getEpsL(1.515,W_ARR_FED[2],Q2_edge_r)/getEpsL(1.515,W_ARR_FED[2],Q2gen);

////sigma_l
//cout << "qqq2 "<< sigma_wl_fed[1] << "\n";
sigma_wl_fed[1] = sigma_wl_fed[1]*pol2_fed(Q2gen,0,1)/pol2_fed(Q2_edge_l,0,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*getEpsL(1.515,W_ARR_FED[2],Q2_edge_l)/getEpsL(1.515,W_ARR_FED[2],Q2gen);
};

if (((Wgen>=1.3625)&&(Wgen<=1.4875)&&(Q2gen>=0.000001)&&(Q2gen<=0.275))||((Wgen>=1.4875)&&(Wgen<1.5875)&&(Q2gen>=0.000001)&&(Q2gen<=0.225))){
//sigma_t
sigma_wr_fed[0] = sigma_wr_fed[0]*func_sigma_t_fed(Q2gen,Wright_bin_fed - 2)/func_sigma_t_fed(Q2_edge_r,Wright_bin_fed - 2);
sigma_wl_fed[0] = sigma_wl_fed[0]*func_sigma_t_fed(Q2gen,Wleft_bin_fed - 2)/func_sigma_t_fed(Q2_edge_l,Wleft_bin_fed - 2);

//sigma_l
sigma_wr_fed[1] = sigma_wr_fed[1]*pol2_fed(Q2gen,Wright_bin_fed - 2,1)/pol2_fed(Q2_edge_r,Wright_bin_fed - 2,1);
sigma_wr_fed[1] = sigma_wr_fed[1]*getEpsL(1.515,W_ARR_FED[Wright_bin_fed],Q2_edge_r)/getEpsL(1.515,W_ARR_FED[Wright_bin_fed],Q2gen);

//sigma_l
//cout << "qqq2 "<< sigma_wl_fed[1] << "\n";
sigma_wl_fed[1] = sigma_wl_fed[1]*pol2_fed(Q2gen,Wleft_bin_fed - 2,1)/pol2_fed(Q2_edge_l,Wleft_bin_fed - 2,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*getEpsL(1.515,W_ARR_FED[Wleft_bin_fed],Q2_edge_l)/getEpsL(1.515,W_ARR_FED[Wleft_bin_fed],Q2gen);

};

//sigma_c2f
sigma_wr_fed[2] = sigma_wr_fed[2]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_r); 
sigma_wl_fed[2] = sigma_wl_fed[2]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_l); 

//sigma_s2f
sigma_wr_fed[3] = sigma_wr_fed[3]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_r); 
sigma_wl_fed[3] = sigma_wl_fed[3]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_l); 

//sigma_cf
sigma_wr_fed[4] = sigma_wr_fed[4]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_r); 
sigma_wl_fed[4] = sigma_wl_fed[4]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_l); 

//sigma_sf
sigma_wr_fed[5] = sigma_wr_fed[5]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_r); 
sigma_wl_fed[5] = sigma_wl_fed[5]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_l);




if ((Wgen>=1.4875)&&(Wgen<1.5125)&&(Q2gen>=0.225)&&(Q2gen<=0.275)){

sigma_wl_fed[0] = sigma_wl_fed[0]*func_sigma_t_fed(Q2gen,Wleft_bin_fed - 2)/func_sigma_t_fed(Q2_edge_l,Wleft_bin_fed - 2);

sigma_wl_fed[1] = sigma_wl_fed[1]*pol2_fed(Q2gen,Wleft_bin_fed - 2,1)/pol2_fed(Q2_edge_l,Wleft_bin_fed - 2,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*getEpsL(1.515,W_ARR_FED[Wleft_bin_fed],Q2_edge_l)/getEpsL(1.515,W_ARR_FED[Wleft_bin_fed],Q2gen);

//sigma_c2f
sigma_wl_fed[2] = sigma_wl_fed[2]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_l); 
//sigma_s2f
sigma_wl_fed[3] = sigma_wl_fed[3]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_l); 
//sigma_cf
sigma_wl_fed[4] = sigma_wl_fed[4]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_l); 
//sigma_sf
sigma_wl_fed[5] = sigma_wl_fed[5]*Func_q2_dep(Q2gen)/Func_q2_dep(Q2_edge_l);
};


//We are doing 1dim linear W-interpolation
for (Short_t i=0;i<6;i++){
sigma_final[i] = 1./fabs(W_ARR_FED[Wright_bin_fed]-W_ARR_FED[Wleft_bin_fed]);
sigma_final[i] = sigma_final[i]*(sigma_wr_fed[i]*fabs(W_ARR_FED[Wleft_bin_fed]-Wgen)+sigma_wl_fed[i]*fabs(W_ARR_FED[Wright_bin_fed]-Wgen));
//cout <<i<< " "<< sigma_final[i] << sigma_wr_fed[i] <<" "<< sigma_wl_fed[i] <<"\n";
};

};


//We get explicitly different sigmas from the array
sigma_t_final = sigma_final[0];
sigma_l_final = sigma_final[1];
sigma_c2f_final = sigma_final[2];
sigma_s2f_final = sigma_final[3];
sigma_cf_final = sigma_final[4];
sigma_sf_final = sigma_final[5];


};
