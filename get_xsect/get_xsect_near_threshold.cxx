#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "interpol_fedotov_thresh.h"
#include "interpol_fedotov.h"
#include "get_xsect_fedotov.h"
#include "get_xsect_rip_fed_join.h"
#include "get_xsect_14_18_lowq2_fit.h"
#include "interpol_int.h"
using namespace std;

Short_t getWbin_fed_thresh (Float_t W) {

Short_t bin;
if ((W>=1.2375)&&(W<=1.2625)) bin = 0;
if ((W>=1.2625)&&(W<=1.2875)) bin = 1;
if ((W>=1.2875)&&(W<=1.3125)) bin = 2;

if ((W<1.2375)||(W>1.3125)) {
cout << "Error, wrong W range, threshold "<< "\n";
bin = -100;
};

return bin;
};



void get_xsect_near_threshold(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final ){
Float_t sigma_wr_fed[6];
Float_t sigma_wl_fed[6];

Float_t sigma_wr_both[6];
Float_t sigma_wl_both[6];

Short_t Wleft_bin = getWbin_fed_thresh(Wgen);
Short_t Wright_bin = Wleft_bin+1;

//cout << Wgen<<" "<< Wleft_bin<< " "<< Wright_bin<< endl;

Short_t Q2left_bin,Q2right_bin ;


Short_t  s12left_wleft_bin,s12right_wleft_bin,s12left_wright_bin,s12right_wright_bin,s23left_wleft_bin ,s23right_wleft_bin,s23left_wright_bin,s23right_wright_bin;

 s12left_wleft_bin = getsbin_fed(Wleft_bin, s12gen, S12_ARR_FED_THRESH[9][Wleft_bin], S12_ARR_FED_THRESH[0][Wleft_bin]);
 s12right_wleft_bin = s12left_wleft_bin +1;



 s23left_wleft_bin = getsbin_fed(Wleft_bin, s23gen, S23_ARR_FED_THRESH[9][Wleft_bin], S23_ARR_FED_THRESH[0][Wleft_bin]);
 s23right_wleft_bin = s23left_wleft_bin +1;




Short_t thetaleft_bin = getanglebin_fed(thetagen,THETA_ARR_FED[7]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin_fed(alphagen,ALPHA_ARR_FED[7]);
Short_t alpharight_bin = alphaleft_bin+1;

Float_t sigma_final[6];


if ((Wgen>=1.2375)&&(Wgen<=1.3125)&&(Q2gen>0.275)&&(Q2gen<0.575)){

if ((Wgen>=1.2375)&&(Wgen<=1.2875)){

 Q2left_bin = getQ2bin_fed(Q2gen);
 Q2right_bin = Q2left_bin+1;


 s12left_wright_bin = getsbin_fed(Wright_bin, s12gen, S12_ARR_FED_THRESH[9][Wright_bin], S12_ARR_FED_THRESH[0][Wright_bin]);
 s12right_wright_bin = s12left_wright_bin +1;
 
s23left_wright_bin = getsbin_fed(Wright_bin, s23gen, S23_ARR_FED_THRESH[9][Wright_bin], S23_ARR_FED_THRESH[0][Wright_bin]);
 s23right_wright_bin = s23left_wright_bin +1;

for (Short_t i=0;i<6;i++){
interpol_fedotov_thresh(4,Q2left_bin,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_q2left[i],i);

interpol_fedotov_thresh(4,Q2right_bin,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_q2right[i],i);

interpol_fedotov_thresh(4,Q2right_bin,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_q2right[i],i);

interpol_fedotov_thresh(4,Q2left_bin,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_q2left[i],i);

interpol_fedotov_thresh(2,0,0, Wleft_bin, Wright_bin, Q2left_bin, Q2right_bin, 0,  0, 0,  0, Wgen, Q2gen, 0., 0.,sigma_final[i], i);
};
};


if ((Wgen>=1.2875)&&(Wgen<=1.3125)){
 Q2left_bin = getQ2bin_fed(Q2gen);
 Q2right_bin = Q2left_bin+1;
 s12left_wright_bin = getsbin_fed(0, s12gen, S12_ARR_FED[9][0], S12_ARR_FED[0][0]);
 s12right_wright_bin = s12left_wright_bin +1;
 
 s23left_wright_bin = getsbin_fed(0, s23gen, S23_ARR_FED[9][0], S23_ARR_FED[0][0]);
 s23right_wright_bin = s23left_wright_bin +1;
 
for (Short_t i=0;i<6;i++){
interpol_fedotov(4,Q2left_bin,0,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_q2left[i],i);

interpol_fedotov(4,Q2right_bin,0,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_q2right[i],i);

interpol_fedotov_thresh(4,Q2right_bin,2,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_q2right[i],i);

interpol_fedotov_thresh(4,Q2left_bin,2,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_q2left[i],i);


sigma_wl_both[i] = 1./fabs(Q2_ARR_FED[Q2left_bin]-Q2_ARR_FED[Q2right_bin]);
sigma_wl_both[i] = sigma_wl_both[i]*(sigma_wleft_q2right[i]*fabs(Q2_ARR_FED[Q2left_bin]-Q2gen) + sigma_wleft_q2left[i]*fabs(Q2_ARR_FED[Q2right_bin]-Q2gen));

sigma_wr_both[i] = 1./fabs(Q2_ARR_FED[Q2left_bin]-Q2_ARR_FED[Q2right_bin]);
sigma_wr_both[i] =sigma_wr_both[i]*(sigma_wright_q2right[i]*fabs(Q2_ARR_FED[Q2left_bin]-Q2gen) + sigma_wright_q2left[i]*fabs(Q2_ARR_FED[Q2right_bin]-Q2gen));

sigma_final[i] = 1./0.025;
sigma_final[i] = sigma_final[i]*(sigma_wr_both[i]*fabs(1.2875-Wgen)+sigma_wl_both[i]*fabs(1.3125-Wgen));
};
};


};//end   ((Wgen>=1.2375)&&(Wgen<=1.3125)&&(Q2gen>0.275)&&(Q2gen<0.575)){

//---------------------------------------------------------------------------------------------------------------------


if ((Wgen>=1.2375)&&(Wgen<=1.3125)&&(Q2gen>=0.0002)&&(Q2gen<=0.275)){

if ((Wgen>=1.2375)&&(Wgen<=1.2875)){
for (Short_t i=0;i<6;i++){
s12left_wright_bin = getsbin_fed(Wright_bin, s12gen, S12_ARR_FED_THRESH[9][Wright_bin], S12_ARR_FED_THRESH[0][Wright_bin]);
 s12right_wright_bin = s12left_wright_bin +1;
 
s23left_wright_bin = getsbin_fed(Wright_bin, s23gen, S23_ARR_FED_THRESH[9][Wright_bin], S23_ARR_FED_THRESH[0][Wright_bin]);
 s23right_wright_bin = s23left_wright_bin +1;
 
interpol_fedotov_thresh(4,1,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed[i],i);

interpol_fedotov_thresh(4,1,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);

};
};


if ((Wgen>=1.2875)&&(Wgen<=1.3125)){
s12left_wright_bin = getsbin_fed(0, s12gen, S12_ARR_FED[9][0], S12_ARR_FED[0][0]);
 s12right_wright_bin = s12left_wright_bin +1;
 
 s23left_wright_bin = getsbin_fed(0, s23gen, S23_ARR_FED[9][0], S23_ARR_FED[0][0]);
 s23right_wright_bin = s23left_wright_bin +1;

for (Short_t i=0;i<6;i++){
interpol_fedotov(4,1,0,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed[i],i);

interpol_fedotov_thresh(4,1,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);
};


};

sigma_wr_fed[0] = sigma_wr_fed[0]*func_sigma_t_fed(Q2gen,0)/func_sigma_t_fed(0.275,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*func_sigma_t_fed(Q2gen,0)/func_sigma_t_fed(0.275,0);

//sigma_l
sigma_wr_fed[1] = sigma_wr_fed[1]*pol2_fed(Q2gen,0,1)/pol2_fed(0.275,0,1);
sigma_wr_fed[1] = sigma_wr_fed[1]*getEpsL(1.515,W_ARR_FED[2],0.275)/getEpsL(1.515,W_ARR_FED[2],Q2gen);

////sigma_l
//cout << "qqq2 "<< sigma_wl_fed[1] << "\n";
sigma_wl_fed[1] = sigma_wl_fed[1]*pol2_fed(Q2gen,0,1)/pol2_fed(0.275,0,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*getEpsL(1.515,W_ARR_FED[2],0.275)/getEpsL(1.515,W_ARR_FED[2],Q2gen);

//sigma_c2f
sigma_wr_fed[2] = sigma_wr_fed[2]*Func_q2_dep(Q2gen)/Func_q2_dep(0.275); 
sigma_wl_fed[2] = sigma_wl_fed[2]*Func_q2_dep(Q2gen)/Func_q2_dep(0.275); 

//sigma_s2f
sigma_wr_fed[3] = sigma_wr_fed[3]*Func_q2_dep(Q2gen)/Func_q2_dep(0.275); 
sigma_wl_fed[3] = sigma_wl_fed[3]*Func_q2_dep(Q2gen)/Func_q2_dep(0.275); 

//sigma_cf
sigma_wr_fed[4] = sigma_wr_fed[4]*Func_q2_dep(Q2gen)/Func_q2_dep(0.275); 
sigma_wl_fed[4] = sigma_wl_fed[4]*Func_q2_dep(Q2gen)/Func_q2_dep(0.275); 

//sigma_sf
sigma_wr_fed[5] = sigma_wr_fed[5]*Func_q2_dep(Q2gen)/Func_q2_dep(0.275); 
sigma_wl_fed[5] = sigma_wl_fed[5]*Func_q2_dep(Q2gen)/Func_q2_dep(0.275);

for (Short_t i=0;i<6;i++){
sigma_final[i] = 1./0.025;
if ((Wgen>=1.2875)&&(Wgen<=1.3125)) sigma_final[i] = sigma_final[i]*(sigma_wr_fed[i]*fabs(1.2875-Wgen)+sigma_wl_fed[i]*fabs(1.3125-Wgen));

if ((Wgen>=1.2375)&&(Wgen<=1.2875)) sigma_final[i] = sigma_final[i]*(sigma_wr_fed[i]*fabs(W_ARR_FED_THRESH[Wleft_bin]-Wgen)+sigma_wl_fed[i]*fabs(W_ARR_FED_THRESH[Wright_bin]-Wgen));
};




};// end if ((Wgen>=1.2375)&&(Wgen<=1.3125)&&(Q2gen>0.002)&&(Q2gen<0.275)){

//-------------------------------------------------------------------------------


if ((Wgen>=1.2375)&&(Wgen<1.3125)&&(Q2gen>=0.575)&&(Q2gen<=1.3)){



if ((Wgen>=1.2375)&&(Wgen<=1.2875)){
for (Short_t i=0;i<6;i++){
s12left_wright_bin = getsbin_fed(Wright_bin, s12gen, S12_ARR_FED_THRESH[9][Wright_bin], S12_ARR_FED_THRESH[0][Wright_bin]);
 s12right_wright_bin = s12left_wright_bin +1;
 
s23left_wright_bin = getsbin_fed(Wright_bin, s23gen, S23_ARR_FED_THRESH[9][Wright_bin], S23_ARR_FED_THRESH[0][Wright_bin]);
 s23right_wright_bin = s23left_wright_bin +1;

interpol_fedotov_thresh(4,6,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed[i],i);

interpol_fedotov_thresh(4,6,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);
};
};

if ((Wgen>=1.2875)&&(Wgen<=1.3125)){
s12left_wright_bin = getsbin_fed(0, s12gen, S12_ARR_FED[9][0], S12_ARR_FED[0][0]);
 s12right_wright_bin = s12left_wright_bin +1;
 
 s23left_wright_bin = getsbin_fed(0, s23gen, S23_ARR_FED[9][0], S23_ARR_FED[0][0]);
 s23right_wright_bin = s23left_wright_bin +1;

for (Short_t i=0;i<6;i++){
interpol_fedotov(4,6,0,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wr_fed[i],i);

interpol_fedotov_thresh(4,6,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wl_fed[i],i);
};


};










//Here we Q2-scale sigma_t, sigma_t, sigma_c2f and sigma_cf with the fit_functions (corresponding pol1) for Wright_bin and Wleft_bin
if ((Q2gen>=0.575)&&(Q2gen<=0.65)){
//sigma_t
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_0575_065(Q2gen,Wright_bin,0)/pol1_fed_0575_065(0.575,Wright_bin,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_0575_065(Q2gen,Wleft_bin,0)/pol1_fed_0575_065(0.575,Wleft_bin,0);
////sigma_l
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_0575_065(Q2gen,Wright_bin,1)/pol1_fed_0575_065(0.575,Wright_bin,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_0575_065(Q2gen,Wleft_bin,1)/pol1_fed_0575_065(0.575,Wleft_bin,1);
//sigma_c2f
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_0575_065(Q2gen,Wright_bin,2)/pol1_fed_0575_065(0.575,Wright_bin,2);
sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_0575_065(Q2gen,Wleft_bin,2)/pol1_fed_0575_065(0.575,Wleft_bin,2);
//sigma_cf
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_0575_065(Q2gen,Wright_bin,4)/pol1_fed_0575_065(0.575,Wright_bin,4);
sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_0575_065(Q2gen,Wleft_bin,4)/pol1_fed_0575_065(0.575,Wleft_bin,4);


};

if ((Q2gen>=0.65)&&(Q2gen<=0.95)){
//sigma_t
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_0575_065(0.65,Wright_bin,0)/pol1_fed_0575_065(0.575,Wright_bin,0);
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_065_095(Q2gen,Wleft_bin,0)/pol1_fed_065_095(0.65,Wleft_bin,0);


sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_0575_065(0.65,Wleft_bin,0)/pol1_fed_0575_065(0.575,Wleft_bin,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_065_095(Q2gen,Wleft_bin,0)/pol1_fed_065_095(0.65,Wleft_bin,0);

////sigma_l
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_0575_065(0.65,Wright_bin,1)/pol1_fed_0575_065(0.575,Wright_bin,1);
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_065_095(Q2gen,Wright_bin,1)/pol1_fed_065_095(0.65,Wright_bin,1);

sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_0575_065(0.65,Wleft_bin,1)/pol1_fed_0575_065(0.575,Wleft_bin,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_065_095(Q2gen,Wleft_bin,1)/pol1_fed_065_095(0.65,Wleft_bin,1);

//sigma_c2f
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_0575_065(0.65,Wright_bin,2)/pol1_fed_0575_065(0.575,Wright_bin,2);
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_065_095(Q2gen,Wright_bin,2)/pol1_fed_065_095(0.65,Wright_bin,2);

sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_0575_065(0.65,Wleft_bin,2)/pol1_fed_0575_065(0.575,Wleft_bin,2);
sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_065_095(Q2gen,Wleft_bin,2)/pol1_fed_065_095(0.65,Wleft_bin,2);

//sigma_cf
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_0575_065(0.65,Wright_bin,4)/pol1_fed_0575_065(0.575,Wright_bin,4);
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_065_095(Q2gen,Wright_bin,4)/pol1_fed_065_095(0.65,Wright_bin,4);

sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_0575_065(0.65,Wleft_bin,4)/pol1_fed_0575_065(0.575,Wleft_bin,4);
sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_065_095(Q2gen,Wleft_bin,4)/pol1_fed_065_095(0.65,Wleft_bin,4);


};


if ((Q2gen>=0.95)&&(Q2gen<=1.3)){

//sigma_t
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_0575_065(0.65,Wright_bin,0)/pol1_fed_0575_065(0.575,Wright_bin,0);
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_065_095(0.95,Wleft_bin,0)/pol1_fed_065_095(0.65,Wleft_bin,0);
sigma_wr_fed[0] = sigma_wr_fed[0]*pol1_fed_095_130(Q2gen,Wright_bin,0)/pol1_fed_095_130(0.95,Wright_bin,0);


sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_0575_065(0.65,Wleft_bin,0)/pol1_fed_0575_065(0.575,Wleft_bin,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_065_095(0.95,Wleft_bin,0)/pol1_fed_065_095(0.65,Wleft_bin,0);
sigma_wl_fed[0] = sigma_wl_fed[0]*pol1_fed_095_130(Q2gen,Wleft_bin,0)/pol1_fed_095_130(0.95,Wleft_bin,0);

////sigma_l
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_0575_065(0.65,Wright_bin,1)/pol1_fed_0575_065(0.575,Wright_bin,1);
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_065_095(0.95,Wright_bin,1)/pol1_fed_065_095(0.65,Wright_bin,1);
sigma_wr_fed[1] = sigma_wr_fed[1]*pol1_fed_095_130(Q2gen,Wright_bin,1)/pol1_fed_095_130(0.95,Wright_bin,1);

sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_0575_065(0.65,Wleft_bin,1)/pol1_fed_0575_065(0.575,Wleft_bin,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_065_095(0.95,Wleft_bin,1)/pol1_fed_065_095(0.65,Wleft_bin,1);
sigma_wl_fed[1] = sigma_wl_fed[1]*pol1_fed_095_130(Q2gen,Wleft_bin,1)/pol1_fed_095_130(0.95,Wleft_bin,1);
//sigma_c2f
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_0575_065(0.65,Wright_bin,2)/pol1_fed_0575_065(0.575,Wright_bin,2);
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_065_095(0.95,Wright_bin,2)/pol1_fed_065_095(0.65,Wright_bin,2);
sigma_wr_fed[2] = sigma_wr_fed[2]*pol1_fed_095_130(Q2gen,Wright_bin,2)/pol1_fed_095_130(0.95,Wright_bin,2);

sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_0575_065(0.65,Wleft_bin,2)/pol1_fed_0575_065(0.575,Wleft_bin,2);
sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_065_095(0.95,Wleft_bin,2)/pol1_fed_065_095(0.65,Wleft_bin,2);
sigma_wl_fed[2] = sigma_wl_fed[2]*pol1_fed_095_130(Q2gen,Wleft_bin,2)/pol1_fed_095_130(0.95,Wleft_bin,2);
//sigma_cf
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_0575_065(0.65,Wright_bin,4)/pol1_fed_0575_065(0.575,Wright_bin,4);
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_065_095(0.95,Wright_bin,4)/pol1_fed_065_095(0.65,Wright_bin,4);
sigma_wr_fed[4] = sigma_wr_fed[4]*pol1_fed_095_130(Q2gen,Wright_bin,4)/pol1_fed_095_130(0.95,Wright_bin,4);

sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_0575_065(0.65,Wleft_bin,4)/pol1_fed_0575_065(0.575,Wleft_bin,4);
sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_065_095(0.95,Wleft_bin,4)/pol1_fed_065_095(0.65,Wleft_bin,4);
sigma_wl_fed[4] = sigma_wl_fed[4]*pol1_fed_095_130(Q2gen,Wleft_bin,4)/pol1_fed_095_130(0.95,Wleft_bin,4);



};


//sigma_s2f
sigma_wr_fed[3] = sigma_wr_fed[3]*Func_q2_dep(Q2gen)/Func_q2_dep(0.575); 
sigma_wl_fed[3] = sigma_wl_fed[3]*Func_q2_dep(Q2gen)/Func_q2_dep(0.575); 

//sigma_sf
sigma_wr_fed[5] = sigma_wr_fed[5]*Func_q2_dep(Q2gen)/Func_q2_dep(0.575); 
sigma_wl_fed[5] = sigma_wl_fed[5]*Func_q2_dep(Q2gen)/Func_q2_dep(0.575);




//We are doing 1dim linear W-interpolation
for (Short_t i=0;i<6;i++){
sigma_final[i] = 1./0.025;
//sigma_final[i] = sigma_final[i]*(sigma_wr_fed[i]*fabs(W_ARR_FED_THRESH[Wleft_bin]-Wgen)+sigma_wl_fed[i]*fabs(W_ARR_FED_THRESH[Wright_bin]-Wgen));

if ((Wgen>=1.2875)&&(Wgen<=1.3125)) sigma_final[i] = sigma_final[i]*(sigma_wr_fed[i]*fabs(1.2875-Wgen)+sigma_wl_fed[i]*fabs(1.3125-Wgen));

if ((Wgen>=1.2375)&&(Wgen<=1.2875)) sigma_final[i] = sigma_final[i]*(sigma_wr_fed[i]*fabs(W_ARR_FED_THRESH[Wleft_bin]-Wgen)+sigma_wl_fed[i]*fabs(W_ARR_FED_THRESH[Wright_bin]-Wgen));
};

//sigma_final[3] = 0.;
//sigma_final[5] = 0.;

};//end if ((Wgen>=1.3125)&&(Wgen<1.4125)&&(Q2gen>=0.575)&&(Q2gen<=1.3))











//We get explicitly different sigmas from the array


sigma_t_final = sigma_final[0];
sigma_l_final = sigma_final[1];
sigma_c2f_final = sigma_final[2];
sigma_s2f_final = sigma_final[3];
sigma_cf_final = sigma_final[4];
sigma_sf_final = sigma_final[5];


};
