#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "interpol_q2_13_wgt_3.h"
#include "interpol_phot_wgt_3.h"
#include "interpol_rip3.h"
#include "get_xsect_ripani.h"
#include "get_xsect_w16_18_lowq2_fit.h"
#include "get_xsect_golovach.h"
#include "interpol_int.h"
using namespace std;

Short_t getWbin_q2_13_wgt_3 (Float_t W) {

Short_t bin;
if ((W>=3.1375)&&(W<=3.2375)) bin = 0;
if ((W>=3.2375)&&(W<=3.3375)) bin = 1;
if ((W>=3.3375)&&(W<=3.4375)) bin = 2;
if ((W>=3.4375)&&(W<=3.5375)) bin = 3;
if ((W>=3.5375)&&(W<=3.6375)) bin = 4;
if ((W>=3.6375)&&(W<=3.7375)) bin = 5;
if ((W>=3.7375)&&(W<=3.8375)) bin = 6;
if ((W>=3.8375)&&(W<=3.9375)) bin = 7;
if ((W>=3.9375)&&(W<=4.0375)) bin = 8;
if ((W>=4.0375)&&(W<=4.1375)) bin = 9;
if ((W>=4.1375)&&(W<=4.2375)) bin = 10;
if ((W>=4.2375)&&(W<=4.3375)) bin = 11;
if ((W>=4.3375)&&(W<=4.4375)) bin = 12;
if ((W>=4.4375)&&(W<=4.5375)) bin = 13;

if ((W <3.1374 )||(W > 4.5376)) {
cout << "Error, wrong W range, Wgtr3" << "\n";
bin = -100;
};

return bin;
};

Short_t getsbin_wgt_3 (Float_t sgen, Float_t Smax, Float_t Smax2, Float_t Smin) {

Short_t bin;
if ((sgen>=Smin)&&(sgen<=Smax2)) bin = int((sgen-Smin)/((Smax2 - Smin)/11.));
if ((sgen>=Smax2)&&(sgen<=Smax)) bin = 11 + int((sgen-Smax2)/((Smax - Smax2)/4.));

if (sgen<Smin) bin = 0;
if (sgen>Smax) bin = 14;

return bin;
};


//This subroutine estimates cross section in the region (W>=3.1375)&&(W<=4.5375)&&(Q2>around 0)&&(Q2<1.3)
//STEP 1-2. The subroutine takes cross sections for Q2 = 1.3 GeV^2 and then for Q2 = 0 GeV^2 (both pure model-based) 
//This is the grid and xsect array for both cross section sets:
//W_ARR_gt_3[15];
//S12_ARR_gt_3[16][15];
//S23_ARR_gt_3[16][15];
//THETA_ARR[6]; 
//ALPHA_ARR[6];
//SIGMA_ARR_gt_3[6][15][16][16][6][6]; <---SIGMA_ARR_RIP3[6][10][16][16][6][6];
//SIGMA_ARR_phot_gt_3[15][16][16][6][6];,<---SIGMA_ARR_GOL2[10][16][16][6][6];<---SIGMA_ARR_RIP3[6][10][16][16][6][6];


//STEP 3. The aforementioned two cross section sets are scaled and mixed together. Scaling and mixing is specific to different Q2 ranges and different sigmas.

void get_xsect_w30_45(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final){

//STEP 1-2. The subroutine takes cross sections for Q2 = 1.3 GeV^2 and then for Q2 = 0 GeV^2 (both pure model-based) 
Short_t Wleft_bin = getWbin_q2_13_wgt_3(Wgen);
Short_t Wright_bin = Wleft_bin+1;

//cout <<Wgen <<" "<< Wleft_bin<<  "\n";

Short_t  s12left_wleft_bin = getsbin_wgt_3(s12gen, S12_ARR_gt_3[15][Wleft_bin],S12_ARR_gt_3[11][Wleft_bin], S12_ARR_gt_3[0][Wleft_bin]);
Short_t  s12right_wleft_bin = s12left_wleft_bin +1;

Short_t  s12left_wright_bin = getsbin_wgt_3(s12gen, S12_ARR_gt_3[15][Wright_bin],S12_ARR_gt_3[11][Wright_bin], S12_ARR_gt_3[0][Wright_bin]);
Short_t  s12right_wright_bin = s12left_wright_bin +1;
 
Short_t  s23left_wleft_bin = getsbin_wgt_3(s23gen, S23_ARR_gt_3[15][Wleft_bin],S23_ARR_gt_3[11][Wleft_bin], S23_ARR_gt_3[0][Wleft_bin]);
Short_t  s23right_wleft_bin = s23left_wleft_bin +1;

Short_t  s23left_wright_bin = getsbin_wgt_3(s23gen, S23_ARR_gt_3[15][Wright_bin],S23_ARR_gt_3[11][Wright_bin], S23_ARR_gt_3[0][Wright_bin]);
Short_t  s23right_wright_bin = s23left_wright_bin +1;

//cout << s23left_wleft_bin <<" "<< Wgen << " "<<Wleft_bin<< " " << s23gen << " "<< S23_ARR_gt_3[s23left_wleft_bin][Wleft_bin] << " "<< S23_ARR_gt_3[s23right_wleft_bin][Wleft_bin] << "\n";

Short_t thetaleft_bin = getanglebin(thetagen,THETA_ARR[5]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin(alphagen,ALPHA_ARR[5]);
Short_t alpharight_bin = alphaleft_bin+1;

//cout << thetaleft_bin << " "<< alphaleft_bin<< "\n";

Float_t sigma_final[6];
Float_t sigma_t_fin_ph;
Float_t sigma_wright[6];
Float_t sigma_wleft[6];
Float_t sigma_wright_ph;
Float_t sigma_wleft_ph,sigma_t_rip2;
Float_t factor;

//for Q2 = 1.3 GeV^2
for (Short_t j=0;j<6;j++){
interpol_q2_13_wgt_3(4,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright[j],j);

interpol_q2_13_wgt_3(4,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft[j],j);
};

//for Q2 = 0 GeV^2
interpol_phot_wgt_3(4,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright_ph);

interpol_phot_wgt_3(4,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft_ph);

//STEP 3. The aforementioned two cross section sets are scaled and mixed together. Scaling and mixing is specific to different Q2 ranges and different sigmas.

//SIGMA_T
if ((Q2gen > 0.0002)&&(Q2gen <=0.65)){

//1dim W-interpolation for photon point
sigma_t_fin_ph = 1./fabs(W_ARR_gt_3[Wright_bin]-W_ARR_gt_3[Wleft_bin]);
sigma_t_fin_ph = sigma_t_fin_ph*(sigma_wright_ph*fabs(W_ARR_gt_3[Wleft_bin]-Wgen)+sigma_wleft_ph*fabs(W_ARR_gt_3[Wright_bin]-Wgen));

//Q2-scaling
sigma_t_fin_ph = sigma_t_fin_ph*func_sigma_t(Q2gen,16)/func_sigma_t(0.003,16);

//1dim W-interpolation for Q2 = 1.3 GeV^2
sigma_final[0] = 1./fabs(W_ARR_gt_3[Wright_bin]-W_ARR_gt_3[Wleft_bin]);
sigma_final[0] = sigma_final[0]*(sigma_wright[0]*fabs(W_ARR_gt_3[Wleft_bin]-Wgen)+sigma_wleft[0]*fabs(W_ARR_gt_3[Wright_bin]-Wgen));

//mixing
sigma_t_final =(Q2gen-0.0003)/1.3*3.*sigma_final[0] + (1.3-Q2gen)/1.3*sigma_t_fin_ph;

};


if ((Q2gen>=0.65)&&(Q2gen<=1.3)){

//Q2-scale for photon point
sigma_wright_ph = sigma_wright_ph*func_sigma_t(0.65,16)/func_sigma_t(0.003,16);
sigma_wleft_ph = sigma_wleft_ph*func_sigma_t(0.65,16)/func_sigma_t(0.003,16);

//1dim W-interpolation for photon point
sigma_t_fin_ph = 1./fabs(W_ARR_gt_3[Wright_bin]-W_ARR_gt_3[Wleft_bin]);
sigma_t_fin_ph = sigma_t_fin_ph*(sigma_wright_ph*fabs(W_ARR_gt_3[Wleft_bin]-Wgen)+sigma_wleft_ph*fabs(W_ARR_gt_3[Wright_bin]-Wgen));

//1dim W-interpolation for Q2 = 1.3 GeV^2
sigma_final[0] = 1./fabs(W_ARR_gt_3[Wright_bin]-W_ARR_gt_3[Wleft_bin]);
sigma_final[0] = sigma_final[0]*(sigma_wright[0]*fabs(W_ARR_gt_3[Wleft_bin]-Wgen)+sigma_wleft[0]*fabs(W_ARR_gt_3[Wright_bin]-Wgen));

sigma_t_rip2 = sigma_final[0];

//mixing
sigma_t_fin_ph = (1.3-0.65)/1.3*sigma_t_fin_ph +(0.65-0.0003)/1.3*3.*sigma_final[0] ;


//1dim Q2-interpolation
sigma_t_final = 1./0.65;
sigma_t_final = sigma_t_final*(sigma_t_rip2*fabs(0.65 - Q2gen) + sigma_t_fin_ph*fabs(1.3-Q2gen) );

};

//SIGMA_L
if ((Q2gen>=0.95)&&(Q2gen<=1.3)){

sigma_wright[1] = sigma_wright[1]*(139.083 - Q2gen*96.2661)/(139.083 - 1.3*96.2661);
sigma_wleft[1] = sigma_wleft[1]*(139.083 - Q2gen*96.2661)/(139.083 - 1.3*96.2661);

sigma_l_final = 1./fabs(W_ARR_gt_3[Wright_bin]-W_ARR_gt_3[Wleft_bin]);
sigma_l_final = sigma_l_final*(sigma_wright[1]*fabs(W_ARR_gt_3[Wleft_bin]-Wgen)+sigma_wleft[1]*fabs(W_ARR_gt_3[Wright_bin]-Wgen));
};

if ((Q2gen>=0.65)&&(Q2gen<=0.95)){

sigma_wright[1] = sigma_wright[1]*(139.083 - 0.95*96.2661)/(139.083 - 1.3*96.2661);
sigma_wleft[1] = sigma_wleft[1]*(139.083 - 0.95*96.2661)/(139.083 - 1.3*96.2661);

sigma_wright[1] = sigma_wright[1]*(69.6255 - Q2gen*23.1535)/(69.6255 - 0.95*23.1535);
sigma_wleft[1] = sigma_wleft[1]*(69.6255 - Q2gen*23.1535)/(69.6255 - 0.95*23.1535);

sigma_l_final = 1./fabs(W_ARR_gt_3[Wright_bin]-W_ARR_gt_3[Wleft_bin]);
sigma_l_final = sigma_l_final*(sigma_wright[1]*fabs(W_ARR_gt_3[Wleft_bin]-Wgen)+sigma_wleft[1]*fabs(W_ARR_gt_3[Wright_bin]-Wgen));
//cout << sigma_l_final << "\n";
};

if ((Q2gen>=0.00002)&&(Q2gen<=0.65)){

sigma_wright[1] = sigma_wright[1]*(139.083 - 0.95*96.2661)/(139.083 - 1.3*96.2661);
sigma_wleft[1] = sigma_wleft[1]*(139.083 - 0.95*96.2661)/(139.083 - 1.3*96.2661);

sigma_wright[1] = sigma_wright[1]*(69.6255 - 0.65*23.1535)/(69.6255 - 0.95*23.1535);
sigma_wleft[1] = sigma_wleft[1]*(69.6255 - 0.65*23.1535)/(69.6255 - 0.95*23.1535);

//sigma_l
sigma_wright[1] = sigma_wright[1]*pol2(Q2gen,16,1)/pol2(0.65,16,1);
sigma_wright[1] = sigma_wright[1]*getEpsL(2.445,1.8125,0.65)/getEpsL(2.445,1.8125,Q2gen);

//sigma_l
sigma_wleft[1] = sigma_wleft[1]*pol2(Q2gen,16,1)/pol2(0.65,16,1);
sigma_wleft[1] = sigma_wleft[1]*getEpsL(2.445,1.8125,0.65)/getEpsL(2.445,1.8125,Q2gen);

sigma_l_final = 1./fabs(W_ARR_gt_3[Wright_bin]-W_ARR_gt_3[Wleft_bin]);
sigma_l_final = sigma_l_final*(sigma_wright[1]*fabs(W_ARR_gt_3[Wleft_bin]-Wgen)+sigma_wleft[1]*fabs(W_ARR_gt_3[Wright_bin]-Wgen));
};

//Other sigmas
//sigma_c2f
sigma_wright[2] = sigma_wright[2]*Func_q2_dep(Q2gen)/Func_q2_dep(1.3); 
sigma_wleft[2] = sigma_wleft[2]*Func_q2_dep(Q2gen)/Func_q2_dep(1.3); 

//sigma_s2f
sigma_wright[3] = sigma_wright[3]*Func_q2_dep(Q2gen)/Func_q2_dep(1.3); 
sigma_wleft[3] = sigma_wleft[3]*Func_q2_dep(Q2gen)/Func_q2_dep(1.3); 

//sigma_cf
sigma_wright[4] = sigma_wright[4]*Func_q2_dep(Q2gen)/Func_q2_dep(1.3); 
sigma_wleft[4] = sigma_wleft[4]*Func_q2_dep(Q2gen)/Func_q2_dep(1.3); 

//sigma_sf
sigma_wright[5] = sigma_wright[5]*Func_q2_dep(Q2gen)/Func_q2_dep(1.3); 
sigma_wleft[5] = sigma_wleft[5]*Func_q2_dep(Q2gen)/Func_q2_dep(1.3);


for (Short_t i=2;i<6;i++){
sigma_final[i] = 1./fabs(W_ARR_gt_3[Wright_bin]-W_ARR_gt_3[Wleft_bin]);
sigma_final[i] = sigma_final[i]*(sigma_wright[i]*fabs(W_ARR_gt_3[Wleft_bin]-Wgen)+sigma_wleft[i]*fabs(W_ARR_gt_3[Wright_bin]-Wgen));
};

sigma_c2f_final =  sigma_final[2]; 
sigma_s2f_final = sigma_final[3];  
sigma_cf_final =  sigma_final[4]; 
sigma_sf_final =  sigma_final[5];

};
