#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "interpol_rip2.h"
#include "get_xsect_ripani.h"
#include "get_xsect_w16_18_lowq2_fit.h"
#include "get_xsect_golovach.h"
#include "interpol_int.h"
using namespace std;


Short_t getWbin_rip2 (Float_t W) {

Short_t bin;
if ((W>=1.837499)&&(W<=1.8625)) bin = 0;
if ((W>=1.8625)&&(W<=1.8875)) bin = 1;
if ((W>=1.8875)&&(W<=1.9125)) bin = 2;
if ((W>=1.9125)&&(W<=1.9375)) bin = 3;

if ((W>=1.9375)&&(W<=1.9625)) bin = 4;
if ((W>=1.9625)&&(W<=1.9875)) bin = 5;
if ((W>=1.9875)&&(W<=2.0125)) bin = 6;
if ((W>=2.0125)&&(W<=2.0375)) bin = 7;

if ((W>=2.0375)&&(W<=2.0625)) bin = 8;
if ((W>=2.0625)&&(W<=2.0875)) bin = 9;

if ((W>=2.0875)&&(W<=2.1125)) bin = 10;
if ((W>=2.1125)&&(W<=2.1375)) bin = 11;

if ((W>=2.1375)&&(W<=2.1875)) bin = 12;
if ((W>=2.1875)&&(W<=2.2375)) bin = 13;
if ((W>=2.2375)&&(W<=2.2875)) bin = 14;
if ((W>=2.2875)&&(W<=2.3375)) bin = 15;
if ((W>=2.3375)&&(W<=2.3875)) bin = 16;
if ((W>=2.3875)&&(W<=2.4375)) bin = 17;
if ((W>=2.4375)&&(W<=2.4875)) bin = 18;
if ((W>=2.4875)&&(W<=2.5375)) bin = 19;

if ((W < 1.8125)||(W > 2.538)) {
cout << "Error, wrong W range, Rip2 " << "\n";
bin = -100;
}

return bin;
};

//This subroutine estimates cross section in the region (W>=1.8375)&&(W<=2.5375)&&(Q2>around 0)&&(Q2<1.3)
//STEP 1. The subroutine takes cross sections for Q2 = 1.3 GeV^2 (data-based up to W = 2.1 GeV and pure model for W from 2.1 to 2.5 GeV) 
//This is the grid and xsect array for cross sections at Q2 = 1.3 GeV^2:
//W_ARR_RIP2[21];
//S12_ARR_RIP2[12][21];
//S23_ARR_RIP2[12][21];
//THETA_ARR[6]; 
//ALPHA_ARR[6];
//SIGMA_ARR_RIP2[6][21][12][12][6][6];

//STEP 2. The subroutine takes Golovach cross sections for W from 1.8375 to 2.5375 GeV
//STEP 3. The aforementioned two cross section sets are scaled and mixed together. Scaling and mixing is specific to different Q2 ranges and different sigmas.

void get_xsect_w18_25(Float_t E_beam,Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final){

Float_t A_tmp[11];
Float_t sigma_t_wright_gol,sigma_t_wleft_gol, sigma_t_rip2,sigma_l_rip2;
Float_t sigma_t_gol;
Short_t w_left_bin_gol;

//STEP 1. The subroutine takes cross sections for Q2 = 1.3 GeV^2 (data-based up to W = 2.1 GeV and pure model for W from 2.1 to 2.5 GeV) 

Short_t Wleft_bin = getWbin_rip2(Wgen);
Short_t Wright_bin = Wleft_bin+1;
//cout << Wgen <<" "<< Wleft_bin << " " <<Wright_bin <<"\n";

Short_t s12left_wleft_bin = getsbin(Wleft_bin, s12gen, S12_ARR_RIP2[11][Wleft_bin], S12_ARR_RIP2[0][Wleft_bin]);
Short_t s12right_wleft_bin = s12left_wleft_bin +1;

Short_t s12left_wright_bin = getsbin(Wright_bin, s12gen, S12_ARR_RIP2[11][Wright_bin], S12_ARR_RIP2[0][Wright_bin]);
Short_t s12right_wright_bin = s12left_wright_bin +1;

Short_t s23left_wleft_bin = getsbin(Wleft_bin, s23gen, S23_ARR_RIP2[11][Wleft_bin], S23_ARR_RIP2[0][Wleft_bin]);
Short_t s23right_wleft_bin = s23left_wleft_bin +1;

Short_t s23left_wright_bin = getsbin(Wright_bin, s23gen, S23_ARR_RIP2[11][Wright_bin], S23_ARR_RIP2[0][Wright_bin]);
Short_t s23right_wright_bin = s23left_wright_bin +1;

//cout << s23left_wleft_bin <<" "<< Wgen << " "<<Wleft_bin<< " " << s23gen << " "<< S23_ARR_RIP2[s23left_wleft_bin][Wleft_bin] << " "<< S23_ARR_RIP2[s23right_wleft_bin][Wleft_bin] << "\n";

Short_t thetaleft_bin = getanglebin(thetagen,THETA_ARR[5]);
Short_t thetaright_bin = thetaleft_bin+1;

Short_t alphaleft_bin = getanglebin(alphagen,ALPHA_ARR[5]);
Short_t alpharight_bin = alphaleft_bin+1;

//cout << alphaleft_bin << " "<< alphagen<< " "<< ALPHA_ARR[alphaleft_bin] << " "<< ALPHA_ARR[ alpharight_bin] << "\n";

Float_t sigma_final[6];
Float_t sigma_wright[6];
Float_t sigma_wleft[6];
Float_t sigma_wright_rip[6];
Float_t sigma_wleft_rip[6];
Float_t factor;


//then we are doing 4d-interpolation for each (Wleft_bin, Q2_left_bin),  (Wright_bin, Q2_left_bin), (Wright_bin, Q2_right_bin) and (Wleft_bin, Q2_right_bin) and obtain cross-secton in that points (4 GLOBAL 6dim arrays)
//0 - sigma_t, 1 - sigma_l, 2 - sigma_c2f, 3 - sigma_s2f, 4 - sigma_cf, 5 - sigma_sf
for (Short_t i=0;i<6;i++){
interpol_rip2(4,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright[i],i);

interpol_rip2(4,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft[i],i);
};
 
//STEP 2. The subroutine takes Golovach cross sections for W from 1.8375 to 2.5375 GeV
get_xsect_gol_datamod(Wgen, s12gen,s23gen, thetagen, alphagen, w_left_bin_gol, sigma_t_wright_gol,sigma_t_wleft_gol );

//if ((isnan(sigma_wright[0]))||(isnan(sigma_wleft[0]))) cout << sigma_wright[0]<<" "<<sigma_wleft[0]<<" "<<s12gen<<" "<<S12_ARR_RIP2[][Wleft_bin] <<s23gen<<" "<<s12left_wright_bin<<" "<< s23left_wright_bin<<" "<<s12left_wleft_bin<<" "<<s23left_wleft_bin<<" subr\n";

//STEP 3. The aforementioned two cross section sets are scaled and mixed together. Scaling and mixing is specific to different Q2 ranges and different sigmas.

//SIGMA_T
if ((Q2gen > 0.0002)&&(Q2gen <=0.65)){

//1dim W-interpolation for golovach
sigma_t_gol = 1./fabs(W_ARR_RIP2[Wright_bin]-W_ARR_RIP2[Wleft_bin]);
sigma_t_gol = sigma_t_gol*(sigma_t_wright_gol*fabs(W_ARR_RIP2[Wleft_bin]-Wgen)+sigma_t_wleft_gol*fabs(W_ARR_RIP2[Wright_bin]-Wgen));

sigma_t_gol = sigma_t_gol*func_sigma_t(Q2gen,16)/func_sigma_t(0.003,16);

//1dim W-interpolation for Q2 = 1.3 GeV^2
sigma_final[0] = 1./fabs(W_ARR_RIP2[Wright_bin]-W_ARR_RIP2[Wleft_bin]);
sigma_final[0] = sigma_final[0]*(sigma_wright[0]*fabs(W_ARR_RIP2[Wleft_bin]-Wgen)+sigma_wleft[0]*fabs(W_ARR_RIP2[Wright_bin]-Wgen));

//mixing
sigma_t_final =(Q2gen-0.0003)/1.3*3.*sigma_final[0] + (1.3-Q2gen)/1.3*sigma_t_gol;

//adjusting
if ((Q2gen > 0.4)&&(Q2gen <=0.65)&&(E_beam<2.45)) sigma_t_final =(-0.8*Q2gen+1.32)*sigma_t_final;
if ((Q2gen > 0.4)&&(Q2gen <=0.65)&&(E_beam>2.45)&&(E_beam<3.)) sigma_t_final =(-0.4*Q2gen+1.16)*sigma_t_final;

};


if ((Q2gen>=0.65)&&(Q2gen<=1.3)){

//Q2-scale for golovach
sigma_t_wright_gol = sigma_t_wright_gol*func_sigma_t(0.65,16)/func_sigma_t(0.003,16);
sigma_t_wleft_gol = sigma_t_wleft_gol*func_sigma_t(0.65,16)/func_sigma_t(0.003,16);

//1dim W-interpolation for golovach
sigma_t_gol = 1./fabs(W_ARR_RIP2[Wright_bin]-W_ARR_RIP2[Wleft_bin]);
sigma_t_gol = sigma_t_gol*(sigma_t_wright_gol*fabs(W_ARR_RIP2[Wleft_bin]-Wgen)+sigma_t_wleft_gol*fabs(W_ARR_RIP2[Wright_bin]-Wgen));

//1dim W-interpolation for Q2 = 1.3 GeV^2
sigma_final[0] = 1./fabs(W_ARR_RIP2[Wright_bin]-W_ARR_RIP2[Wleft_bin]);
sigma_final[0] = sigma_final[0]*(sigma_wright[0]*fabs(W_ARR_RIP2[Wleft_bin]-Wgen)+sigma_wleft[0]*fabs(W_ARR_RIP2[Wright_bin]-Wgen));

sigma_t_rip2 = sigma_final[0];

//mixing
sigma_t_gol = (1.3-0.65)/1.3*sigma_t_gol +(0.65-0.0003)/1.3*3.*sigma_final[0] ;


//1dim Q2-interpolation
sigma_t_final = 1./0.65;
sigma_t_final = sigma_t_final*(sigma_t_rip2*fabs(0.65 - Q2gen) + sigma_t_gol*fabs(1.3-Q2gen) );

//adjusting
//sigma_t_final = (4./13.*Q2gen+0.6)*sigma_t_final;
//sigma_t_final = (5./6.*Q2gen+5./24.)*sigma_t_final; 	
if ((E_beam>2.45)&&(E_beam<3)) sigma_t_final = (1./13.*Q2gen+0.85)*sigma_t_final;  
//if (E_beam<2.4) sigma_t_final = (4./13.*Q2gen+0.6)*sigma_t_final; 
if (E_beam<2.45) sigma_t_final = (2./13.*Q2gen+0.7)*sigma_t_final; 
//if (E_beam<2.5) sigma_t_final = 0.9*sigma_t_final; 
//sigma_t_final = 0.75*sigma_t_final;

};


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//SIGMA_L
if ((Q2gen>=0.95)&&(Q2gen<=1.3)){

sigma_wright[1] = sigma_wright[1]*(139.083 - Q2gen*96.2661)/(139.083 - 1.3*96.2661);
sigma_wleft[1] = sigma_wleft[1]*(139.083 - Q2gen*96.2661)/(139.083 - 1.3*96.2661);

sigma_l_final = 1./fabs(W_ARR_RIP2[Wright_bin]-W_ARR_RIP2[Wleft_bin]);
sigma_l_final = sigma_l_final*(sigma_wright[1]*fabs(W_ARR_RIP2[Wleft_bin]-Wgen)+sigma_wleft[1]*fabs(W_ARR_RIP2[Wright_bin]-Wgen));
};


if ((Q2gen>=0.65)&&(Q2gen<=0.95)){

sigma_wright[1] = sigma_wright[1]*(139.083 - 0.95*96.2661)/(139.083 - 1.3*96.2661);
sigma_wleft[1] = sigma_wleft[1]*(139.083 - 0.95*96.2661)/(139.083 - 1.3*96.2661);

sigma_wright[1] = sigma_wright[1]*(69.6255 - Q2gen*23.1535)/(69.6255 - 0.95*23.1535);
sigma_wleft[1] = sigma_wleft[1]*(69.6255 - Q2gen*23.1535)/(69.6255 - 0.95*23.1535);

sigma_l_final = 1./fabs(W_ARR_RIP2[Wright_bin]-W_ARR_RIP2[Wleft_bin]);
sigma_l_final = sigma_l_final*(sigma_wright[1]*fabs(W_ARR_RIP2[Wleft_bin]-Wgen)+sigma_wleft[1]*fabs(W_ARR_RIP2[Wright_bin]-Wgen));
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

sigma_l_final = 1./fabs(W_ARR_RIP2[Wright_bin]-W_ARR_RIP2[Wleft_bin]);
sigma_l_final = sigma_l_final*(sigma_wright[1]*fabs(W_ARR_RIP2[Wleft_bin]-Wgen)+sigma_wleft[1]*fabs(W_ARR_RIP2[Wright_bin]-Wgen));

};


sigma_l_final = 0.5*sigma_l_final;

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
sigma_final[i] = 1./fabs(W_ARR_RIP2[Wright_bin]-W_ARR_RIP2[Wleft_bin]);
sigma_final[i] = sigma_final[i]*(sigma_wright[i]*fabs(W_ARR_RIP2[Wleft_bin]-Wgen)+sigma_wleft[i]*fabs(W_ARR_RIP2[Wright_bin]-Wgen));
};

sigma_c2f_final =  sigma_final[2]; 
sigma_s2f_final = sigma_final[3];  
sigma_cf_final =  sigma_final[4]; 
sigma_sf_final =  sigma_final[5]; 


};
