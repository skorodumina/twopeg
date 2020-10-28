#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "interpol_rip3.h"
#include "get_xsect_ripani.h"
#include "get_xsect_14_18_lowq2_fit.h"
#include "get_xsect_gol2.h"
#include "get_xsect_golovach.h"
#include "interpol_int.h"
using namespace std;



Short_t getWbin_rip3 (Float_t W) {

Short_t bin;
if ((W>=2.5875)&&(W<=2.6375)) bin = 0;
if ((W>=2.6375)&&(W<=2.6875)) bin = 1;
if ((W>=2.6875)&&(W<=2.7375)) bin = 2;
if ((W>=2.7375)&&(W<=2.7875)) bin = 3;
if ((W>=2.7875)&&(W<=2.8375)) bin = 4;
if ((W>=2.8375)&&(W<=2.8875)) bin = 5;
if ((W>=2.8875)&&(W<=2.9375)) bin = 6;
if ((W>=2.9375)&&(W<=2.9875)) bin = 7;
if ((W>=2.9875)&&(W<=3.0375)) bin = 8;


if ((W <2.58749 )||(W > 3.038)) {
cout << "Error, wrong W range, WRIP3 " << "\n";
bin = -100;
}
return bin;
};
//This subrouting is doing the following:
//1 - using auxiliary functions (getWbin, getQ2bin, getsbin, getanglebin) we identify the number of left and right point 
//2 - then we are doing 4d-interpolation for each (Wleft_bin, Q2_left_bin),  (Wright_bin, Q2_left_bin), (Wright_bin, Q2_right_bin) and (Wleft_bin, Q2_right_bin) and obtain cross-secton in that points (4 GLOBAL 6dim arrays)
//3 - then we are doing 2d-interpolation for the points written above and obtain sigma_final[6]


void get_xsect_25_30(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final){

Short_t Wleft_bin = getWbin_rip3(Wgen);
Short_t Wright_bin = Wleft_bin+1;
Float_t A_tmp[11];
//cout << getWbin_rip2(2.0875) << "\n";

//cout << Wgen <<" "<< Wleft_bin << " " <<Wright_bin <<"\n";

Float_t sigma_t_wright_gol,sigma_t_wleft_gol, sigma_t_rip2,sigma_l_rip2;
Float_t sigma_t_gol;
Short_t w_left_bin_gol;

//cout << Wgen << " "<< Wleft_bin << " t \n"; 



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
interpol_rip3(4,Wright_bin,s12left_wright_bin,s12right_wright_bin,s23left_wright_bin,s23right_wright_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wright[i],i);

interpol_rip3(4,Wleft_bin,s12left_wleft_bin,s12right_wleft_bin,s23left_wleft_bin,s23right_wleft_bin,thetaleft_bin,thetaright_bin,alphaleft_bin,alpharight_bin,s12gen,s23gen,thetagen,alphagen,sigma_wleft[i],i);

};

get_xsect_gol2(Wgen, s12gen,s23gen, thetagen, alphagen, w_left_bin_gol, sigma_t_wright_gol,sigma_t_wleft_gol );

//for sigma transferse

if ((Q2gen > 0.0002)&&(Q2gen <=0.65)){


//1dim W-interpolation for golovach
sigma_t_gol = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_t_gol = sigma_t_gol*(sigma_t_wright_gol*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_t_wleft_gol*fabs(W_ARR_RIP3[Wright_bin]-Wgen));

sigma_t_gol = sigma_t_gol*func_sigma_t(Q2gen,16)/func_sigma_t(0.003,16);



//for (Short_t i=0;i<6;i++){
sigma_final[0] = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_final[0] = sigma_final[0]*(sigma_wright[0]*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_wleft[0]*fabs(W_ARR_RIP3[Wright_bin]-Wgen));

//};
//sigma_t_rip2 = sigma_final[0];

sigma_t_final =(Q2gen-0.0003)/1.3*3.*sigma_final[0] + (1.3-Q2gen)/1.3*sigma_t_gol;

//if ((Q2gen > 0.4)&&(Q2gen <=0.65)&&(E_beam<2.45)) sigma_t_final =(-0.8*Q2gen+1.32)*sigma_t_final;
//if ((Q2gen > 0.4)&&(Q2gen <=0.65)&&(E_beam>2.45)&&(E_beam<3.)) sigma_t_final =(-0.4*Q2gen+1.16)*sigma_t_final;

//sigma_t_final = sigma_t_gol;
};


if ((Q2gen>=0.65)&&(Q2gen<=1.3)){

//Q2-scale for golovach
sigma_t_wright_gol = sigma_t_wright_gol*func_sigma_t(0.65,16)/func_sigma_t(0.003,16);
sigma_t_wleft_gol = sigma_t_wleft_gol*func_sigma_t(0.65,16)/func_sigma_t(0.003,16);

//1dim W-interpolation for golovach
sigma_t_gol = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_t_gol = sigma_t_gol*(sigma_t_wright_gol*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_t_wleft_gol*fabs(W_ARR_RIP3[Wright_bin]-Wgen));


sigma_final[0] = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_final[0] = sigma_final[0]*(sigma_wright[0]*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_wleft[0]*fabs(W_ARR_RIP3[Wright_bin]-Wgen));


sigma_t_rip2 = sigma_final[0];


sigma_t_gol = (1.3-0.65)/1.3*sigma_t_gol +(0.65-0.0003)/1.3*3.*sigma_final[0] ;


//1dim Q2-interpolation
sigma_t_final = 1./0.65;
sigma_t_final = sigma_t_final*(sigma_t_rip2*fabs(0.65 - Q2gen) + sigma_t_gol*fabs(1.3-Q2gen) );

//if ((E_beam>2.45)&&(E_beam<3)) sigma_t_final = (1./13.*Q2gen+0.85)*sigma_t_final;  
//if (E_beam<2.45) sigma_t_final = (2./13.*Q2gen+0.7)*sigma_t_final; 
};


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//for sigma longitudinal


if ((Q2gen>=0.95)&&(Q2gen<=1.3)){

sigma_wright[1] = sigma_wright[1]*(139.083 - Q2gen*96.2661)/(139.083 - 1.3*96.2661);
sigma_wleft[1] = sigma_wleft[1]*(139.083 - Q2gen*96.2661)/(139.083 - 1.3*96.2661);


sigma_l_final = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_l_final = sigma_l_final*(sigma_wright[1]*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_wleft[1]*fabs(W_ARR_RIP3[Wright_bin]-Wgen));

};


if ((Q2gen>=0.65)&&(Q2gen<=0.95)){

sigma_wright[1] = sigma_wright[1]*(139.083 - 0.95*96.2661)/(139.083 - 1.3*96.2661);
sigma_wleft[1] = sigma_wleft[1]*(139.083 - 0.95*96.2661)/(139.083 - 1.3*96.2661);

sigma_wright[1] = sigma_wright[1]*(69.6255 - Q2gen*23.1535)/(69.6255 - 0.95*23.1535);
sigma_wleft[1] = sigma_wleft[1]*(69.6255 - Q2gen*23.1535)/(69.6255 - 0.95*23.1535);


sigma_l_final = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_l_final = sigma_l_final*(sigma_wright[1]*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_wleft[1]*fabs(W_ARR_RIP3[Wright_bin]-Wgen));
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


sigma_l_final = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_l_final = sigma_l_final*(sigma_wright[1]*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_wleft[1]*fabs(W_ARR_RIP3[Wright_bin]-Wgen));

};

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
sigma_final[i] = 1./fabs(W_ARR_RIP3[Wright_bin]-W_ARR_RIP3[Wleft_bin]);
sigma_final[i] = sigma_final[i]*(sigma_wright[i]*fabs(W_ARR_RIP3[Wleft_bin]-Wgen)+sigma_wleft[i]*fabs(W_ARR_RIP3[Wright_bin]-Wgen));

};

sigma_c2f_final =  sigma_final[2]; 
sigma_s2f_final = sigma_final[3];  
sigma_cf_final =  sigma_final[4]; 
sigma_sf_final =  sigma_final[5]; 

};
