#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"



using namespace std;






Short_t getWbin_int (Float_t W) {
if ((W>1.2625)&&(W<3.0125)) return int((W-1.2625)/0.025);
if ((W<1.2625)||(W>3.0125)) {
cout << W<< "   WRONG W-VALUE. INTERPOLATION ERROR \n";
return -100.; 


};
};

Short_t getQ2bin_int (Float_t Q2) {

if ((Q2>0.0005)&&(Q2<0.0255)) return 0;

if ((Q2>0.0255)&&(Q2<1.2255)) return int((Q2-0.0255)/0.05)+1;

if ((Q2>1.2255)&&(Q2<1.3)) return 25;

if  ((Q2<0.0005)||(Q2>1.3)) {
cout << Q2<< " WRONG Q2-VALUE. INTERPOLATION ERROR \n";
return -100.; 

};

};



Float_t Func_q2_dep(Float_t Q2){

return 1./(1.+Q2/0.7)/(1.+Q2/0.7);
};





void interpol_int(Float_t Q2gen, Float_t Wgen,Float_t &sigma_t_int, Float_t &sigma_l_int){
Float_t Q2_old; 
Float_t W_old; 

Q2_old = Q2gen;
W_old = Wgen;
if (Q2gen>1.29) Q2gen = 1.29;
if (Q2gen<0.0006) Q2gen = 0.0006;
if (Wgen<1.2626) Wgen = 1.2626;
if (Wgen > 3.0125) Wgen = 3.0124;

Short_t Wleft_bin = getWbin_int(Wgen);
Short_t Wright_bin = Wleft_bin+1;



Short_t Q2left_bin = getQ2bin_int(Q2gen);
Short_t Q2right_bin = Q2left_bin+1;

//cout << Q2gen << " "<< Q2left_bin << "\n";
Float_t factor;
factor = 1./fabs(W_ARR_2pi_INT[Wright_bin]-W_ARR_2pi_INT[Wleft_bin]);
factor = factor/fabs(Q2_ARR_2pi_INT[Q2right_bin]-Q2_ARR_2pi_INT[Q2left_bin]);

sigma_t_int = SIGMA_T_ARR_2pi_INT[Q2left_bin][Wleft_bin]*fabs(W_ARR_2pi_INT[Wright_bin]-Wgen)*fabs(Q2_ARR_2pi_INT[Q2right_bin]-Q2gen);
sigma_t_int = sigma_t_int + SIGMA_T_ARR_2pi_INT[Q2left_bin][Wright_bin]*fabs(W_ARR_2pi_INT[Wleft_bin]-Wgen)*fabs(Q2_ARR_2pi_INT[Q2right_bin]-Q2gen);
sigma_t_int = sigma_t_int + SIGMA_T_ARR_2pi_INT[Q2right_bin][Wright_bin]*fabs(W_ARR_2pi_INT[Wleft_bin]-Wgen)*fabs(Q2_ARR_2pi_INT[Q2left_bin]-Q2gen);
sigma_t_int = sigma_t_int + SIGMA_T_ARR_2pi_INT[Q2right_bin][Wleft_bin]*fabs(W_ARR_2pi_INT[Wright_bin]-Wgen)*fabs(Q2_ARR_2pi_INT[Q2left_bin]-Q2gen);
sigma_t_int = sigma_t_int*factor;



sigma_l_int = SIGMA_L_ARR_2pi_INT[Q2left_bin][Wleft_bin]*fabs(W_ARR_2pi_INT[Wright_bin]-Wgen)*fabs(Q2_ARR_2pi_INT[Q2right_bin]-Q2gen);
sigma_l_int = sigma_l_int + SIGMA_L_ARR_2pi_INT[Q2left_bin][Wright_bin]*fabs(W_ARR_2pi_INT[Wleft_bin]-Wgen)*fabs(Q2_ARR_2pi_INT[Q2right_bin]-Q2gen);
sigma_l_int = sigma_l_int + SIGMA_L_ARR_2pi_INT[Q2right_bin][Wright_bin]*fabs(W_ARR_2pi_INT[Wleft_bin]-Wgen)*fabs(Q2_ARR_2pi_INT[Q2left_bin]-Q2gen);
sigma_l_int = sigma_l_int + SIGMA_L_ARR_2pi_INT[Q2right_bin][Wleft_bin]*fabs(W_ARR_2pi_INT[Wright_bin]-Wgen)*fabs(Q2_ARR_2pi_INT[Q2left_bin]-Q2gen);
sigma_l_int = sigma_l_int*factor;



//cout << SIGMA_T_ARR_2pi_INT[Q2left_bin][Wleft_bin] << " "<< SIGMA_T_ARR_2pi_INT[Q2left_bin][Wright_bin] << " "<< SIGMA_T_ARR_2pi_INT[Q2right_bin][Wright_bin] << " "<< SIGMA_T_ARR_2pi_INT[Q2right_bin][Wleft_bin] <<" " <<sigma_t_int << "\n";

if (Q2_old >1.29){
sigma_t_int = sigma_t_int*Func_q2_dep(Q2_old)/Func_q2_dep(1.29);
sigma_l_int = sigma_l_int*Func_q2_dep(Q2_old)/Func_q2_dep(1.29);

};
if (W_old<1.2375){
sigma_t_int = 0.;
sigma_l_int = 0.;
};


Q2gen=Q2_old;
Wgen=W_old;

};



