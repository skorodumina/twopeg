#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"

using namespace std;

//This subroutine performs the cross section interpolation for Q2 = 0 GeV^2 and W from 2.5875 to 3.0375 GeV (pure-model based)
//This is the grid and xsect arrays for this cross sections:
//W_ARR_RIP3[10];
//S12_ARR_RIP3[16][10];
//S23_ARR_RIP3[16][10];
//THETA_ARR[6]; 
//ALPHA_ARR[6];
//SIGMA_ARR_GOL2[10][16][16][6][6]; (the same as sigma_t from SIGMA_ARR_RIP3[6][10][16][16][6][6], but with scaling factors, which are W, s12/s23, and angles dependent)


void interpol_gol2(Short_t Wbin, Short_t a_l_bin, Short_t a_r_bin, Short_t b_l_bin, Short_t b_r_bin, Short_t c_l_bin, Short_t c_r_bin,Short_t d_l_bin, Short_t d_r_bin,  Float_t a, Float_t b, Float_t c, Float_t d, Float_t &sigma_inter){


Short_t  s12_left_bin = a_l_bin;
Short_t  s12_right_bin =  a_r_bin;
Short_t  s23_left_bin = b_l_bin;
Short_t  s23_right_bin = b_r_bin;
Short_t  theta_left_bin = c_l_bin;
Short_t  theta_right_bin = c_r_bin;
Short_t  alpha_left_bin = d_l_bin;
Short_t  alpha_right_bin = d_r_bin;
 
Float_t  s12 = a;
Float_t  s23 = b;
Float_t  theta = c;
Float_t  alpha = d;


Float_t factor;
factor = 1./fabs(S12_ARR_RIP3[s12_right_bin][Wbin]-S12_ARR_RIP3[s12_left_bin][Wbin]);
factor = factor/fabs(S23_ARR_RIP3[s23_right_bin][Wbin]-S23_ARR_RIP3[s23_left_bin][Wbin]);
factor = factor/fabs(ALPHA_ARR[alpha_right_bin]-ALPHA_ARR[alpha_left_bin]);
factor = factor/fabs(THETA_ARR[theta_right_bin]-THETA_ARR[theta_left_bin]);

//1 -------------------
sigma_inter = SIGMA_ARR_GOL2[Wbin][s23_left_bin][s12_left_bin][theta_left_bin][alpha_left_bin];
sigma_inter = sigma_inter*fabs(S12_ARR_RIP3[s12_right_bin][Wbin]-s12);
sigma_inter = sigma_inter*fabs(S23_ARR_RIP3[s23_right_bin][Wbin]-s23);
sigma_inter = sigma_inter*fabs(THETA_ARR[theta_right_bin]-theta);
sigma_inter = sigma_inter*fabs(ALPHA_ARR[alpha_right_bin]-alpha);

//2 ---------------------------

Float_t sigma_tmp;

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_left_bin][s12_left_bin][theta_left_bin][alpha_right_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_right_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_right_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_right_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_left_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

//3 ---------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_left_bin][s12_left_bin][theta_right_bin][alpha_right_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_right_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_right_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_left_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_left_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

//4 ---------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_right_bin][s12_left_bin][theta_right_bin][alpha_right_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_right_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_left_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_left_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_left_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

//5 ---------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_right_bin][s12_right_bin][theta_right_bin][alpha_right_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_left_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_left_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_left_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_left_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 6 --------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_right_bin][s12_right_bin][theta_right_bin][alpha_left_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_left_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_left_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_left_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_right_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 7 -------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_right_bin][s12_right_bin][theta_left_bin][alpha_left_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_left_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_left_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_right_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_right_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 8 --------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_left_bin][s12_right_bin][theta_left_bin][alpha_left_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_left_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_right_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_right_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_right_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 9 --------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_left_bin][s12_left_bin][theta_right_bin][alpha_left_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_right_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_right_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_left_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_right_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 10 -------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_right_bin][s12_left_bin][theta_left_bin][alpha_left_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_right_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_left_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_right_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_right_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 11 --------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_right_bin][s12_right_bin][theta_left_bin][alpha_right_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_left_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_left_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_right_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_left_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 12 --------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_left_bin][s12_right_bin][theta_right_bin][alpha_right_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_left_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_right_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_left_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_left_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 13 --------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_right_bin][s12_left_bin][theta_right_bin][alpha_left_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_right_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_left_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_left_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_right_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 14 --------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_left_bin][s12_right_bin][theta_left_bin][alpha_right_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_left_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_right_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_right_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_left_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 15 --------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_right_bin][s12_left_bin][theta_left_bin][alpha_right_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_right_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_left_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_right_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_left_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

// 16 ---------------------------

sigma_tmp = SIGMA_ARR_GOL2[Wbin][s23_left_bin][s12_right_bin][theta_right_bin][alpha_left_bin];
sigma_tmp = sigma_tmp*fabs(S12_ARR_RIP3[s12_left_bin][Wbin]-s12);
sigma_tmp = sigma_tmp*fabs(S23_ARR_RIP3[s23_right_bin][Wbin]-s23);
sigma_tmp = sigma_tmp*fabs(THETA_ARR[theta_left_bin]-theta);
sigma_tmp = sigma_tmp*fabs(ALPHA_ARR[alpha_right_bin]-alpha);
sigma_inter = sigma_inter+sigma_tmp;

sigma_inter = sigma_inter*factor;

};
