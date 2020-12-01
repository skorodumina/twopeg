#include <iomanip>
#include <string>
#include <stdio.h> 
#include <stdlib.h>  
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
#include "global.h"
#include <sstream>

using namespace std;

//This subroutine reads cross sections from files and fills out the following GLOBAL arrays:
//FOR RIPANI CROSS SECTIONS:
//W[17], Q2[3]
//S12_ARR[12][17],  S23_ARR[12][17],  THETA_ARR[6], ALPHA_ARR[6]
//SIGMA_ARR[6][3][17][12][12][6][6] <-> SIGMA_ARR[flag_sigma][q2_bin][w_binj][s23_bin][s12_bin][theta_bin][alpha_bin],
//where the flag_sigma corresponds to 0 - sigma_t, 1 - sigma_l, 2 - sigma_c2f, 3 - sigma_s2f, 4 - sigma_cf, 5 - sigma_sf
//----------------
//FOR GOLOVACH PHOTOPRODUCTION CROSS SECTIONS 
//W_ARR_GOL[30];
//S12_ARR_GOL[16][30];
//S23_ARR_GOL[16][30];
//THETA_ARR_GOL[14]; 
//ALPHA_ARR_GOL[14];
//SIGMA_ARR_GOL[30][16][16][14][14]; - correspond for sigma_t

void  read_xsect_files(){

string file_names[51];
string file_names_gol[30];
string file_names_fed[56];

string file_names_rip2[21];
string file_names_rip3[10];

string file_names_int[27];


Short_t wbin, q2bin;


//files' names (Ripani_cr_sects)
//Ripani cross sections (1.4125 < W < 1.8125 GeV, 0.65 < Q2 < 1.3 GeV2) \n";
for (Short_t ww=0; ww<17;ww++) {
PATH << data_dir_2pi.str() << "data/ripani_cr_sect/rip_4diffsec_065_" <<14125 + 250*ww << ".dat";
file_names[ww] = PATH.str();
PATH.str("");
};
for (Short_t ww=0; ww<17;ww++) {
PATH << data_dir_2pi.str() << "data/ripani_cr_sect/rip_4diffsec_095_" <<14125 + 250*ww << ".dat";
file_names[ww+17] = PATH.str();
PATH.str("");
};
for (Short_t ww=0; ww<17;ww++) {
PATH << data_dir_2pi.str() << "data/ripani_cr_sect/rip_4diffsec_130_" <<14125 + 250*ww << ".dat";
file_names[ww+34] = PATH.str();
PATH.str("");
};



//files' names (Golovach_cr_sects)
// Golovach cross sections (1.6125 < W < 2.5375 GeV at photon point) \n";

PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec161.dat";
file_names_gol[0] =  PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec164.dat";
file_names_gol[1] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec166.dat";
file_names_gol[2] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec169.dat";
file_names_gol[3] =PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec171.dat";
file_names_gol[4] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec174.dat";
file_names_gol[5] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec176.dat";
file_names_gol[6] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec179.dat";
file_names_gol[7] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec181.dat";
file_names_gol[8] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec184.dat";
file_names_gol[9] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec186.dat";
file_names_gol[10] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec189.dat";
file_names_gol[11] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec191.dat";
file_names_gol[12] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec194.dat";
file_names_gol[13] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec196.dat";
file_names_gol[14] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec199.dat";
file_names_gol[15] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec201.dat";
file_names_gol[16] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec204.dat";
file_names_gol[17] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec206.dat";
file_names_gol[18] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec209.dat";
file_names_gol[19] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec211.dat";
file_names_gol[20] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec214.dat";
file_names_gol[21] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec219.dat";
file_names_gol[22] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec224.dat";
file_names_gol[23] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec228.dat";
file_names_gol[24] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec234.dat";
file_names_gol[25] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec239.dat";
file_names_gol[26] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec244.dat";
file_names_gol[27] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec249.dat";
file_names_gol[28] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/golovach_cr_sect/5diffsec254.dat";
file_names_gol[29] = PATH.str();

//files' names (Fedotov_cr_sects)
//Fedotov cross sections (1.3125 < W < 1.5875 GeV, 0.225 < Q2 < 0.575 GeV2) \n";

for (Short_t ww=0; ww<4;ww++) {
PATH.str("");
PATH << data_dir_2pi.str() <<  "data/fedotov_cr_sect/fedotov_4diffsec_0225_" << 15125 + 250*ww << ".dat";
file_names_fed[ww] = PATH.str();
PATH.str("");
};
for (Short_t ww=0; ww<11;ww++) {
PATH.str("");
PATH << data_dir_2pi.str() <<  "data/fedotov_cr_sect/fedotov_4diffsec_0275_" << 13125 + 250*ww << ".dat";
file_names_fed[ww+4] =  PATH.str();
PATH.str("");
};
for (Short_t ww=0; ww<10;ww++) {
PATH.str("");
PATH << data_dir_2pi.str() <<  "data/fedotov_cr_sect/fedotov_4diffsec_0325_" << 13125 + 250*ww << ".dat";
file_names_fed[ww+15] =  PATH.str();
PATH.str("");
};
for (Short_t ww=0; ww<9;ww++) {
PATH.str("");
PATH << data_dir_2pi.str() <<  "data/fedotov_cr_sect/fedotov_4diffsec_0425_" << 13125 + 250*ww << ".dat";
file_names_fed[ww+25] = PATH.str(); 
PATH.str("");
};
for (Short_t ww=0; ww<8;ww++) {
PATH.str("");
PATH << data_dir_2pi.str() <<  "data/fedotov_cr_sect/fedotov_4diffsec_0475_" << 13125 + 250*ww << ".dat";
file_names_fed[ww+34] = PATH.str();
PATH.str("");
};
for (Short_t ww=0; ww<8;ww++) {
PATH.str("");
PATH << data_dir_2pi.str() <<  "data/fedotov_cr_sect/fedotov_4diffsec_0525_" << 13125 + 250*ww << ".dat";
file_names_fed[ww+42] =  PATH.str();
PATH.str("");
};
for (Short_t ww=0; ww<6;ww++) {
PATH.str("");
PATH << data_dir_2pi.str() << "data/fedotov_cr_sect/fedotov_4diffsec_0575_" << 13125 + 250*ww << ".dat"; 
file_names_fed[ww+50] =  PATH.str();
PATH.str("");
};


//model cross sections (1.8375 < W < 2 GeV, Q2 = 1.3 GeV2) \n";
PATH.str("");
for (Short_t ww=0; ww<13; ww++) {
PATH << data_dir_2pi.str() << "data/rip_q2_130_w_18_21/wgt18_4diffsec_130_"<< 18375 + 250*ww << ".dat";
file_names_rip2[ww] = PATH.str();
PATH.str("");
};
for (Short_t ww=13; ww<21; ww++) {
PATH << data_dir_2pi.str() << "data/rip_q2_130_w_18_21/wgt18_4diffsec_130_"<<  21875+ 500*(ww-13) << ".dat";
file_names_rip2[ww] = PATH.str(); 
PATH.str("");
};

for (Short_t ww=0; ww<10; ww++) {
PATH.str("");
PATH << data_dir_2pi.str() <<   "data/rip_q2_130_w_18_21/wgt18_4diffsec_130_" << 25875 + ww*500 << ".dat";
file_names_rip3[ww] = PATH.str(); 
PATH.str("");
};

//Integral cr sect
for (Short_t ww=0; ww<27; ww++) {
PATH.str("");
PATH << data_dir_2pi.str() <<  "data/int_sec_new/intsec_q2_" << ww << ".dat";
file_names_int[ww] = PATH.str();
PATH.str("");
};




cout<< "Reading integral cross sections for rad eff \n";
for (Short_t i=0; i<27; i++) {


 Q2_ARR_2pi_INT[0] = 0.0005;
 Q2_ARR_2pi_INT[26] = 1.3;
if ((i>0)&&(i<26)) Q2_ARR_2pi_INT[i] = 0.0255 +0.05*(i-1);


string dummy,xsect;

string file=file_names_int[i];
ifstream input(file.c_str());
if (input.fail()) cout << "ALARM! File # "<< i<< " with integral cross sections FAILED to open! \n"; 
if(input.is_open()){
for (Int_t iwint = 1; iwint <=71; iwint++) {
getline(input,xsect);
W_ARR_2pi_INT[iwint-1] = atof(xsect.c_str());


getline(input,xsect);
SIGMA_T_ARR_2pi_INT[i][iwint-1]= atof(xsect.c_str());
getline(input,xsect);
SIGMA_L_ARR_2pi_INT[i][iwint-1]= atof(xsect.c_str());

};
};
input.close();

};








//Define theta, alpha and Q2 arrays for Ripani

THETA_ARR[0] = 0.;
THETA_ARR[1] = 0.6343185;
THETA_ARR[2] = 1.258637;
THETA_ARR[3] = 1.882956;
THETA_ARR[4] = 2.507274;
THETA_ARR[5] = 3.141593;


ALPHA_ARR[0] = 0.;
ALPHA_ARR[1] = 1.262637;
ALPHA_ARR[2] = 2.515274;
ALPHA_ARR[3] = 3.767911;
ALPHA_ARR[4] = 5.020548;
ALPHA_ARR[5] = 6.283186;

Q2_ARR[0] = 0.65;
Q2_ARR[1] = 0.95;
Q2_ARR[2] = 1.30;



cout<<"Reading Ripani   cross sections \n";
//loop over files for Ripani
for (Short_t i=0; i<=50; i++) {

wbin = i - Int_t(i/17)*17;
q2bin = Int_t(i/17);
//Define W array
W_ARR[i - Int_t(i/17)*17] = 1.4125 + 0.025*(i - Int_t(i/17)*17);



string dummy,xsect;

string file=file_names[i];
ifstream input(file.c_str());
if (input.fail()) cout << "ALARM! File for q2bin = "<< q2bin<<", W = "<<  W_ARR[i - Int_t(i/17)*17]<<" GeV with Ripani cr sect FAILED to open! \n";
if(input.is_open()){

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
getline(input,xsect);
//Define 2dim s12 and s23 arrays
S12_ARR[is12-1][wbin] = atof(xsect.c_str());
getline(input,xsect);
S23_ARR[is23-1][wbin] = atof(xsect.c_str());
getline(input,dummy);
getline(input,dummy);

getline(input,xsect);
//sigma_t
SIGMA_ARR[0][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_l
SIGMA_ARR[1][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_c2f
SIGMA_ARR[2][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_s2f
SIGMA_ARR[3][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_cf
SIGMA_ARR[4][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_sf
SIGMA_ARR[5][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
EPSILON_L_RIPANI[q2bin][wbin] = atof(xsect.c_str());
getline(input,dummy);



};
};
};
};
};
input.close();

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {



for (Short_t j=0;j<6;j++){


SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR[11][wbin]-S12_ARR[0][wbin])*(S23_ARR[11][wbin]-S23_ARR[0][wbin]);

SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*M_PI*2.*M_PI;


if ((q2bin==2)&&(wbin==15)) SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.02;
if ((q2bin==2)&&(wbin==16)) SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.95;
};

};
};
};
};

};//end loop-ripani  - i
///---------------------------------------------------------------------------------------
//Define theta, alpha arrays for Golovach

THETA_ARR_GOL[0] = 0.;
THETA_ARR_GOL[1] = 0.2501225;
THETA_ARR_GOL[2] =0.4902450;
THETA_ARR_GOL[3] =0.7303675;
THETA_ARR_GOL[4] =0.9704900;
THETA_ARR_GOL[5] =1.210613;
THETA_ARR_GOL[6] =1.450735;
THETA_ARR_GOL[7] =1.690858;
THETA_ARR_GOL[8] =1.930980;
THETA_ARR_GOL[9] =2.171103;
THETA_ARR_GOL[10] =2.411225;
THETA_ARR_GOL[11] =2.651348;
THETA_ARR_GOL[12] =2.891470;
THETA_ARR_GOL[13] =3.141593;

ALPHA_ARR_GOL[0] = 0.;
ALPHA_ARR_GOL[1] = 0.4917835;
ALPHA_ARR_GOL[2] = 0.9735669;
ALPHA_ARR_GOL[3] = 1.455350;
ALPHA_ARR_GOL[4] = 1.937134;
ALPHA_ARR_GOL[5] = 2.418917;
ALPHA_ARR_GOL[6] = 2.900701;
ALPHA_ARR_GOL[7] = 3.382484;
ALPHA_ARR_GOL[8] = 3.864268;
ALPHA_ARR_GOL[9] = 4.346051;
ALPHA_ARR_GOL[10] = 4.827835;
ALPHA_ARR_GOL[11] = 5.309618;
ALPHA_ARR_GOL[12] = 5.791402;
ALPHA_ARR_GOL[13] = 6.283185;

cout<<"Reading Golovach cross sections \n";
//loop over files for Golovach
for (Short_t i=0; i<=29; i++) {


//Define W array
if ((i>=0)&&(i<=21)) W_ARR_GOL[i] = 1.6125 + 0.025*i;

if ((i>=22)&&(i<=29)) W_ARR_GOL[i] = 2.1875 + 0.05*(i-22);



string dummy,xsect;

string file=file_names_gol[i];
ifstream input(file.c_str());
if (input.fail()) cout << "ALARM! File for W = "<< W_ARR_GOL[i] <<"GeV with Golovach photoproduction cr sect FAILED to open! \n";
if(input.is_open()){

for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=14; itheta++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {
getline(input,xsect);
//Define 2dim s12 and s23 arrays
S12_ARR_GOL[is12-1][i] = atof(xsect.c_str());
getline(input,xsect);
S23_ARR_GOL[is23-1][i] = atof(xsect.c_str());
getline(input,dummy);
getline(input,dummy);
getline(input,xsect);
//sigma_t
SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
getline(input,dummy);



};
};
};
};
};
input.close();





for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=14; itheta++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {



SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_GOL[15][i]-S12_ARR_GOL[0][i])*(S23_ARR_GOL[15][i]-S23_ARR_GOL[0][i]);

SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*M_PI*2.*M_PI;
//SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*2.*M_PI;

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.932;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.928;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.98;

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.92*1.006*1.000051212;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.98*1.00773713;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.91*0.993;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.97*1.0129949;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.986;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.03;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.03*1.019;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.1*1.05*0.9696;



//----------------------------------------
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07*0.98;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07*0.95;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.05*0.95;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.04*0.95;

if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.06*0.95;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.04*0.95;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07*0.95;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.13*0.95;


if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.01*0.95;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.06*0.95;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*0.99*0.95;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.0105*0.95;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.02*0.95;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.08*0.95;
//---------------------

if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.03*0.95;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.03*1.01;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.16*1.02;

if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.02*0.95;

if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.01*1.04*0.95;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.1*1.0*0.954;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07*1.07*0.95;

if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][ialpha-1]*1.07;
};
};
};
};

for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.8;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.91;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.94;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.1;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.3;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.3;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.15;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.15;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.6;


if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.91;
//if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.1;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.3;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.35;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.15;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.15;
if (i==1) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.2;


if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.78;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.1;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.3;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.3;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.2;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.2;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.2;


if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.62;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.97;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1]*1.05;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.2;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.35;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.3;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.;


if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*1.25;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.9;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*1.25;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*1.35;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.55;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.5;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.5;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.5;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*1.3;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*1.1;
//if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.9;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.5;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.4;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*1.6;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*2.7;

if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.85;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.9;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.8;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.05;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.1;

if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.87;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.83;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.9;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.85;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.85;
if (i==7) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*0.9;
//-----------------------------------------------
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.75;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.81;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.9;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.7;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.03;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.8;



if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.52;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.58;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.64;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1]*0.88;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.95;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.8;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.03;

if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.85;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.85;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*0.9;


if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.515;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.54;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.69;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1]*0.92;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.8;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.07;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.9;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.9;




if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.43;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.6;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.8;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1]*0.97;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.76;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.71;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.07;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.84;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.78;




if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.43;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.57;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.83;

if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][7][ialpha-1]*0.97;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.73;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.05;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.96;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.96;






if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.43;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.6;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.8;

if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][8][ialpha-1]*0.9;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.75;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*1.07;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.94;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.96;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.68;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.58;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.74;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.85;








if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.48;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.4;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.65;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.7;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*0.9;




if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.25;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.44;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.55;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][3][ialpha-1]*0.85;
//if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][9][ialpha-1]*0.78;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*1.2;
//if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.96;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.96;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.2;



if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.83;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.53;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.87;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.6;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.1;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.47;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.83;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.7;




if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.233;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.48;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;


if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.21;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.48;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.78;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;


if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.2;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.3;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.65;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.3;


//-------------------------------------------------


if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.28;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.55;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;

if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.75;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.2;



if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.4;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.58;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;

if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][10][ialpha-1]*0.8;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.8;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.5;


if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.25;



if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.5;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.75;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.8;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.6;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][13][ialpha-1]*1.2;


if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.76;
//if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.97;
//if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.8;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.4;




if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*0.95;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.9;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.95;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.4;


if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*1.07;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.95;
//if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.9;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.3;


if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.81;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][2][ialpha-1]*0.85;

if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][11][ialpha-1]*0.6;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.3;



if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][0][ialpha-1]*1.23;
//if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][1][ialpha-1]*0.93;


if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][is12-1][12][ialpha-1]*0.27;

};
};
};


for (Int_t is23 = 1; is23 <=16; is23++) {

for (Int_t itheta = 1; itheta <=14; itheta++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {

SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.94;


if (i==0) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.9;
if (i==0) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.95;
if (i==0) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.95;
if (i==0) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.75;

if (i==1) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.9;
if (i==1) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.92;
if (i==1) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.92;
if (i==1) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.9;


if (i==2) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.8;
if (i==2) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;

if (i==3) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.75;
if (i==3) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;
if (i==3) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.1;
if (i==3) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.15;
if (i==3) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.15;
if (i==3) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==3) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.05;



if (i==4) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.75;
if (i==4) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;
if (i==4) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.8;
if (i==4) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.05;
if (i==4) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.1;


if (i==5) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.75;
if (i==5) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;
if (i==5) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.8;
if (i==5) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.94;
if (i==5) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.98;
if (i==5) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.85;

if (i==6) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.55;
if (i==6) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.85;
if (i==6) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.8;
if (i==6) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.1;
if (i==6) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.15;


if (i==7) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.55;
if (i==7) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.85;
if (i==7) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.78;
if (i==7) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*0.95;
if (i==7) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.8;


//if ((i>=8)&&(i<=12)) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.95;

//--------------------------------------------------------
if ((i>=8)&&(i<=12)) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.95;

if (i==8) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.25;
if (i==8) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==8) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.1;
if (i==8) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.3;
if (i==8) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.3;
if (i==8) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.34;

if (i==8) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.34;
if (i==8) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.26;
if (i==8) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*1.1;
if (i==8) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==8) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.5;

if (i==9) SIGMA_ARR_GOL[i][is23-1][0][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][0][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.5;
if (i==9) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.56;
if (i==9) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;

if (i==9) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.13;
if (i==9) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==9) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.55;
if (i==9) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.55;
if (i==9) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.54;

if (i==9) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.4;
if (i==9) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.12;
if (i==9) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.9;
if (i==9) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.9;
if (i==9) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.52;
if (i==9) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.37;


if (i==10) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==10) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.55;
if (i==10) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==10) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.15;
if (i==10) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.27;
if (i==10) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.27;
if (i==10) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.29;
if (i==10) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;

if (i==10) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.95;
if (i==10) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.95;
if (i==10) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.95;
if (i==10) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==10) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==11) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==11) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.68;
if (i==11) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==11) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.14;
if (i==11) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.32;
if (i==11) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.37;
if (i==11) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.36;
if (i==11) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.17;



if (i==11) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.9;
if (i==11) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==11) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==12) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==12) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==12) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==12) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==12) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.39;
if (i==12) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.42;
if (i==12) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==12) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.17;




if (i==12) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==12) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


//if (i==13) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.93;

if (i==13) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==13) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==13) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==13) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==13) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.43;
if (i==13) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.47;
if (i==13) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.42;
if (i==13) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.13;

if (i==13) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==13) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.6;
if (i==13) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.5;
if (i==13) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.45;
if (i==13) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.27;


if (i==14) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.7;
if (i==14) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][0][ialpha-1]*0.93;
if (i==14) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][0][ialpha-1]*0.97;


if (i==14) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==14) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==14) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==14) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==14) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.39;
if (i==14) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.42;
if (i==14) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==14) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.14;

if (i==14) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.05;
//if (i==14) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==14) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.6;
if (i==14) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.9;
if (i==14) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.3;
if (i==14) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.2;



if (i==15) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.65;
if (i==15) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.75;


if (i==15) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==15) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==15) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==15) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==15) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.39;
if (i==15) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.42;
if (i==15) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==15) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.17;

if (i==15) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.05;
//if (i==15) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==15) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.8;
if (i==15) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.6;
if (i==15) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==15) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.2;






if (i==16) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.8;

if (i==16) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==16) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==16) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==16) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.23;
if (i==16) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==16) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.52;
if (i==16) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==16) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==16) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.2;
//if (i==15) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==16) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.8;
if (i==16) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.6;
if (i==16) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==16) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.7;




if (i==17) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.83;
if (i==17) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][0][ialpha-1]*0.83;


if (i==17) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==17) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==17) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==17) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.23;
if (i==17) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==17) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.49;
if (i==17) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==17) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==17) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.2;
//if (i==15) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==17) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.8;
if (i==17) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.6;
if (i==17) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==17) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.5;



if (i==18) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.85;
if (i==18) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][0][ialpha-1]*0.85;

if (i==18) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==18) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.7;
if (i==18) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.7;
if (i==18) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.7;



if (i==18) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==18) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==18) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;

if (i==18) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.23;
if (i==18) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==18) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.5;
if (i==18) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==18) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==18) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.2;
if (i==18) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.9;
if (i==18) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.6;
if (i==18) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.25;
if (i==18) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.9;



if (i==19) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==19) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.5;
if (i==19) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.8;

if (i==19) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.23;
if (i==19) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==19) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==19) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==19) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==19) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.72;
if (i==19) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.55;
if (i==19) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.4;
if (i==19) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.37;


if (i==20) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==20) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.4;
if (i==20) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;

if (i==20) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.35;
if (i==20) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.45;
if (i==20) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==20) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==20) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==20) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.72;
if (i==20) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.55;
if (i==20) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.4;
if (i==20) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.37;







if (i==21) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][0][ialpha-1]*0.85;
if (i==21) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][0][ialpha-1]*0.85;
if (i==21) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][0][ialpha-1]*0.95;

if ((i>=14)&&(i<=21)) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.9;


if (i==21) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.85;
if (i==21) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.9;
if (i==21) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==21) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.85;


if (i==21) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.3;
if (i==21) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.4;
if (i==21) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.9;

if (i==21) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.35;
if (i==21) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.4;
if (i==21) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.5;
if (i==21) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.3;

if (i==21) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==21) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.95;
//if (i==21) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.72;
if (i==21) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*1.25;
if (i==21) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.3;
if (i==21) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*1.1;


//---------------------------------------------------------------------------


if (i==22) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==22) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.35;

if (i==22) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.45;
if (i==22) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;
if (i==22) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.4;
if (i==22) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.3;
if (i==22) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.4;
if (i==22) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.2;
if (i==22) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.1;
if (i==22) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==22) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.15;
if (i==22) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.1;





if (i==23) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==23) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==23) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.45;
if (i==23) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.6;

if (i==23) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.6;
if (i==23) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.8;


if (i==23) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.6;
if (i==23) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==23) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.2;
if (i==23) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==23) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.15;
if (i==23) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.1;



if (i==24) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.8;
if (i==24) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==24) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==24) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==24) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.7;

if (i==24) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.15;
if (i==24) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.48;
if (i==24) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.8;


if (i==24) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.7;
if (i==24) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.4;
if (i==24) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.2;
if (i==24) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.15;
if (i==24) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.15;
if (i==24) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.1;




if (i==25) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==25) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==25) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==25) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.25;


if (i==25) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==25) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==25) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.1;
if (i==25) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.7;
if (i==25) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.9;


if (i==25) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.7;
if (i==25) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.3;
if (i==25) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.1;
if (i==25) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*0.82;
if (i==25) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.05;
if (i==25) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.85;
if (i==25) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.65;




if (i==26) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==26) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==26) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==26) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.15;

if (i==26) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==26) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==26) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.1;
if (i==26) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.9;
if (i==26) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.8;


if (i==26) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.5;
if (i==26) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.32;
if (i==26) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==26) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.05;
if (i==26) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.1;
if (i==26) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.16;
if (i==26) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.98;

if (i==26) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.15;



if (i==27) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==27) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==27) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==27) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;

if (i==27) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==27) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==27) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.67;
if (i==27) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*2.;
if (i==27) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.75;


if (i==27) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.5;
if (i==27) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.42;
if (i==27) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.15;
if (i==27) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.1;
if (i==27) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.1;
if (i==27) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.12;
if (i==27) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*0.98;

if (i==27) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.08;


if (i==28) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.75;

if (i==28) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==28) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==28) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
if (i==28) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;

if (i==28) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==28) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==28) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*0.95;
if (i==28) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.7;
if (i==28) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*2.05;


if (i==28) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.75;
if (i==28) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.73;
if (i==28) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.49;
if (i==28) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.73;
if (i==28) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*1.6;
if (i==28) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*1.1;
if (i==28) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*1.45;

if (i==28) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*1.5;



//if (i==27) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==29) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.5;
if (i==29) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==29) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;


if (i==29) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==29) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==29) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*0.7;
//if (i==29) SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][14][itheta-1][ialpha-1]*0.3;





if (i==29) SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][1][itheta-1][ialpha-1]*0.35;
if (i==29) SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][2][itheta-1][ialpha-1]*0.87;

if (i==29) SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][3][itheta-1][ialpha-1]*1.;
if (i==29) SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][4][itheta-1][ialpha-1]*1.7;
if (i==29) SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][5][itheta-1][ialpha-1]*1.9;


if (i==29) SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][6][itheta-1][ialpha-1]*1.65;
if (i==29) SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][7][itheta-1][ialpha-1]*1.65;
if (i==29) SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][8][itheta-1][ialpha-1]*1.46;
if (i==29) SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][9][itheta-1][ialpha-1]*1.54;
if (i==29) SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][10][itheta-1][ialpha-1]*3.02;
if (i==29) SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][11][itheta-1][ialpha-1]*0.85;
if (i==29) SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][12][itheta-1][ialpha-1]*1.3;

if (i==29) SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][is23-1][13][itheta-1][ialpha-1]*2.;




};
};
};




for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=14; itheta++) {


if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.1;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.1;


if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.2;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.2;

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.1;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.1;


if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.1;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.11;

if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.05;
if (i==0) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.05;


if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.1;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.1;


if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.15;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.15;

if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;


if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.3;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.3;

if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.2;
if (i==2) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.2;



if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.1;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.1;

if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;


if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.3;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.3;

if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.2;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.2;


if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.9;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.9;

if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==3) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.9;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.9;

if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;

if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.1;
if (i==4) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.1;




if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.9;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.9;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.8;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.8;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.75;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.75;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.7;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.7;

if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==5) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.9;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.9;

if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.85;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.85;

if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;

if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.75;
if (i==6) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.75;

//--------------------------------------------
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.45;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.45;

if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.35;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.35;




if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.96;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.78;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;

if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;


if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.78;
if (i==8) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.96;


if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.65;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.65;

if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.45;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.45;

if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.03;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.78;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.75;


if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;


if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.75;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.79;
if (i==9) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.03;



if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;

if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;

if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.9;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.9;


if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.85;
if (i==10) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.85;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.92;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.89;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.89;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.92;


if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;

if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==11) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;


if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.92;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.92;


if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;


if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==12) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;



if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.05;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.05;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.95;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.95;



if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;


if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.2;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.2;

if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.9;
if (i==13) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.9;




//if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.95;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.92;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.03;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.92;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.92;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.03;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.92;
//if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.95;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;

if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.1;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.1;

if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.9;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.9;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.9;
if (i==14) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.9;



if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.03;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.03;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.08;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.08;


if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.85;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.85;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.95;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.95;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;


if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.65;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.65;


if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.65;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.65;

if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.75;
if (i==15) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.75;






if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.08;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.25;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.87;


if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.08;

if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.85;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.85;

if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.25;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.87;


if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.7;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.7;

if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.55;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.55;


if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.75;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.75;

if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==16) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;



if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.08;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.08;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.2;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.2;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.9;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.9;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.7;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.7;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.7;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.7;

if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.9;
if (i==17) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.9;


if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.9;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.9;



if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.2;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.2;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==18) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;





if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.1;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.1;



if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.5;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.5;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==19) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;





if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.1;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.1;




if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.5;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.5;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==20) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;



if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*1.1;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*1.1;


if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*2.5;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*2.5;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*2.;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*2.;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.45;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.45;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*1.05;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*1.05;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.7;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.7;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.8;
if (i==21) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.8;


//----------------------------------------------------------

if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.7;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.7;

if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==22) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;



if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.25;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.3;

if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.4;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.4;

if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==23) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;



if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.35;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.4;

if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.5;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.5;

if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==24) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;




if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.35;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.35;

if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.5;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.5;

if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==25) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;




if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.35;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.35;

if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.5;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.5;

if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==26) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;



if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.8;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.8;


if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.55;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.55;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.6;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.6;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.65;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.65;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][3]*0.75;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][10]*0.75;

if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*0.9;
if (i==27) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*0.9;


if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.35;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.35;


if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.5;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.5;

if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.8;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.8;

if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][0]*1.3;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][13]*1.3;

if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][1]*1.15;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][12]*1.15;


if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][2]*1.15;
if (i==28) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][11]*1.15;



if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][6]*0.55;
if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][7]*0.55;

if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][5]*0.75;
if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][8]*0.75;

if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][4]*0.9;
if (i==29) SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9] = SIGMA_ARR_GOL[i][is23-1][is12-1][itheta-1][9]*0.9;



};
};
};


for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=14; itheta++) {
for (Int_t ialpha = 1; ialpha <=14; ialpha++) {



if (i==0) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.1;
if (i==0) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.2;
if (i==0) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.7;
if (i==0) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.8;
if (i==0) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.05;

if (i==1) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.9;

if (i==1) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.2;
if (i==1) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.1;
if (i==1) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.85;
if (i==1) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.7;
if (i==1) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.8;
//if (i==1) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.05;
if (i==1) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.2;
if (i==1) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.2;



if (i==2) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==2) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==2) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*0.9;

if (i==2) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.8;
if (i==2) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.9;
if (i==2) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.1;


if (i==3) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.8;
if (i==3) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==3) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==3) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*0.9;
if (i==3) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.8;
if (i==3) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.93;
if (i==3) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.05;
if (i==3) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.1;

if (i==4) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==4) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==4) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*0.75;
if (i==4) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*0.95;


if (i==5) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.95;
if (i==5) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.05;
if (i==5) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.28;
if (i==5) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.06;
if (i==5) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.04;
if (i==5) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*0.96;
if (i==5) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.17;


if (i==6) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.85;
if (i==6) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.25;
if (i==6) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.25;



if (i==7) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.8;
if (i==7) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.15;
if (i==7) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.25;
if (i==7) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.25;



//----------------------------------------------

if (i==8) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.9;
if (i==8) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.7;

if (i==8) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;


if (i==8) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==8) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.17;
if (i==8) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.3;
if (i==8) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.35;
if (i==8) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.12;


if (i==9) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.8;
//if (i==9) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.97;
if (i==9) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.71;
if (i==9) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.84;

if (i==9) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.15;
if (i==9) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.45;
if (i==9) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.54;
if (i==9) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.57;
if (i==9) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.38;
if (i==9) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.2;
if (i==9) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.08;
if (i==9) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.8;




if (i==10) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.8;
if (i==10) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.89;
if (i==10) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.9;

if (i==10) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.13;
if (i==10) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==10) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.42;
if (i==10) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.45;
if (i==10) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.36;
if (i==10) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.27;
if (i==10) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.1;

//if (i==10) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.9;
if (i==10) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*0.9;

if (i==11) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.85;
if (i==11) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.92;
if (i==11) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.8;

if (i==11) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.8;
//if (i==11) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.13;
if (i==11) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==11) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.42;
if (i==11) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.48;
if (i==11) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.44;
if (i==11) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.3;
if (i==11) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.15;


if (i==12) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.4;
if (i==12) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.63;
if (i==12) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.85;
if (i==12) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.75;


//if (i==12) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.1;
if (i==12) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.18;
if (i==12) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.4;
if (i==12) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.45;
if (i==12) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.4;
if (i==12) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.3;
if (i==12) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.15;



if (i==13) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.6;
if (i==13) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.94;
if (i==13) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.83;
if (i==13) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.83;

if (i==13) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*0.99;
if (i==13) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.35;
if (i==13) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.6;
if (i==13) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.8;
if (i==13) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.8;
if (i==13) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.6;
if (i==13) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.4;
if (i==13) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.2;




if (i==14) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.3;
if (i==14) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.91;
if (i==14) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.91;
if (i==14) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.83;

if (i==14) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==14) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.45;
if (i==14) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.9;
if (i==14) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*2.17;
if (i==14) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*2.2;
if (i==14) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.6;
if (i==14) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.35;
if (i==14) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.2;


if (i==15) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*0.8;



if (i==15) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.3;
if (i==15) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.74;
if (i==15) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.9;
if (i==15) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.97;

if (i==15) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.24;
if (i==15) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.63;
if (i==15) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*2.02;
if (i==15) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*2.1;
if (i==15) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*2.0;
if (i==15) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.89;
if (i==15) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.95;
if (i==15) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.6;



if (i==16) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.7;
if (i==16) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==16) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.09;
if (i==16) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.15;

if (i==16) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.37;
if (i==16) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.53;
if (i==16) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.62;
if (i==16) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.7;
if (i==16) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.72;
if (i==16) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.63;
if (i==16) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.5;
if (i==16) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.4;






if (i==17) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.5;
if (i==17) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.5;
if (i==17) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.63;
if (i==17) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.81;
if (i==17) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*0.8;
if (i==17) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;



if (i==17) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*1.2;
if (i==17) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.38;
if (i==17) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.12;
//if (i==17) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.96;


if (i==17) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.22;
if (i==17) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==17) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.4;
if (i==17) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.5;
if (i==17) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.51;
if (i==17) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.41;
if (i==17) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.41;
if (i==17) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.4;

if (i==17) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.2;


if (i==18) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.8;
if (i==18) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*1.2;
if (i==18) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.25;
if (i==18) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.07;
if (i==18) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.09;


if (i==18) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.33;
if (i==18) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.58;
if (i==18) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.83;
if (i==18) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.9;
if (i==18) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.76;
if (i==18) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.51;
if (i==18) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.4;
if (i==18) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.24;

if (i==18) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.2;


if (i==19) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.1;
if (i==19) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.33;
if (i==19) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.42;

if (i==19) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.4;
if (i==19) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.35;
if (i==19) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.27;
if (i==19) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.26;

if (i==19) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.18;
if (i==19) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.1;
if (i==19) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*0.6;






if (i==20) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.7;
if (i==20) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.7;
if (i==20) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.8;
if (i==20) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*0.8;
if (i==20) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==20) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.9;





if (i==20) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*1.15;
if (i==20) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.2;
if (i==20) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.18;
if (i==20) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.22;



if (i==20) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.41;
if (i==20) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.54;
if (i==20) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.52;

if (i==20) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.55;
if (i==20) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.47;
if (i==20) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.37;
if (i==20) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.32;

if (i==20) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.2;
if (i==20) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.07;
if (i==20) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*0.6;





if (i==21) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.5;
if (i==21) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.7;
if (i==21) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.75;
if (i==21) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*0.8;
if (i==21) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*0.9;
if (i==21) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==21) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.9;





if (i==21) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*1.15;
if (i==21) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.15;
if (i==21) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.05;
if (i==21) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.18;


if (i==21) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.19;
if (i==21) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.4;
if (i==21) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.48;

if (i==21) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.56;
if (i==21) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.52;
if (i==21) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.45;
if (i==21) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.3;

if (i==21) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.1;
if (i==21) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.07;
if (i==21) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*0.6;

//------------------------------------

//if (i==22) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.95;
//if (i==22) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.02;
if (i==22) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.97;

if (i==22) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.25;
if (i==22) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.42;
if (i==22) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.55;

if (i==22) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.7;
if (i==22) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.55;
if (i==22) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.45;
if (i==22) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.3;



if (i==23) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.67;
if (i==23) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.87;
if (i==23) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.82;

if (i==23) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.15;

if (i==23) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.25;
if (i==23) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==23) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.23;

if (i==23) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.3;
if (i==23) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.65;
if (i==23) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*2.1;
if (i==23) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*2.4;

if (i==23) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.8;


//if (i==23) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.85;


if (i==24) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.58;
if (i==24) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.91;
if (i==24) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*0.97;

if (i==24) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.15;

if (i==24) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.26;
if (i==24) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.27;
if (i==24) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.2;

if (i==24) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.13;
if (i==24) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.33;
if (i==24) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.55;
if (i==24) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.65;

if (i==24) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.6;


//if (i==24) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*0.85;


if (i==25) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.45;
if (i==25) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.95;
if (i==25) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.05;

if (i==25) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.45;

if (i==25) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.45;
if (i==25) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.45;
if (i==25) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.3;

if (i==25) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.15;
if (i==25) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.25;
if (i==25) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.6;
if (i==25) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*2.0;

if (i==25) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.85;
if (i==25) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.3;



if (i==26) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.92;
if (i==26) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.8;
if (i==26) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.1;

if (i==26) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.53;

if (i==26) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.6;
if (i==26) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.7;
if (i==26) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.63;

if (i==26) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.4;
if (i==26) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.35;
if (i==26) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.6;
if (i==26) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.7;

if (i==26) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.6;
if (i==26) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.3;



if (i==27) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.91;
if (i==27) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.9;
if (i==27) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.04;

if (i==27) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.3;

if (i==27) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.28;
if (i==27) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.3;
if (i==27) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.3;

if (i==27) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.35;
if (i==27) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.3;
if (i==27) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.4;
if (i==27) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*1.7;

if (i==27) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*1.67;
if (i==27) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.3;



if (i==28) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.78;
//if (i==28) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*0.95;
if (i==28) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.13;

if (i==28) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.3;

if (i==28) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.5;
if (i==28) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*1.67;
if (i==28) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.75;

if (i==28) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.67;
if (i==28) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.6;
if (i==28) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*1.8;
if (i==28) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*2.6;

if (i==28) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*2.9;
if (i==28) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*1.75;

if (i==28) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*1.7;



if (i==29) SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][1][is12-1][itheta-1][ialpha-1]*0.64;
if (i==29) SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][2][is12-1][itheta-1][ialpha-1]*1.02;
if (i==29) SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][3][is12-1][itheta-1][ialpha-1]*1.33;

if (i==29) SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][4][is12-1][itheta-1][ialpha-1]*1.5;

if (i==29) SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][5][is12-1][itheta-1][ialpha-1]*1.65;
if (i==29) SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][6][is12-1][itheta-1][ialpha-1]*2.;
if (i==29) SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][7][is12-1][itheta-1][ialpha-1]*1.99;

if (i==29) SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][8][is12-1][itheta-1][ialpha-1]*1.9;
if (i==29) SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][9][is12-1][itheta-1][ialpha-1]*1.79;
if (i==29) SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][10][is12-1][itheta-1][ialpha-1]*2.15;
if (i==29) SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][11][is12-1][itheta-1][ialpha-1]*2.95;

if (i==29) SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][12][is12-1][itheta-1][ialpha-1]*3.0;
if (i==29) SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][13][is12-1][itheta-1][ialpha-1]*2.;

if (i==29) SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL[i][14][is12-1][itheta-1][ialpha-1]*1.9;



};
};
};


};//end loop-golovach


for (Short_t i=0;i<=9;i++){
S12_ARR_FED_THRESH[i][0] = 0.07791914 + i*(0.089538593-0.07791914)/9.;
S12_ARR_FED_THRESH[i][1] = 0.07791914 + i*(0.105125093-0.07791914)/9.;
S12_ARR_FED_THRESH[i][2] = 0.07791914 + i*(0.121961593-0.07791914)/9.;

S23_ARR_FED_THRESH[i][0] = 1.161744 + i*(1.205450285 -1.161744 )/9.;
S23_ARR_FED_THRESH[i][1] = 1.161744 + i*(1.260971785 -1.161744 )/9.;
S23_ARR_FED_THRESH[i][2] = 1.161744 + i*(1.317743285 -1.161744 )/9.;

};

W_ARR_FED_THRESH[0] = 1.2375;
W_ARR_FED_THRESH[1] = 1.2625;
W_ARR_FED_THRESH[2] = 1.2875;

//FEDOTOV CROSS SECTION READING 
THETA_ARR_FED[0] = 0.;
THETA_ARR_FED[1] =0.4559418;
THETA_ARR_FED[2] =0.9018836;
THETA_ARR_FED[3] =1.347825;
THETA_ARR_FED[4] =1.793767;
THETA_ARR_FED[5] =2.239709;
THETA_ARR_FED[6] =2.685651;
THETA_ARR_FED[7] =3.141593;

ALPHA_ARR_FED[0] = 0.;
ALPHA_ARR_FED[1] = 0.9047408;
ALPHA_ARR_FED[2] = 1.799482;
ALPHA_ARR_FED[3] = 2.694222;
ALPHA_ARR_FED[4] = 3.588963;
ALPHA_ARR_FED[5] = 4.483704;
ALPHA_ARR_FED[6] = 5.378445;
ALPHA_ARR_FED[7] = 6.283185;


for (Short_t i=0; i<=5; i++) {
for (Int_t is23 = 1; is23 <=10; is23++) {
for (Int_t is12 = 1; is12 <=10; is12++) {
for (Int_t itheta = 1; itheta <=8; itheta++) {
for (Int_t ialpha = 1; ialpha <=8; ialpha++) {

SIGMA_ARR_FED[i][0][0][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][1][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][2][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][3][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][4][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][5][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][6][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][0][7][is23-1][is12-1][itheta-1][ialpha-1] = 0.;

SIGMA_ARR_FED[i][1][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;

SIGMA_ARR_FED[i][2][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][2][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;


SIGMA_ARR_FED[i][3][9][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][3][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][3][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;

SIGMA_ARR_FED[i][4][8][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][4][9][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][4][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][4][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;


SIGMA_ARR_FED[i][5][8][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][5][9][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][5][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][5][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;



SIGMA_ARR_FED[i][6][6][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][7][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][8][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][9][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][10][is23-1][is12-1][itheta-1][ialpha-1] = 0.;
SIGMA_ARR_FED[i][6][11][is23-1][is12-1][itheta-1][ialpha-1] = 0.;


};
};
};
};
};


cout<<"Reading Fedotov  cross sections \n";

for (Short_t i=0; i<=55; i++) {

if ((i>=0)&&(i<=3)){
q2bin = 0;
wbin = 8+i;

};

if ((i>=4)&&(i<=14)){
q2bin = 1;
wbin = i-4;
};

if ((i>=15)&&(i<=24)){
q2bin = 2;
wbin = i-15;
};

if ((i>=25)&&(i<=33)){
q2bin = 3;
wbin = i-25;
};

if ((i>=34)&&(i<=41)){
q2bin = 4;
wbin = i-34;
};

if ((i>=42)&&(i<=49)){
q2bin = 5;
wbin = i-42;
};

if ((i>=50)&&(i<=55)){
q2bin = 6;
wbin = i-50;
};


if ((q2bin >=0)&&(q2bin <=2)) Q2_ARR_FED[q2bin] = 0.225+0.05*q2bin;

if ((q2bin >=3)&&(q2bin <=6)) Q2_ARR_FED[q2bin] = 0.225+0.05*(q2bin+1);

W_ARR_FED[wbin] = 1.3125+0.025*wbin;
string dummy,xsect;

string file=file_names_fed[i];
ifstream input(file.c_str());
if (input.fail()) cout << "ALARM! File for Q2 = "<< Q2_ARR_FED[q2bin]<<" GeV2 and W = "<<W_ARR_FED[wbin]<<" GeV with Fedotov cr sect FAILED to open! \n";
if(input.is_open()){

for (Int_t is23 = 1; is23 <=10; is23++) {
for (Int_t is12 = 1; is12 <=10; is12++) {
for (Int_t itheta = 1; itheta <=8; itheta++) {
for (Int_t ialpha = 1; ialpha <=8; ialpha++) {
getline(input,xsect);
//Define 2dim s12 and s23 arrays
S12_ARR_FED[is12-1][wbin] = atof(xsect.c_str());
//cout << S12_ARR_FED[is12-1][wbin] <<"\n";
getline(input,xsect);
S23_ARR_FED[is23-1][wbin] = atof(xsect.c_str());
//cout << S23_ARR_FED[is12-1][wbin] <<"\n";
getline(input,dummy);
getline(input,dummy);

getline(input,xsect);
//sigma_t
SIGMA_ARR_FED[0][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_l
SIGMA_ARR_FED[1][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_c2f
SIGMA_ARR_FED[2][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_s2f
SIGMA_ARR_FED[3][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_cf
SIGMA_ARR_FED[4][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_sf
SIGMA_ARR_FED[5][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,dummy);

getline(input,dummy);



};
};
};
};








};
input.close();

for (Int_t is23 = 1; is23 <=10; is23++) {
for (Int_t is12 = 1; is12 <=10; is12++) {
for (Int_t itheta = 1; itheta <=8; itheta++) {
for (Int_t ialpha = 1; ialpha <=8; ialpha++) {



for (Short_t j=0;j<6;j++){

SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_FED[9][wbin]-S12_ARR_FED[0][wbin])*(S23_ARR_FED[9][wbin]-S23_ARR_FED[0][wbin]);



SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*0.1;
SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*0.4;
SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*0.7;




SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*M_PI*2.*M_PI;



if ((q2bin==0)&&(wbin==10)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.97;
if ((q2bin==0)&&(wbin==9)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.06;

if ((q2bin==1)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.80;
if ((q2bin==1)&&(wbin==7)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.93;
if ((q2bin==1)&&(wbin==8)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.93;
if ((q2bin==1)&&(wbin==9)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.93;
if ((q2bin==1)&&(wbin==10)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.98;



if ((q2bin==2)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.25;
if ((q2bin==2)&&(wbin==8)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.05;
if ((q2bin==2)&&(wbin==9)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.05;


if ((q2bin==3)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.1;
if ((q2bin==3)&&(wbin==7)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.96;

if ((q2bin==4)&&(wbin==7)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.97;

if ((q2bin==6)&&(wbin==5)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.92;


if ((q2bin==2)&&(wbin==4)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.02;
if ((q2bin==2)&&(wbin==2)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.97;


if ((q2bin==3)&&(wbin==4)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.98;
if ((q2bin==3)&&(wbin==2)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.08;

if ((q2bin==4)&&(wbin==3)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.05;


if ((q2bin==5)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.05;
if ((q2bin==5)&&(wbin==3)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.96;
if ((q2bin==5)&&(wbin==2)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.87;


//if ((q2bin==6)&&(wbin==0)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.95;
//if ((q2bin==6)&&(wbin==3)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.07;
if ((q2bin==6)&&(wbin==4)) SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED[j][q2bin][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.08;

//---------------------------------------------------

//SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_FED_THRESH[9][0]-S12_ARR_FED_THRESH[0][0])*(S23_ARR_FED_THRESH[9][0]-S23_ARR_FED_THRESH[0][0]);

//SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*M_PI*2.*M_PI;



//SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_FED_THRESH[9][1]-S12_ARR_FED_THRESH[0][1])*(S23_ARR_FED_THRESH[9][1]-S23_ARR_FED_THRESH[0][1]);

//SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][1][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*M_PI*2.*M_PI;


//SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_FED_THRESH[9][2]-S12_ARR_FED_THRESH[0][2])*(S23_ARR_FED_THRESH[9][2]-S23_ARR_FED_THRESH[0][2]);

//SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][2][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*M_PI*2.*M_PI;





//SIGMA_ARR_FED_THRESH[j][q2bin][0][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]*1.1;






//cout << SIGMA_ARR_FED_THRESH[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]  << " "<< SIGMA_ARR_FED[j][q2bin][0][is23-1][is12-1][itheta-1][ialpha-1]<<"\n";
};

};
};
};
};



};


cout<< "Reading model    cross sections at Q2 = 1.3 GeV2 \n";
for (Short_t i=0; i<=20; i++) {
wbin = i;

if ((i>=0)&&(i<=12))  W_ARR_RIP2[i] = 1.8375 + 0.025*i;
if ((i>=13)&&(i<=20)) W_ARR_RIP2[i] = 2.1875 + 0.05*(i-13);
 

string dummy,xsect;

string file=file_names_rip2[i];
ifstream input(file.c_str());
if (input.fail()) cout << "ALARM! File for W = "<< W_ARR_RIP2[i]<< " GeV with model cr sect at Q2 = 1.3 GeV FAILED to open! \n";
if(input.is_open()){

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
getline(input,xsect);
//Define 2dim s12 and s23 arrays
S12_ARR_RIP2[is12-1][wbin] = atof(xsect.c_str());
getline(input,xsect);
S23_ARR_RIP2[is23-1][wbin] = atof(xsect.c_str());
getline(input,dummy);
getline(input,dummy);
getline(input,xsect);
//sigma_t
SIGMA_ARR_RIP2[0][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_l
SIGMA_ARR_RIP2[1][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_c2f
SIGMA_ARR_RIP2[2][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_s2f
SIGMA_ARR_RIP2[3][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_cf
SIGMA_ARR_RIP2[4][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_sf
SIGMA_ARR_RIP2[5][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
EPS_L_RIP2[wbin] = atof(xsect.c_str());
getline(input,dummy);



};
};
};
};
};
input.close();







for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {



for (Short_t j=0;j<6;j++){


SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_RIP2[11][wbin]-S12_ARR_RIP2[0][wbin])*(S23_ARR_RIP2[11][wbin]-S23_ARR_RIP2[0][wbin]);

SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*M_PI*2.*M_PI;

//SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*2.*M_PI;

SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.1;

if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.15*1.02*1.05*0.98*1.075;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.007*1.105;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.1*1.01*0.95*0.9*1.09;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.96*1.03*1.05*0.95;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.94*1.02*0.85;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.02*0.82;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.03*1.02*0.94;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.06*1.02*0.9;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.07*1.04*0.9;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.03*0.83*0.95;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.16*1.07*1.05;


if ((i==3)&&(j==1)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.8;
if ((i==4)&&(j==1)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.1;

if ((i>=11)&&(i<=20)&&(j==0)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.69;
if ((i>=11)&&(i<=20)&&(j==1)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*8.9;

if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.02;

if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.15;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.8;

if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.03;
if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.1;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.1;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.25;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.25;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.35;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.67;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*2.1;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*2.9;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*3.5;

//adjusting c2f 11/27/2020
if ((i==0)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.6;
if ((i==1)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.6;
if ((i==2)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.4;

if ((i==3)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.7;
if ((i==5)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*3.7;
if ((i==6)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.8;
if ((i==7)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.2;

if ((i==8)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.5*SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if ((i==9)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.15*SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if ((i==10)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.15*SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if ((i==11)&&(j==2)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.15*SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
};

//if ((i>=11)&&(i<=20)) SIGMA_ARR_RIP2[0][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[0][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.69;
//if ((i>=11)&&(i<=20)) SIGMA_ARR_RIP2[1][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[1][wbin][is23-1][is12-1][itheta-1][ialpha-1]*8.9;



};
};
};
};


for (Short_t j=0;j<6;j++){

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {


if ((i>=0)&&(i<=10)) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.5;

if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*1.1;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*1.05;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.3;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.35;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.1;


if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.6;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*0.55;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*0.55;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.55;
//if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*1.3;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.85;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.8;

if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.85;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.2;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.8;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.8;


if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.3;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.3;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.8;
//if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.7;
//if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.5;


if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*1.2;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.2;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.2;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.2;

if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.4;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.4;
//if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.05;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.85;

//if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.6;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.3;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.3;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.9;





if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.85;
//if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*16;
//if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.4;
//if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.7;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.75;

if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.3;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*1.2;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.75;

if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.5;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.1;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.9;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*1.1;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.9;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.75;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.9;

if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.75;
if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.65;


if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.7;
if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.6;

if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.75;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.6;

if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.9;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.9;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.9;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.9;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.9;



if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*0.65;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.9;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.75;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.75;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.45;

if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.8;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.8;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.7;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.7;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.7;





if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*1.2;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.95;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.8;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.75;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.45;

if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.9;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.9;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.9;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.9;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.9;



if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.7;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*1.1;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.7;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.9;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.8;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.9;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.9;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.9;


if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.9;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.9;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*1.1;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.9;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.95;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.8;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.8;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.8;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.9;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.8;


if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.8;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.8;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.7;
//if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.65;

if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.6;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.8;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.7;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.6;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.6;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.4;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.4;




//if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.25;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.9;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.8;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.8;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.7;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.6;


if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*0.7;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.85;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.6;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.5;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.4;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.4;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.3;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.3;






if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][2][itheta-1][ialpha-1]*0.8;
//if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.2;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.9;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.8;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.8;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.7;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.6;

if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.85;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.6;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.5;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.4;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.4;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.3;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.3;


};
};
};


for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {

if ((i>=0)&&(i<=10)) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.7;



if (i==0) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.2;
if (i==0) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.05;
if (i==0) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.75;
if (i==0) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.95;
//if (i==0) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*1.2;
if (i==0) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.15;
if (i==0) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.34;
if (i==0) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.25;


//if (i==1) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==1) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.9;
if (i==1) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.8;
if (i==1) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.8;
//if (i==1) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.3;
if (i==1) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.25;


//if (i==2) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.95;
//if (i==2) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.3;
if (i==2) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.5;
if (i==2) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==2) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.85;
if (i==2) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*1.25;
if (i==2) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.35;
if (i==2) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.55;
if (i==2) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.35;


if (i==3) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.1;
if (i==3) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.7;
if (i==3) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.6;

if (i==3) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.7;
if (i==3) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.76;
if (i==3) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.3;
if (i==3) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.7;
if (i==3) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.5;

if (i==4) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.23;
if (i==4) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.13;
if (i==4) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.75;
if (i==4) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.7;
if (i==4) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.95;
if (i==4) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*1.03;
if (i==4) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.15;
if (i==4) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.25;
if (i==4) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.25;

if (i==5) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.3;
if (i==5) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.7;
if (i==5) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.55;
if (i==5) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.7;
if (i==5) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.8;

if (i==5) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.45;
if (i==5) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.65;
if (i==5) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.55;



if (i==6) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
//if (i==6) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.1;
if (i==6) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.6;
if (i==6) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==6) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.75;
if (i==6) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.8;
if (i==6) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.35;
if (i==6) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.5;


if (i==7) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.65;
if (i==7) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.55;
if (i==7) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==7) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.8;
if (i==7) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.2;
if (i==7) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.5;
if (i==7) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.5;

if (i==8) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.95;
//if (i==8) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.1;
if (i==8) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.7;
if (i==8) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==8) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.7;
if (i==8) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.87;
if (i==8) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.03;
if (i==8) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.27;
if (i==8) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.6;


if (i==9) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.9;
if (i==9) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.6;
if (i==9) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.55;
if (i==9) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.7;
if (i==9) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.87;
if (i==9) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.4;
if (i==9) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.75;


if (i==10) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.75;
if (i==10) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.8;
if (i==10) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.6;
if (i==10) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.57;
if (i==10) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.7;
if (i==10) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.87;
if (i==10) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.13;
if (i==10) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.55;
if (i==10) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.75;



if (i==11) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.95;
if (i==11) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.75;
if (i==11) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*0.7;
if (i==11) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.6;
if (i==11) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.6;
if (i==11) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.6;
if (i==11) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.7;
if (i==11) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.7;
if (i==11) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.8;

if (i==12) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.95;
if (i==12) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*0.95;
if (i==12) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.05;
if (i==12) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==12) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.9;




if (i==13) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.95;
if (i==13) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.9;
if (i==13) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.9;
if (i==13) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.85;
if (i==13) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.95;








if (i==14) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==14) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.9;
if (i==14) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.85;
if (i==14) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.8;
//if (i==14) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.9;
//if (i==14) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.2;



if (i==15) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;

if (i==15) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.9;
if (i==15) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.85;
if (i==15) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.7;
if (i==15) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*1.05;
if (i==15) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.1;


if (i==16) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==16) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.9;
if (i==16) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.82;
if (i==16) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.7;
if (i==16) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.97;
if (i==16) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.2;
if (i==16) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.2;


if (i==17) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==17) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*0.9;
if (i==17) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.8;
if (i==17) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.7;
if (i==17) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.9;
if (i==17) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.2;
if (i==17) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.2;




if (i==18) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==18) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.35;
if (i==18) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.25;

if (i==18) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.05;
if (i==18) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.8;
if (i==18) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.65;
if (i==18) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.75;
if (i==18) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.2;
if (i==18) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.2;



if (i==19) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==19) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.45;
if (i==19) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.3;
if (i==19) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.1;
if (i==19) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.85;
if (i==19) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.65;
if (i==19) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.6;
if (i==19) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.05;
if (i==19) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.05;

if (i==20) SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;

if (i==20) SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.55;
if (i==20) SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.3;
if (i==20) SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.1;
if (i==20) SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.85;
if (i==20) SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.65;
if (i==20) SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.6;
//if (i==20) SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][8][is12-1][itheta-1][ialpha-1]*1.05;
if (i==20) SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][9][is12-1][itheta-1][ialpha-1]*1.15;


};
};
};


for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
if ((i>=0)&&(i<=10)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.2;
if ((i>=0)&&(i<=10)) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*0.75;


if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*0.85;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.1;


if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*2.5;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.855;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.7;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*1.05;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.6;

if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*1.25;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.05;
//if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.95;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.05;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.8;


if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.1;
//if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.9;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.1;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.8;


if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.1;
//if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.85;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.05;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.7;


if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*1.5;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.2;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.1;
//if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.9;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.8;


if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.1;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.9;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.1;


if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*1.5;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.25;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.1;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.9;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.1;

if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.15;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.1;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.95;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.9;


if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.25;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.85;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.2;


if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.2;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.05;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.85;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][4][ialpha-1]*1.2;



if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.25;
if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.25;
if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.75;


if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*0.95;
if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.25;
if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.25;
if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.95;
if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.65;


if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*0.85;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.65;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.7;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.95;
if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.65;

if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*0.85;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.8;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.7;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.95;
if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.65;


if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*0.85;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.9;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.7;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.95;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.55;

if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*0.75;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.95;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.7;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][3][ialpha-1]*0.95;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.5;


if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*0.6;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.95;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.75;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.5;



if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*0.57;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*1.95;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.75;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.5;

if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*0.55;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*2.05;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.75;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.5;

if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][0][ialpha-1]*0.55;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][1][ialpha-1]*2.05;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][2][ialpha-1]*1.75;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][5][ialpha-1]*0.5;
};
};
};





for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {


SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0]*0.75;
SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5]*0.75;

SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.05;
SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.05;


SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.15;
SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.15;


if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0]*0.85;
if (i==0) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5]*0.85;


if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0]*0.7;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5]*0.7;

//if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*0.8;
//if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*0.8;

if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==1) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;



if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;

if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0]*0.7;
if (i==2) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5]*0.7;

if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==3) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;


if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.6;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.6;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==4) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;

if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.6;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.6;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==5) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;


if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==6) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;

if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==7) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;


if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.5;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.5;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==8) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;


if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.7;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.7;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.4;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.4;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][0]*0.8;
if (i==9) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][5]*0.8;

if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*1.4;
if (i==10) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*1.4;


if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][2]*0.75;
if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][3]*0.75;

if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.2;
if (i==11) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.2;


//if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*0.9;
//if (i==12) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*0.9;


//if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*0.9;
//if (i==13) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*0.9;

//if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*0.95;
//if (i==14) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*0.92;

if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.05;
if (i==15) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.05;


if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.05;
if (i==16) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.05;

if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.05;
if (i==17) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.05;

if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.05;
if (i==18) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.05;

if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.05;
if (i==19) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.05;


if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][1]*1.05;
if (i==20) SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP2[j][wbin][is23-1][is12-1][itheta-1][4]*1.05;


};
};
};




};




};//end lop rip2 - i



for (Short_t i=0; i<=9; i++) {
wbin = i;


 W_ARR_RIP3[i] = 2.5875 + 0.05*i;
 
 



//cout << "Reading RIPANI diff cross sections for Q^2 = 1.3 GeV^2, W = "<< W_ARR_RIP3[i] << " GeV \n";

string dummy,xsect;

string file=file_names_rip3[i];
ifstream input(file.c_str());
if (input.fail()) cout << "ALARM! File for W = " << W_ARR_RIP3[i]<<" GeV with model cr sect at Q2 = 1.3 GeV2 FAILED to open! \n";
if(input.is_open()){

for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
getline(input,xsect);
//Define 2dim s12 and s23 arrays
S12_ARR_RIP3[is12-1][wbin] = atof(xsect.c_str());

getline(input,xsect);
S23_ARR_RIP3[is23-1][wbin] = atof(xsect.c_str());
getline(input,dummy);
getline(input,dummy);
getline(input,xsect);
//sigma_t
SIGMA_ARR_RIP3[0][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
// SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_l
SIGMA_ARR_RIP3[1][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_c2f
SIGMA_ARR_RIP3[2][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_s2f
SIGMA_ARR_RIP3[3][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_cf
SIGMA_ARR_RIP3[4][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);
//sigma_sf
SIGMA_ARR_RIP3[5][wbin][is23-1][is12-1][itheta-1][ialpha-1] = atof(xsect.c_str());
getline(input,xsect);

getline(input,dummy);



};
};
};
};
};
input.close();







for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {



for (Short_t j=0;j<6;j++){


SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*(S12_ARR_RIP3[15][wbin]-S12_ARR_RIP3[0][wbin])*(S23_ARR_RIP3[15][wbin]-S23_ARR_RIP3[0][wbin]);

SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*M_PI*M_PI*2.*M_PI;



if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.95;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.87;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.85;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.8;
if (i==6) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.6;
if (i==7) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.6*1.02;
if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.7*1.05;
if (i==9) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*0.7*1.18;

if (j==0) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.69*0.87;//mod 11/27/2020
if (j==1) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*8.9*0.45;//mod 11/27/2020

//adjusting sigma l 11/27/2020
if ((i==5)&&(j==1)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.15;
if ((i==6)&&(j==1)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*1.15;

//Further integral adjustment is on line 4260
};

};
};
};
};






for (Short_t j=0;j<6;j++){

for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {


if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.9;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.75;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.7;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.7;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.65;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.6;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.6;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.5;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*0.7;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1]*0.6;



if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.9;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.7;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.7;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.65;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.6;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.6;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.5;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*0.7;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*0.6;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1]*0.6;


if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][2][itheta-1][ialpha-1]*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.85;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.7*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.65*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.65*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.6*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.55*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.55*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.45*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*0.65*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*0.5*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1]*0.55;




if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][2][itheta-1][ialpha-1]*0.9;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.8;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.65*0.85;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.6*0.85;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.6*0.85;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.55*0.85;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.5*0.85;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.45*0.85;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.4*0.85;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*0.6*0.85;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*0.45*0.85;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1]*0.5;



if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][2][itheta-1][ialpha-1]*0.9;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.85;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.7;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.65;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.65;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.6;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.55;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.55;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.45;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*0.65;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1]*0.55;


if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][1][itheta-1][ialpha-1]*0.7;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][2][itheta-1][ialpha-1]*0.9;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][4][itheta-1][ialpha-1]*0.85;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][5][itheta-1][ialpha-1]*0.7;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][6][itheta-1][ialpha-1]*0.65;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][7][itheta-1][ialpha-1]*0.65;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][8][itheta-1][ialpha-1]*0.6;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][9][itheta-1][ialpha-1]*0.55;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][10][itheta-1][ialpha-1]*0.55;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.45;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*0.65;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*0.5;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1]*0.55;



if (i==6) SIGMA_ARR_RIP3[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.15;
if (i==6) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.95;
if (i==6) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*1.5;
if (i==6) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*1.15;

if (i==7) SIGMA_ARR_RIP3[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.15;
if (i==7) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.95;
if (i==7) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*1.5;
if (i==7) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*1.15;


if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.15;
if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.95;
if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*1.5;
if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*1.15;
if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1]*0.8;


if (i==9) SIGMA_ARR_RIP3[j][wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][3][itheta-1][ialpha-1]*1.15;
if (i==9) SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][11][itheta-1][ialpha-1]*0.95;
if (i==9) SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][12][itheta-1][ialpha-1]*1.5;
if (i==9) SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][13][itheta-1][ialpha-1]*1.15;
if (i==9) SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][14][itheta-1][ialpha-1]*0.8;



};
};
};


for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {

SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.8;

if (i==0) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.3;
if (i==0) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.3;
if (i==0) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.15;
//if (i==0) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==0) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.85;
if (i==0) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.7;
if (i==0) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.6;
if (i==0) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.5;
if (i==0) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.5;

if (i==0) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.75;
if (i==0) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.6;
if (i==0) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.4;



if (i==1) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.3;
if (i==1) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.3;
if (i==1) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.15;
//if (i==1) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==1) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.85;
if (i==1) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.7;
if (i==1) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.6;
if (i==1) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.5;
if (i==1) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.5;

if (i==1) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.75;
if (i==1) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.6;
if (i==1) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.4;

if (i==2) SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.8;
if (i==2) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.4;
if (i==2) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.35;
if (i==2) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.2;
//if (i==2) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==2) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==2) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.75;
if (i==2) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.65;
if (i==2) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.55;
if (i==2) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.45;

if (i==2) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.55;
if (i==2) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.45;
if (i==2) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.4;


if (i==3) SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.8;
if (i==3) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.4;
if (i==3) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.35;
if (i==3) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.2;
//if (i==3) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==3) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==3) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.75;
if (i==3) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.65;
if (i==3) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.55;
if (i==3) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.45;

if (i==3) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.55;
if (i==3) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.45;
if (i==3) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.7;
if (i==3) SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1]*1.4;

if (i==4) SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.75;
if (i==4) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.4;
if (i==4) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.35;
if (i==4) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.2;
//if (i==4) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==4) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==4) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.75;
if (i==4) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.65;
if (i==4) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.55;
if (i==4) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.45;

if (i==4) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.55;
if (i==4) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.3;
if (i==4) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.7;
if (i==4) SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1]*1.4;






if (i==5) SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.85;
if (i==5) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.4;
if (i==5) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.35;
if (i==5) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.2;
//if (i==5) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*1.05;
if (i==5) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.9;
if (i==5) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.75;
if (i==5) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.65;
if (i==5) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.55;
if (i==5) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.45;
if (i==5) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.4;
if (i==5) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.1;
if (i==5) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.7;
if (i==5) SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1]*1.4;



if (i==6) SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.85;
if (i==6) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.4;
if (i==6) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.35;
if (i==6) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.2;
if (i==6) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.9;
if (i==6) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.8;
if (i==6) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.75;
if (i==6) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.65;
if (i==6) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.55;
if (i==6) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.45;
if (i==6) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.4;
if (i==6) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.1;
if (i==6) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.7;
if (i==6) SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1]*1.4;


if (i==7) SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==7) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.4;
if (i==7) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.35;
if (i==7) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.2;
if (i==7) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.9;
if (i==7) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.8;
if (i==7) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.75;
if (i==7) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.65;
if (i==7) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.55;
if (i==7) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.5;
if (i==7) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.5;
if (i==7) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.1;
if (i==7) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.7;
if (i==7) SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1]*1.4;

if (i==8) SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==8) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.4;
if (i==8) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.35;
if (i==8) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.2;
if (i==8) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.9;
if (i==8) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.8;
if (i==8) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.75;
if (i==8) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.65;
if (i==8) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.55;
if (i==8) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.5;
if (i==8) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.5;
if (i==8) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.1;
if (i==8) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.7;
if (i==8) SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1]*1.4;



if (i==9) SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][1][is12-1][itheta-1][ialpha-1]*0.9;
if (i==9) SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][2][is12-1][itheta-1][ialpha-1]*1.4;
if (i==9) SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][3][is12-1][itheta-1][ialpha-1]*1.35;
if (i==9) SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][4][is12-1][itheta-1][ialpha-1]*1.2;
if (i==9) SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][5][is12-1][itheta-1][ialpha-1]*0.9;
if (i==9) SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][6][is12-1][itheta-1][ialpha-1]*0.8;
if (i==9) SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][7][is12-1][itheta-1][ialpha-1]*0.75;
if (i==9) SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][8][is12-1][itheta-1][ialpha-1]*0.65;
if (i==9) SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][9][is12-1][itheta-1][ialpha-1]*0.55;
if (i==9) SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][10][is12-1][itheta-1][ialpha-1]*0.5;
if (i==9) SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][11][is12-1][itheta-1][ialpha-1]*0.5;
if (i==9) SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][12][is12-1][itheta-1][ialpha-1]*1.1;
if (i==9) SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][13][is12-1][itheta-1][ialpha-1]*1.7;
if (i==9) SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][14][is12-1][itheta-1][ialpha-1]*1.4;





};
};
};



for (Int_t is12 = 1; is12 <=12; is12++) {
for (Int_t is23 = 1; is23 <=12; is23++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {




SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1]*0.5;
SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1]*2.4;
SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][4][ialpha-1]*1.4;
SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][5][ialpha-1]*0.5;



if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*1.8;

if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*1.8;

if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1]*1.05;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1]*1.3;
if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*1.8;

if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1]*1.05;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1]*1.5;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*2.1;


if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1]*1.05;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1]*1.7;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*2.1;

if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1]*1.05;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1]*1.7;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*2.1;


if (i==6) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1]*1.05;
if (i==6) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1]*1.7;
if (i==6) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*1.9;

if (i==7) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1]*1.1;
if (i==7) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1]*2.3;
if (i==7) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*1.9;


if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1]*1.1;
if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1]*2.3;
if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*1.9;


if (i==9) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][0][ialpha-1]*0.9;
if (i==9) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][1][ialpha-1]*2.4;
if (i==9) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][2][ialpha-1]*1.9;



};
};
};


for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {

SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][2]*1.15;
SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][3]*1.15;

SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][1]*1.05;
SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][4]*1.05;

SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][0]*0.75;
SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][5]*0.75;

if (i>2) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][0]*0.9;
if (i>2) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][5]*0.9;
};
};
};


};



for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {

SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.805*SIGMA_ARR_RIP3[0][wbin][is23-1][is12-1][itheta-1][ialpha-1];//11/27/2020

if (wbin==0) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*11.3*0.95;
if (wbin==1) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*14.1*0.95;
if (wbin==2) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*12.7*0.95;
if (wbin==3) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*17.9*0.95;
if (wbin==4) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*14.*0.95;
if (wbin==5) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*15.*0.95;
if (wbin==6) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*16.9*0.95;
if (wbin==7) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*14.8*0.95;
if (wbin==8) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*15.4*0.95;
if (wbin==9) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][ialpha-1]*19.*0.9;

};
};
};
};


for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {

//Initial integral adjustment is on line 3800
for (Short_t j=0;j<6;j++){
if (i==0) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*3.45*1.1*0.8;
if (i==1) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*3.8*1.1*0.8;

if (i==2) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.3*1.1*0.92;
if (i==3) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.6*1.1;
if (i==4) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.5*1.1;
if (i==5) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.5*1.1;
if (i==6) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.5*1.1;

//mod 11/27/2020
if ((i==7)&&(j==0)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.5;
if ((i==7)&&!(j==0)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.5*1.1;


if (i==8) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.33*1.1;
//mod 11/27/2020
if ((i==9)&&(j==0)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.2*1.4;
if ((i==9)&&!(j==0)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1]*4.2*1.1;

//adjusting sigma c2f 11/27/2020
if ((i==2)&&(j==2)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.25*SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if ((i==3)&&(j==2)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.25*SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if ((i==4)&&(j==2)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.25*SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if ((i==5)&&(j==2)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.25*SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if ((i==6)&&(j==2)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.1*SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if ((i==7)&&(j==2)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = -SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if ((i==8)&&(j==2)) SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = -SIGMA_ARR_RIP3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];

};
};
};
};
};



for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {


SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*1.5;
SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*3.;
SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*3.;
SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*3.3;
SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*3.5;
SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*3.7;
SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*3.9;
SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*4.1;
SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*4.2;
SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*4.;
SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.45;


if (wbin==0) SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.31;
if (wbin==0) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.59;
if (wbin==0) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.6;
if (wbin==0) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.1;
if (wbin==0) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.33;
if (wbin==0) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==0) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.7;
if (wbin==0) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.92;
if (wbin==0) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*1.97;
if (wbin==0) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*1.9;
if (wbin==0) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.07;
if (wbin==0) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==0) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;



if (wbin==1) SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.32;
if (wbin==1) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.55;
if (wbin==1) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.62;
if (wbin==1) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.1;
if (wbin==1) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.35;
if (wbin==1) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==1) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.8;
if (wbin==1) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.92;
if (wbin==1) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*1.97;
if (wbin==1) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*1.9;
if (wbin==1) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.07;
if (wbin==1) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==1) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;


if (wbin==2) SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.32;
if (wbin==2) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.55;
if (wbin==2) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.62;
if (wbin==2) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.1;
if (wbin==2) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.35;
if (wbin==2) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==2) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.8;
if (wbin==2) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.92;
if (wbin==2) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*1.97;
if (wbin==2) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*1.9;
if (wbin==2) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.15;
if (wbin==2) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.45;
if (wbin==2) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;




if (wbin==3) SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.32;
if (wbin==3) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.55;
if (wbin==3) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.62;
if (wbin==3) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.1;
if (wbin==3) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.35;
if (wbin==3) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==3) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.8;
if (wbin==3) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.92;
if (wbin==3) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*2.0;
if (wbin==3) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*2.2;
if (wbin==3) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==3) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==3) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;




if (wbin==4) SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.32;
if (wbin==4) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.55;
if (wbin==4) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.62;
if (wbin==4) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.1;
if (wbin==4) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.35;
if (wbin==4) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==4) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.8;
if (wbin==4) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.92;
if (wbin==4) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*2.0;
if (wbin==4) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*2.2;
if (wbin==4) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==4) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==4) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;



if (wbin==5) SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.32;
if (wbin==5) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.55;
if (wbin==5) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.62;
if (wbin==5) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.1;
if (wbin==5) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.4;
if (wbin==5) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==5) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.75;
if (wbin==5) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.85;
if (wbin==5) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*1.98;
if (wbin==5) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*2.;
if (wbin==5) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==5) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==5) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;


if (wbin==6) SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.32;
if (wbin==6) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.55;
if (wbin==6) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.62;
if (wbin==6) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.05;
if (wbin==6) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==6) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.65;
if (wbin==6) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.65;
if (wbin==6) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.85;
if (wbin==6) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*1.9;
if (wbin==6) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*1.9;
if (wbin==6) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==6) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==6) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;



if (wbin==7)  SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.32;
if (wbin==7) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.55;
if (wbin==7) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.62;
if (wbin==7) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.05;
if (wbin==7) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==7) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.65;
if (wbin==7) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.65;
if (wbin==7) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.87;
if (wbin==7) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*1.95;
if (wbin==7) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*1.9;
if (wbin==7) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==7) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.75;
if (wbin==7) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;

if (wbin==8)  SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.32;
if (wbin==8) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.55;
if (wbin==8) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.62;
if (wbin==8) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.05;
if (wbin==8) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==8) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.65;
if (wbin==8) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.65;
if (wbin==8) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.87;
if (wbin==8) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*1.95;
if (wbin==8) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*1.8;
if (wbin==8) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==8) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.75;
if (wbin==8) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;


if (wbin==9)  SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][1][is12-1][itheta-1][ialpha-1]*0.32;
if (wbin==9) SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][2][is12-1][itheta-1][ialpha-1]*0.55;
if (wbin==9) SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][3][is12-1][itheta-1][ialpha-1]*0.62;
if (wbin==9) SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][4][is12-1][itheta-1][ialpha-1]*1.05;
if (wbin==9) SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][5][is12-1][itheta-1][ialpha-1]*1.5;
if (wbin==9) SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][6][is12-1][itheta-1][ialpha-1]*1.65;
if (wbin==9) SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][7][is12-1][itheta-1][ialpha-1]*1.65;
if (wbin==9) SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][8][is12-1][itheta-1][ialpha-1]*1.87;
if (wbin==9) SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][9][is12-1][itheta-1][ialpha-1]*1.95;
if (wbin==9) SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][10][is12-1][itheta-1][ialpha-1]*1.8;
if (wbin==9) SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][11][is12-1][itheta-1][ialpha-1]*1.55;
if (wbin==9) SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][12][is12-1][itheta-1][ialpha-1]*1.75;
if (wbin==9) SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][13][is12-1][itheta-1][ialpha-1]*1.4;
};
};
};




for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {

if (wbin==0) SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1]*0.2;

if (wbin==1) SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1]*0.55;
if (wbin==1) SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1]*0.2;
if (wbin==1) SIGMA_ARR_GOL2[wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][4][itheta-1][ialpha-1]*0.75;

if (wbin==2) SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1]*0.75;
if (wbin==2) SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1]*0.55;


if (wbin==3) SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1]*0.4;
if (wbin==3) SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1]*0.4;
if (wbin==3) SIGMA_ARR_GOL2[wbin][is23-1][4][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][4][itheta-1][ialpha-1]*0.75;


if (wbin==4) SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1]*0.4;

if (wbin==5) SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1]*0.4;
if (wbin==6) SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1]*0.4;

if (wbin==7) SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1]*0.4;


if (wbin==8) SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1]*0.4;
if (wbin==8) SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][3][itheta-1][ialpha-1]*0.94;
if (wbin==8) SIGMA_ARR_GOL2[wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][9][itheta-1][ialpha-1]*0.9;
if (wbin==8) SIGMA_ARR_GOL2[wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][12][itheta-1][ialpha-1]*0.9;


if (wbin==9) SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][2][itheta-1][ialpha-1]*0.5;
if (wbin==9) SIGMA_ARR_GOL2[wbin][is23-1][9][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][9][itheta-1][ialpha-1]*0.75;
if (wbin==9) SIGMA_ARR_GOL2[wbin][is23-1][10][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][10][itheta-1][ialpha-1]*0.9;
if (wbin==9) SIGMA_ARR_GOL2[wbin][is23-1][11][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][11][itheta-1][ialpha-1]*0.9;

if (wbin==9) SIGMA_ARR_GOL2[wbin][is23-1][12][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][12][itheta-1][ialpha-1]*0.75;
};
};
};


for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {

SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0]*5.5;
SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5]*5.5;

SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][1]*1.5;
SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][4]*1.5;



SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][2] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][2]*0.75;
SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][3] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][3]*0.75;


//if (wbin==1) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][1]*1.05;
//if (wbin==1) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][4] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][4]*1.05;

if (wbin==6) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0]*0.8;
if (wbin==6) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5]*0.8;


if (wbin==7) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0]*0.85;
if (wbin==7) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5]*0.85;

if (wbin==8) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0]*0.85;
if (wbin==8) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5]*0.85;

if (wbin==9) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][0]*0.85;
if (wbin==9) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][itheta-1][5]*0.85;

};
};
};




for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {

for (Int_t ialpha = 1; ialpha <=6; ialpha++) {

if (wbin==0) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][4][ialpha-1]*0.8;
if (wbin==0) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][5][ialpha-1]*1.2;

if (wbin==2) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][4][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][4][ialpha-1]*0.9;
if (wbin==2) SIGMA_ARR_GOL2[wbin][is23-1][is12-1][5][ialpha-1] = SIGMA_ARR_GOL2[wbin][is23-1][is12-1][5][ialpha-1]*1.1;


};
};
};


};//end lop rip3 - i



Float_t a[6][15] = {{0.7, 0.6,0.55,0.54,0.53,0.52,0.51,0.5,0.4,0.3,0.3,0.3,0.3,0.2,0.2},//sigma_t
{5., 4.8, 4.6, 4.4 ,4.2, 4.,3.75,3.5,3.25,3.,2.9,2.8,2.7,2.6,2.5},//sigma_l
{0.07, 0.06,0.055,0.054,0.053,0.05,0.045,0.04,0.035,0.03,0.03,0.03,0.03,0.02,0.02},//sigma_c2f
{0.04, 0.035,0.03,0.025,0.02,0.012,0.008,0.007,0.006,0.005,0.003,0.003,0.003,0.002,0.002},//sigma_s2f
{0.3, 0.27,0.25,0.2,0.15,0.13,0.1,0.07,0.05,0.03,0.03,0.03,0.025,0.02,0.015},//sigma_cf
{0.07, 0.06,0.055,0.054,0.053,0.052,0.051,0.045,0.03,0.025,0.021,0.017,0.013,0.012,0.01}};//sigma_sf





Float_t b[15] = {7, 6, 5.5, 5.4, 5.3, 2.5, 2., 1.7, 1.55, 1.45, 1.3, 1.2, 1.1, 1.05, 1.05};

for (Short_t i=0; i<=14; i++) {
wbin = i;

 W_ARR_gt_3[i] = 3.1375 + 0.1*i;
 
// cout << i<< " "<< W_ARR_gt_3[i]<< "\n";
 
for (Short_t is12=0; is12<=11; is12++) {
S12_ARR_gt_3[is12][wbin] = S12_ARR_RIP3[is12][9];
//cout << is12 <<  " "<< S12_ARR_gt_3[is12][wbin]<<"\n";
};

S12_ARR_gt_3[15][wbin] =  (W_ARR_gt_3[wbin] - MP)*(W_ARR_gt_3[wbin] - MP);

//S12_ARR_gt_3[6][wbin] =  S12_ARR_gt_3[5][wbin]+ (S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[5][wbin])/10.;
//S12_ARR_gt_3[7][wbin] =  S12_ARR_gt_3[6][wbin]+ (S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[5][wbin])/10.;
//S12_ARR_gt_3[8][wbin] =  S12_ARR_gt_3[7][wbin]+ (S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[5][wbin])/10.;
//S12_ARR_gt_3[9][wbin] =  S12_ARR_gt_3[8][wbin]+ (S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[5][wbin])/10.;
//S12_ARR_gt_3[10][wbin] =  S12_ARR_gt_3[9][wbin]+ (S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[5][wbin])/10.;
//S12_ARR_gt_3[11][wbin] =  S12_ARR_gt_3[10][wbin]+ (S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[5][wbin])/10.;
S12_ARR_gt_3[12][wbin] =  S12_ARR_gt_3[11][wbin]+ (S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[11][wbin])/4.;
S12_ARR_gt_3[13][wbin] =  S12_ARR_gt_3[12][wbin]+ (S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[11][wbin])/4.;
S12_ARR_gt_3[14][wbin] =  S12_ARR_gt_3[13][wbin]+ (S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[11][wbin])/4.;

//cout << 14<< " "<< S12_ARR_gt_3[14][wbin]<<"\n";
//cout << 15<< " "<< S12_ARR_gt_3[15][wbin]<<"\n";

//-----------------------------------------------------
for (Short_t is23=0; is23<=11; is23++) {
S23_ARR_gt_3[is23][wbin] = S23_ARR_RIP3[is23][9];
};

S23_ARR_gt_3[15][wbin] =  (W_ARR_gt_3[wbin] - MPIM)*(W_ARR_gt_3[wbin] - MPIM);

//S23_ARR_gt_3[6][wbin] = S23_ARR_gt_3[5][wbin] + (S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[5][wbin])/10.;
//S23_ARR_gt_3[7][wbin] = S23_ARR_gt_3[6][wbin] + (S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[5][wbin])/10.;
//S23_ARR_gt_3[8][wbin] = S23_ARR_gt_3[7][wbin] + (S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[5][wbin])/10.;
//S23_ARR_gt_3[9][wbin] = S23_ARR_gt_3[8][wbin] + (S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[5][wbin])/10.;
//S23_ARR_gt_3[10][wbin] = S23_ARR_gt_3[9][wbin] + (S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[5][wbin])/10.;
//S23_ARR_gt_3[11][wbin] = S23_ARR_gt_3[10][wbin] + (S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[5][wbin])/10.;
S23_ARR_gt_3[12][wbin] = S23_ARR_gt_3[11][wbin] + (S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[11][wbin])/4.;
S23_ARR_gt_3[13][wbin] = S23_ARR_gt_3[12][wbin] + (S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[11][wbin])/4.;
S23_ARR_gt_3[14][wbin] = S23_ARR_gt_3[13][wbin] + (S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[11][wbin])/4.;

for (Int_t j = 0; j <=5; j++) {

for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {

SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 1.5*SIGMA_ARR_RIP3[j][9][is23-1][is12-1][itheta-1][ialpha-1]/(S12_ARR_RIP3[15][9]-S12_ARR_RIP3[0][9])/(S23_ARR_RIP3[15][9]-S23_ARR_RIP3[0][9])*(S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[0][wbin])*(S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[0][wbin]);


if (!(SIGMA_ARR_RIP3[j][9][is23-1][is12-1][itheta-1][ialpha-1] > 0.)&&!(SIGMA_ARR_RIP3[j][9][is23-1][is12-1][itheta-1][ialpha-1] < 0.)) {

SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 1.5*a[j][wbin]/(S12_ARR_RIP3[15][9]-S12_ARR_RIP3[0][9])/(S23_ARR_RIP3[15][9]-S23_ARR_RIP3[0][9])*(S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[0][wbin])*(S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[0][wbin]);
};

if (j==0)  SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.9*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];

//if ((j==2))  SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 15*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];

if (wbin==0) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.94*0.7*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==1) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.94*0.65*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==2) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.93*0.63*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==3) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.87*0.6*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==4) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.80*0.58*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==5) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.80*0.57*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==6) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.80*1.05*1.05*0.6*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==7) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.75*1.05*1.05*0.6*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==8) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.72*1.05*1.05*0.59*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==9) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.68*1.05*1.05*0.58*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==10) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.65*1.05*1.05*0.57*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==11) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.63*1.2*0.56*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==12) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.61*1.2*0.54*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==13) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.61*1.2*0.52*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==14) SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.60*1.2*0.5*SIGMA_ARR_gt_3[j][wbin][is23-1][is12-1][itheta-1][ialpha-1];




SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = SIGMA_ARR_GOL2[9][is23-1][is12-1][itheta-1][ialpha-1]/(S12_ARR_RIP3[15][9]-S12_ARR_RIP3[0][9])/(S23_ARR_RIP3[15][9]-S23_ARR_RIP3[0][9])*(S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[0][wbin])*(S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[0][wbin]);

if (!(SIGMA_ARR_GOL2[9][is23-1][is12-1][itheta-1][ialpha-1] > 0.)&&!(SIGMA_ARR_GOL2[9][is23-1][is12-1][itheta-1][ialpha-1] < 0.)){

SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = b[wbin]/(S12_ARR_RIP3[15][9]-S12_ARR_RIP3[0][9])/(S23_ARR_RIP3[15][9]-S23_ARR_RIP3[0][9])*(S12_ARR_gt_3[15][wbin]-S12_ARR_gt_3[0][wbin])*(S23_ARR_gt_3[15][wbin]-S23_ARR_gt_3[0][wbin]);
};


if (wbin==0) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.95*0.6*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==1) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.91*0.56*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==2) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.88*0.52*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==3) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.86*0.48*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==4) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.84*0.46*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==5) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.88*0.44*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==6) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.86*0.47*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==7) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.85*0.46*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==8) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.85*0.43*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==9) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.85*0.41*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==10) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.85*0.37*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==11) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.85*0.37*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==12) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.85*0.34*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==13) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.84*0.32*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];
if (wbin==14) SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1] = 0.83*0.31*SIGMA_ARR_phot_gt_3[wbin][is23-1][is12-1][itheta-1][ialpha-1];



//cout << SIGMA_ARR_gt_3[0][wbin][is23-1][is12-1][itheta-1][ialpha-1]<< "\n"; 



};
};
};
};
};

for (Int_t j = 0; j <=5; j++) {
for (Int_t is23 = 1; is23 <=16; is23++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {
SIGMA_ARR_gt_3[j][wbin][is23-1][3][itheta-1][ialpha-1] = 1.03*SIGMA_ARR_gt_3[j][wbin][is23-1][3][itheta-1][ialpha-1];
//SIGMA_ARR_gt_3[j][wbin][is23-1][4][itheta-1][ialpha-1] = 0.98*SIGMA_ARR_gt_3[j][wbin][is23-1][4][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][is23-1][5][itheta-1][ialpha-1] = 1.03*SIGMA_ARR_gt_3[j][wbin][is23-1][5][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][is23-1][6][itheta-1][ialpha-1] = 0.96*SIGMA_ARR_gt_3[j][wbin][is23-1][6][itheta-1][ialpha-1];

SIGMA_ARR_gt_3[j][wbin][is23-1][9][itheta-1][ialpha-1] = 0.9*SIGMA_ARR_gt_3[j][wbin][is23-1][9][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][is23-1][10][itheta-1][ialpha-1] = 0.9*SIGMA_ARR_gt_3[j][wbin][is23-1][10][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][is23-1][11][itheta-1][ialpha-1] = 0.88*SIGMA_ARR_gt_3[j][wbin][is23-1][11][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][is23-1][12][itheta-1][ialpha-1] = 0.85*SIGMA_ARR_gt_3[j][wbin][is23-1][12][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][is23-1][13][itheta-1][ialpha-1] = 0.88*SIGMA_ARR_gt_3[j][wbin][is23-1][13][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][is23-1][14][itheta-1][ialpha-1] = 0.6*SIGMA_ARR_gt_3[j][wbin][is23-1][14][itheta-1][ialpha-1];


if (wbin >5) SIGMA_ARR_gt_3[j][wbin][is23-1][2][itheta-1][ialpha-1] = 0.55*SIGMA_ARR_gt_3[j][wbin][is23-1][2][itheta-1][ialpha-1];

//---------------------------------------------------------------
//SIGMA_ARR_phot_gt_3[wbin][is23-1][2][itheta-1][ialpha-1] = 0.75*SIGMA_ARR_phot_gt_3[wbin][is23-1][2][itheta-1][ialpha-1];

if (wbin > 5) SIGMA_ARR_phot_gt_3[wbin][is23-1][2][itheta-1][ialpha-1] = 0.87*SIGMA_ARR_phot_gt_3[wbin][is23-1][2][itheta-1][ialpha-1];
//if (wbin > 10) SIGMA_ARR_phot_gt_3[wbin][is23-1][2][itheta-1][ialpha-1] = 0.67*SIGMA_ARR_phot_gt_3[wbin][is23-1][2][itheta-1][ialpha-1];

//SIGMA_ARR_phot_gt_3[wbin][is23-1][9][itheta-1][ialpha-1] = 0.95*SIGMA_ARR_phot_gt_3[wbin][is23-1][9][itheta-1][ialpha-1];

//if (wbin <= 5) SIGMA_ARR_phot_gt_3[wbin][is23-1][6][itheta-1][ialpha-1] = 0.85*SIGMA_ARR_phot_gt_3[wbin][is23-1][6][itheta-1][ialpha-1];
//if (wbin <= 5) SIGMA_ARR_phot_gt_3[wbin][is23-1][7][itheta-1][ialpha-1] = 0.8*SIGMA_ARR_phot_gt_3[wbin][is23-1][7][itheta-1][ialpha-1];

//  SIGMA_ARR_phot_gt_3[wbin][is23-1][5][itheta-1][ialpha-1] = 1.05*SIGMA_ARR_phot_gt_3[wbin][is23-1][5][itheta-1][ialpha-1];
  SIGMA_ARR_phot_gt_3[wbin][is23-1][6][itheta-1][ialpha-1] = 0.95*SIGMA_ARR_phot_gt_3[wbin][is23-1][6][itheta-1][ialpha-1];
  SIGMA_ARR_phot_gt_3[wbin][is23-1][7][itheta-1][ialpha-1] = 0.93*SIGMA_ARR_phot_gt_3[wbin][is23-1][7][itheta-1][ialpha-1];
  SIGMA_ARR_phot_gt_3[wbin][is23-1][8][itheta-1][ialpha-1] = 0.98*SIGMA_ARR_phot_gt_3[wbin][is23-1][8][itheta-1][ialpha-1];
  SIGMA_ARR_phot_gt_3[wbin][is23-1][9][itheta-1][ialpha-1] = 0.95*SIGMA_ARR_phot_gt_3[wbin][is23-1][9][itheta-1][ialpha-1];
  SIGMA_ARR_phot_gt_3[wbin][is23-1][10][itheta-1][ialpha-1] = 0.95*SIGMA_ARR_phot_gt_3[wbin][is23-1][10][itheta-1][ialpha-1];
 //  SIGMA_ARR_phot_gt_3[wbin][is23-1][11][itheta-1][ialpha-1] = 1.03*SIGMA_ARR_phot_gt_3[wbin][is23-1][11][itheta-1][ialpha-1];
  
};
};
};

for (Int_t is12 = 1; is12 <=16; is12++) {
for (Int_t itheta = 1; itheta <=6; itheta++) {
for (Int_t ialpha = 1; ialpha <=6; ialpha++) {

SIGMA_ARR_gt_3[j][wbin][3][is12-1][itheta-1][ialpha-1] = 1.21*SIGMA_ARR_gt_3[j][wbin][3][is12-1][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][4][is12-1][itheta-1][ialpha-1] = 1.12*SIGMA_ARR_gt_3[j][wbin][4][is12-1][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][5][is12-1][itheta-1][ialpha-1] = 0.92*SIGMA_ARR_gt_3[j][wbin][5][is12-1][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][6][is12-1][itheta-1][ialpha-1] = 0.88*SIGMA_ARR_gt_3[j][wbin][6][is12-1][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][7][is12-1][itheta-1][ialpha-1] = 0.95*SIGMA_ARR_gt_3[j][wbin][7][is12-1][itheta-1][ialpha-1];
SIGMA_ARR_gt_3[j][wbin][11][is12-1][itheta-1][ialpha-1] = 0.9*SIGMA_ARR_gt_3[j][wbin][11][is12-1][itheta-1][ialpha-1];

if (wbin > 5) SIGMA_ARR_gt_3[j][wbin][1][is12-1][itheta-1][ialpha-1] = 0.46*SIGMA_ARR_gt_3[j][wbin][1][is12-1][itheta-1][ialpha-1];
if (wbin > 5) SIGMA_ARR_gt_3[j][wbin][2][is12-1][itheta-1][ialpha-1] = 0.85*SIGMA_ARR_gt_3[j][wbin][2][is12-1][itheta-1][ialpha-1];
if (wbin > 5) SIGMA_ARR_gt_3[j][wbin][5][is12-1][itheta-1][ialpha-1] = 0.7*SIGMA_ARR_gt_3[j][wbin][5][is12-1][itheta-1][ialpha-1];
 SIGMA_ARR_gt_3[j][wbin][6][is12-1][itheta-1][ialpha-1] = 0.72*SIGMA_ARR_gt_3[j][wbin][6][is12-1][itheta-1][ialpha-1];
 SIGMA_ARR_gt_3[j][wbin][7][is12-1][itheta-1][ialpha-1] = 0.85*SIGMA_ARR_gt_3[j][wbin][7][is12-1][itheta-1][ialpha-1];

if (wbin > 10) SIGMA_ARR_gt_3[j][wbin][1][is12-1][itheta-1][ialpha-1] = 0.3*SIGMA_ARR_gt_3[j][wbin][1][is12-1][itheta-1][ialpha-1];
if (wbin > 10) SIGMA_ARR_gt_3[j][wbin][2][is12-1][itheta-1][ialpha-1] = 0.6*SIGMA_ARR_gt_3[j][wbin][2][is12-1][itheta-1][ialpha-1];

//--------------------------------------------

SIGMA_ARR_phot_gt_3[wbin][11][is12-1][itheta-1][ialpha-1] = 1.15*SIGMA_ARR_phot_gt_3[wbin][11][is12-1][itheta-1][ialpha-1] ;
SIGMA_ARR_phot_gt_3[wbin][12][is12-1][itheta-1][ialpha-1] = 1.2*SIGMA_ARR_phot_gt_3[wbin][12][is12-1][itheta-1][ialpha-1] ;
SIGMA_ARR_phot_gt_3[wbin][13][is12-1][itheta-1][ialpha-1] = 1.2*SIGMA_ARR_phot_gt_3[wbin][13][is12-1][itheta-1][ialpha-1] ;
SIGMA_ARR_phot_gt_3[wbin][14][is12-1][itheta-1][ialpha-1] = 1.2*SIGMA_ARR_phot_gt_3[wbin][14][is12-1][itheta-1][ialpha-1] ;


SIGMA_ARR_phot_gt_3[wbin][5][is12-1][itheta-1][ialpha-1] = 0.99*SIGMA_ARR_phot_gt_3[wbin][5][is12-1][itheta-1][ialpha-1] ;
//SIGMA_ARR_phot_gt_3[wbin][6][is12-1][itheta-1][ialpha-1] = 0.95*SIGMA_ARR_phot_gt_3[wbin][6][is12-1][itheta-1][ialpha-1] ;
SIGMA_ARR_phot_gt_3[wbin][7][is12-1][itheta-1][ialpha-1] = 1.05*SIGMA_ARR_phot_gt_3[wbin][7][is12-1][itheta-1][ialpha-1] ;
SIGMA_ARR_phot_gt_3[wbin][8][is12-1][itheta-1][ialpha-1] = 1.1*SIGMA_ARR_phot_gt_3[wbin][8][is12-1][itheta-1][ialpha-1] ;
SIGMA_ARR_phot_gt_3[wbin][9][is12-1][itheta-1][ialpha-1] = 1.15*SIGMA_ARR_phot_gt_3[wbin][9][is12-1][itheta-1][ialpha-1] ;
SIGMA_ARR_phot_gt_3[wbin][10][is12-1][itheta-1][ialpha-1] = 1.11*SIGMA_ARR_phot_gt_3[wbin][10][is12-1][itheta-1][ialpha-1] ;


};
};
};



};







};//end lop w.gt.3 - i


 
};

