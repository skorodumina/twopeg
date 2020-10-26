#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
#include "global.h"
#include <sstream>

using namespace std;



void read_fit_param_files(){

string file_names[10];

Short_t wbin, q2bin;

//files' names
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/sigma_t_fit_param.dat";
file_names[0] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/sigma_l_fit_param.dat";
file_names[1] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/sigma_c2f_fit_param.dat";
file_names[2] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/sigma_cf_fit_param.dat";
file_names[3] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/sigmas_int_rip_065_new.dat";
file_names[4] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/golova_photon.dat";
file_names[5] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/sigma_t_fit_param_fed_lin.dat";
file_names[6] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/sigma_l_fit_param_fed.dat";
file_names[7] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/sigma_c2f_fit_param_fed.dat";
file_names[8] = PATH.str();
PATH.str("");
PATH << data_dir_2pi.str() << "data/fit_param/sigma_cf_fit_param_fed.dat";
file_names[9] = PATH.str();


string dummy,xsect;


string file=file_names[0];
ifstream input_t(file.c_str());
if (input_t.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_t.is_open()){

for (Int_t iw = 0; iw <=16; iw++) {
getline(input_t,dummy);
getline(input_t,xsect);
FIT_PARAM_SIGMA_T[0][iw] = atof(xsect.c_str());
getline(input_t,xsect);
FIT_PARAM_SIGMA_T[1][iw] = atof(xsect.c_str());
getline(input_t,xsect);
FIT_PARAM_SIGMA_T[2][iw] = atof(xsect.c_str());
getline(input_t,xsect);
FIT_PARAM_SIGMA_T[3][iw] = atof(xsect.c_str());
};
};
input_t.close();


file=file_names[1];
ifstream input_l(file.c_str());
if (input_l.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_l.is_open()){

for (Int_t iw = 0; iw <=16; iw++) {
getline(input_l,dummy);
getline(input_l,xsect);
FIT_PARAM_SIGMA_L[0][iw] = atof(xsect.c_str());
getline(input_l,xsect);
FIT_PARAM_SIGMA_L[1][iw] = atof(xsect.c_str());
getline(input_l,xsect);
FIT_PARAM_SIGMA_L[2][iw] = atof(xsect.c_str());
};
};
input_l.close();


file=file_names[2];
ifstream input_c2f(file.c_str());
if (input_c2f.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_c2f.is_open()){

for (Int_t iw = 0; iw <=16; iw++) {
getline(input_c2f,dummy);
getline(input_c2f,xsect);
FIT_PARAM_SIGMA_C2F[0][iw] = atof(xsect.c_str());
getline(input_c2f,xsect);
FIT_PARAM_SIGMA_C2F[1][iw] = atof(xsect.c_str());
getline(input_c2f,xsect);
FIT_PARAM_SIGMA_C2F[2][iw] = atof(xsect.c_str());
};
};
input_c2f.close();


file=file_names[3];
ifstream input_cf(file.c_str());
if (input_cf.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_cf.is_open()){

for (Int_t iw = 0; iw <=16; iw++) {
getline(input_cf,dummy);
getline(input_cf,xsect);
FIT_PARAM_SIGMA_CF[0][iw] = atof(xsect.c_str());
getline(input_cf,xsect);
FIT_PARAM_SIGMA_CF[1][iw] = atof(xsect.c_str());
getline(input_cf,xsect);
FIT_PARAM_SIGMA_CF[2][iw] = atof(xsect.c_str());
};
};
input_cf.close();


file=file_names[4];
ifstream input_int_rip065(file.c_str());
if (input_int_rip065.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_int_rip065.is_open()){

for (Int_t iw = 0; iw <=16; iw++) {
getline(input_int_rip065,dummy);
getline(input_int_rip065,dummy);
getline(input_int_rip065,dummy);
getline(input_int_rip065,xsect);
SIGMA_T_INT_RIPANI[iw] = atof(xsect.c_str());
getline(input_int_rip065,dummy);
getline(input_int_rip065,dummy);
getline(input_int_rip065,dummy);
getline(input_int_rip065,dummy);
getline(input_int_rip065,dummy);
};
};
input_int_rip065.close();

file=file_names[5];
ifstream input_int_gol0(file.c_str());
if (input_int_gol0.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_int_gol0.is_open()){

for (Int_t iw = 0; iw <=16; iw++) {
getline(input_int_gol0,dummy);
getline(input_int_gol0,xsect);
SIGMA_T_INT_GOLOVA[iw] = atof(xsect.c_str());
};
};
input_int_gol0.close();



file=file_names[6];
ifstream input_fed_t(file.c_str());
if (input_fed_t.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_fed_t.is_open()){

for (Int_t iw = 0; iw <=9; iw++) {
getline(input_fed_t,dummy);
getline(input_fed_t,xsect);
FIT_PARAM_SIGMA_T_FED[0][iw] = atof(xsect.c_str());
getline(input_fed_t,xsect);
FIT_PARAM_SIGMA_T_FED[1][iw] = atof(xsect.c_str());
};
};
input_fed_t.close();


file=file_names[7];
ifstream input_fed_l(file.c_str());
if (input_fed_l.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_fed_l.is_open()){

for (Int_t iw = 0; iw <=9; iw++) {
getline(input_fed_l,dummy);
getline(input_fed_l,xsect);
FIT_PARAM_SIGMA_L_FED[0][iw] = atof(xsect.c_str());
getline(input_fed_l,xsect);
FIT_PARAM_SIGMA_L_FED[1][iw] = atof(xsect.c_str());
getline(input_fed_l,xsect);
FIT_PARAM_SIGMA_L_FED[2][iw] = atof(xsect.c_str());


};
};
input_fed_l.close();



file=file_names[8];
ifstream input_fed_c2f(file.c_str());
if (input_fed_c2f.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_fed_c2f.is_open()){

for (Int_t iw = 0; iw <=9; iw++) {
getline(input_fed_c2f,dummy);
getline(input_fed_c2f,xsect);
FIT_PARAM_SIGMA_C2F_FED[0][iw] = atof(xsect.c_str());
getline(input_fed_c2f,xsect);
FIT_PARAM_SIGMA_C2F_FED[1][iw] = atof(xsect.c_str());
getline(input_fed_c2f,xsect);
FIT_PARAM_SIGMA_C2F_FED[2][iw] = atof(xsect.c_str());


};
};
input_fed_c2f.close();


file=file_names[9];
ifstream input_fed_cf(file.c_str());
if (input_fed_cf.fail()) cout<< "ALARM! File with fit parameters FAILED to open! \n";
if(input_fed_cf.is_open()){

for (Int_t iw = 0; iw <=9; iw++) {
getline(input_fed_cf,dummy);
getline(input_fed_cf,xsect);
FIT_PARAM_SIGMA_CF_FED[0][iw] = atof(xsect.c_str());
getline(input_fed_cf,xsect);
FIT_PARAM_SIGMA_CF_FED[1][iw] = atof(xsect.c_str());
getline(input_fed_cf,xsect);
FIT_PARAM_SIGMA_CF_FED[2][iw] = atof(xsect.c_str());

};
};
input_fed_cf.close();






 return;
};


