#include "TFile.h"
#include "TMath.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include <sstream>
 using namespace std;

//This file containts subroutines related to taking the values of the input paramenters.
//input_stream(...) -- takes the input parameters from the cin stream, which flows to the executable from the input file. 
//input_cmd_line(...) -- for supporting the command line input. All the parameters are first taken from "data/inp_cmd_line". Then Nevents can be reassigned according to the command line "--trig" option if the one was specified.
//inp_couts(...) -- for interactive communication with the user. Prints to the screen all the input parameters along with some comments on their compatibility.

void input_stream(Float_t &E_beam) {

cout << endl;
cout << "The cin stream input is used\n";
cout << endl;
 
string qqq;
 
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));    
Nevents = atoi(qqq.c_str());
      
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
E_beam = atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
W_min = atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
W_max = atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Q2_min = atof(qqq.c_str());
   
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Q2_max = atof(qqq.c_str());
   
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Theta_min = atof(qqq.c_str());
   
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Theta_max = atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
E_eprime_min = atof(qqq.c_str());
   
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_rad = atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_len = atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_off = atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_dens= atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_radlen= atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_Z= atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_A= atof(qqq.c_str());
    
getline (cin,qqq);
Twi_thick= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_thick= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 

getline (cin,qqq);
Twi_dens= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_dens= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 

getline (cin,qqq);
Twi_radlen= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_radlen= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 

getline (cin,qqq);
Twi_Z= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_Z= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 

getline (cin,qqq);
Twi_A= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_A= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 
   
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_bos= atof(qqq.c_str());
    
getline (cin,qqq);
out_bos_file = qqq.substr(0, qqq.find(" ",0));
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_lund= atof(qqq.c_str());
 
getline (cin,qqq);
out_lund_file = qqq.substr(0, qqq.find(" ",0));
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_radmod= atof(qqq.c_str());
   
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_fermi= atof(qqq.c_str());
    
getline (cin,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_flux= atof(qqq.c_str());

};

//------------------------

void input_cmd_line(Float_t &E_beam, Int_t &argc,  char *argv[], Short_t &flag_seed, Int_t &seed) {

cout << endl;
cout << "The cmd line input is used \nSee 'data/inp_cmd_line' for input parameters\n";
cout << endl; 

string qqq;
string inp_file_name;
bool check_open_fail;
Int_t n_trig = 0;

//Building the full name of the input file 
PATH << data_dir_2pi.str() << "data/inp_cmd_line";
inp_file_name = PATH.str();
PATH.str("");

string dummy,xsect;

string file=inp_file_name;
ifstream input(file.c_str());

//Chencking if the file managed to open
check_open_fail = input.fail();

if (check_open_fail) cout << "ALARM! 'data/inp_cmd_line' FAILED to open!!! \n\n"; 
if(input.is_open()){ 

getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));    
Nevents = atoi(qqq.c_str());
      
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
E_beam = atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
W_min = atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
W_max = atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Q2_min = atof(qqq.c_str());
   
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Q2_max = atof(qqq.c_str());
   
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Theta_min = atof(qqq.c_str());
   
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Theta_max = atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
E_eprime_min = atof(qqq.c_str());
   
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_rad = atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_len = atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_off = atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_dens= atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_radlen= atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_Z= atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
Targ_A= atof(qqq.c_str());
    
getline (input,qqq);
Twi_thick= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_thick= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 

getline (input,qqq);
Twi_dens= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_dens= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 

getline (input,qqq);
Twi_radlen= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_radlen= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 

getline (input,qqq);
Twi_Z= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_Z= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 

getline (input,qqq);
Twi_A= atof(qqq.substr(0, qqq.find(",",0)).c_str());
Twf_A= atof(qqq.substr(qqq.substr(0, qqq.find(",",0)).length()+1,qqq.rfind(",",0)).c_str()); 
   
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_bos= atof(qqq.c_str());
    
getline (input,qqq);
out_bos_file = qqq.substr(0, qqq.find(" ",0));
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_lund= atof(qqq.c_str());
 
getline (input,qqq);
out_lund_file = qqq.substr(0, qqq.find(" ",0));
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_radmod= atof(qqq.c_str());
   
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_fermi= atof(qqq.c_str());
    
getline (input,qqq);
qqq = qqq.substr(0, qqq.find(" ",0));
flag_flux= atof(qqq.c_str());

};
input.close();

//Reassigment of the input parameters if cmd line value is different from that from the data/inp_cmd_line
//Only number of event can be reassigned via the cmd line argument.
//Other input parameters are taken from  data/inp_cmd_line
if (!check_open_fail){
for (Short_t i=0;i<argc;i++){
std::string arg = argv[i];

if (arg=="--trig") {
n_trig = abs(atoi(argv[i+1]));
if (!(n_trig == 0)){
Nevents = n_trig;
cout << "Accoring to '--trig', Nevents was set to " <<Nevents<<"\n";
};
};

if (arg=="--seed") {
seed = abs(atoi(argv[i+1]));
if (!(seed == 0)){
flag_seed = 1;
cout << "Accoring to '--seed', seed was set to " <<seed<<"\n";
};
};

};
};
cout << endl;

};     
    
//-----------------------------------
   
  
void inp_couts(Float_t &E_beam) {
 
  
cout <<"____________INPUT PARAMETERS:______________\n\n";

cout << "Number of events to be generated is " << Nevents << "\n";
cout << "Beam energy is " << E_beam << " GeV" << "\n";
cout << "W_min is  " << W_min << " GeV" << "\n";
cout << "W_max is  " << W_max << " GeV" << "\n";
cout << "Q2_min is  " << Q2_min << " GeV^2" << "\n";
cout << "Q2_max is  " << Q2_max << " GeV^2" << "\n";
cout << "Theta_min is  " << Theta_min << " deg" << "\n";    
cout << "Theta_max is  " << Theta_max << " deg" << "\n";        
cout << "Minimal energy of scattered electron  " << E_eprime_min << " GeV" << "\n";  
cout << "Target radius is  " << Targ_rad << " cm" << "\n";
cout << "Target length is  " << Targ_len << " cm" << "\n";
cout << "Target offset in z is  " << Targ_off << " cm" << "\n";
cout << "Target density is  " << Targ_dens << " g/cm^3" << "\n";
cout << "Target radiation length is  " << Targ_radlen << " cm" << "\n";
cout << "Target Z is  " << Targ_Z <<  "\n";
cout << "Target A is  " << Targ_A <<  "\n";
printf ("Thickness of the target windows initial, final %f,%f um\n", Twi_thick,Twf_thick);
printf ("Density of the target windows initial, final %f,%f g/cm^3\n", Twi_dens,Twf_dens);
printf ("Radiation length of the target windows initial, final %f,%f cm\n", Twi_radlen,Twf_radlen);
cout << "Z of the target windows initial, final  " << Twi_Z<< ", " << Twf_Z <<  "\n";
cout << "A of the target windows initial, final  " << Twi_A<< ", " << Twf_A <<  "\n";
   
#ifdef BOS    
switch (flag_bos) {
case 0:  cout << "Ouput BOS flag  " << flag_bos << "  - no BOS output" << "\n";
break;
case 1:  cout << "Ouput BOS flag  " << flag_bos << "  -  output with MCTK, MCVX banks" << "\n";
break;
case 2:  cout << "Ouput BOS flag  " << flag_bos << "  -  output with PART bank" << "\n";
break;    
};    
#endif    
  
#ifndef BOS
if (!(flag_bos==0)){
flag_bos = 0;
cout << "Ouput BOS flag  " << flag_bos << "  - no BOS output as 'make nobos' was used (for BOS output compile 'make bos') " << "\n";
};
#endif     
  
    
#ifdef BOS      
if (flag_bos == 1) cout << "BOS output file name " << out_bos_file <<"\n";
#endif   
  
switch (flag_lund) {
case 0:  cout << "Ouput LUND flag  " << flag_lund << "  - no LUND output" << "\n";
break;
case 1:  cout << "Ouput LUND flag  " << flag_lund << "  -  output LUND file" << "\n";
break;
};
 
if (flag_lund == 1) cout << "LUND output file name " << out_lund_file <<"\n";
    
switch (flag_radmod) {
case 0:  cout << "Radiative mode flag  " << flag_radmod << "  - no rad effects" << "\n";
break;
case 1:  cout << "Radiative mode flag  " << flag_radmod << "  -  rad eff with no straggling" << "\n";
break;
case 2:  cout << "Radiative mode flag  " << flag_radmod << "  -  rad eff with straggling" << "\n"; 
break;
};
   
switch (flag_fermi) {
case 0:  cout << "Fermi flag " << flag_fermi << "  - no Fermi smearing" << "\n";
break;
case 1:  cout << "Fermi flag " << flag_fermi << "  -  with Fermi smearing" << "\n";
break;
};
    
    
switch (flag_flux) {
case 0: cout << "Flux flag " << flag_flux << "  - virtual photoproduction (model cross section)" << "\n";
break;
case 1: cout << "Flux flag " << flag_flux << "  -  electroproduction (like exp data)" << "\n";
break;
};
    
if((!(flag_radmod == 0))&&(flag_fermi==1))   cout <<"\nCAUTION! Fermi mode flag = "<<  flag_fermi<<" and Rad mode flag = "<<flag_radmod<<".\nIn this EG these two modes are combined naively.\nIt is recommended to use separate TWOPEG-D version,\nwhere Fermi mode and Rad mode are combined properly.\n";
//See CLAS12-Note-2017-001 and CLAS12-Note-2017-014 for details.
    
cout <<"___________________________________________\n\n";
   
     };   
     
     
     
     
     
     
