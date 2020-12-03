#include "TFile.h"
#include "TMath.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "global.h"
#include <sstream>
 using namespace std;

//This file containts subroutines related to taking the values of the input paramenters.
//input_stream(...) -- takes the input parameters from the cin stream a certain input file. 
//input_cline(...) -- for supporting the command line input. All the parameters are first set to their default values. Then they can be reassigned according to the command line options.
//inp_couts(...) -- for interactive communication with the user. Prints to the screen all the input parameters along with some comments on their compatibility.
//cmdl_help() -- provides the user with the list of the avaliable command line options and the defauls values of the input parameters. Executes via "--help".

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

void input_cline(Float_t &E_beam, Int_t &argc, char *argv[], Short_t &flag_seed, Int_t &seed) {
cout << endl;
cout << "The cmd line input is used"<< endl;
cout << "Input parameters are set to their default values"<< endl;
cout << "Use './twopeg --help' for cmd line options"<< endl;
cout << endl;

Int_t n_trig = 0;
Float_t var_f = 0.;
Short_t var_sh = 0;
string var_str;
std::string arg;
char *end;
std::string end_str;
//THESE ARE INPUT PARAMETERS SET BY DEFAULT

   Nevents = 10000; 	// Number of events to be generated
   E_beam = 10.6;	// Beam energy (GeV)
   W_min = 1.4;		// W minimum (GeV) 
   W_max = 2.3;		// Wmaximum (GeV)
   Q2_min = 0.05; 	// Q2 minimum (GeV^2)
   Q2_max = 3.5; 	// Q2 maximum (GeV^2)
   Theta_min = 1.0; 	// minimal theta of scattered electron (deg)
   Theta_max = 50.0;	// maximal theta of scattered electron (deg)
   E_eprime_min = 0.1;  // minimal energy of scattered electron (GeV)
   Targ_rad = 0.6; 	// Target radius in cm
   Targ_len = 2.0; 	// Target length in cm
   Targ_off = 0.0;	// Target offset in z in cm
   Targ_dens = 0.071;	// Target density (g/cm^3)//not used
   Targ_radlen = 888.03; // Target radiation length (cm)
   Targ_Z = 1; 		// Target Z
   Targ_A = 1;	 	// Target A//not used
   Twi_thick = 15.0;	// Thickness of the target windows initial, final (um)
   Twf_thick = 15.0;	// Thickness of the target windows initial, final (um)
   Twi_dens = 2.699;	// Density of target windows initial,final (g/cm^3)//not used
   Twf_dens = 2.699;	// Density of target windows initial,final (g/cm^3)//not used
   Twi_radlen = 8.897;	// Radiation length of target windows initial,final (cm)
   Twf_radlen = 8.897;	// Radiation length of target windows initial,final (cm)
   Twi_Z = 13;		// Target windows Z initial, final
   Twf_Z = 13;		// Target windows Z initial, final
   Twi_A = 27; 		// Target windows A initial, final//not used
   Twf_A = 27; 		// Target windows A initial, final//not used
#ifdef BOS
   flag_bos = 1; 	// Output BOS file 0 - no, 1 - MCTK,MCVX banks, 2 - PART bank
#endif
#ifndef BOS   
   flag_bos = 0; 	// Output BOS file 0 - no, 1 - MCTK,MCVX banks, 2 - PART bank
#endif   
   out_bos_file = "out.bos"; 	// BOS output file name
   flag_lund = 1;	// Otput LUND file 0 - no, 1 - yes
   out_lund_file = "twopeg.dat"; 	// LUND output file name
   flag_radmod = 0;	// Radiative mode 0 - no rad effects, 1 - rad eff with no straggling,  2 - rad eff with straggling
   flag_fermi = 0; 	//Fermi smearing: 0 - no, 1 - yes
   flag_flux = 1;	//Multiplication by virtual photon flux: 0 - no (virtual photoproduction), 1 - yes (electroproduction)
      
//----General and kinematics-----   
for (Short_t i=0;i<argc;i++){
arg = argv[i];

if (arg=="--trig") {
n_trig = abs(atoi(argv[i+1]));
if (!(n_trig == 0)){
Nevents = n_trig;
cout << "Accoring to '--trig',  Nevents was reset to " <<Nevents<<"\n";
};
};
if (arg=="--seed") {
seed = abs(atoi(argv[i+1]));
if (!(seed == 0)){
flag_seed = 1;
cout << "Accoring to '--seed',  seed    was   set to " <<seed<<"\n";
};
};

if (arg=="--ebeam") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
E_beam = var_f;
cout << "Accoring to '--ebeam', Ebeam   was reset to " <<E_beam<<" GeV\n";
};
var_f = 0.;
};

if (arg=="--wmin") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
W_min = var_f;
cout << "Accoring to '--wmin',  Wmin    was reset to " <<W_min<<" GeV\n";
};
var_f = 0.;
};

if (arg=="--wmax") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
W_max = var_f;
cout << "Accoring to '--wmax',  Wmax    was reset to " <<W_max<<" GeV\n";
};
var_f = 0.;
};

if (arg=="--q2min") {
var_f = strtof(argv[i+1], &end);
end_str = end;
if (end_str.size()==0){
Q2_min = abs(var_f);
cout << "Accoring to '--q2min', Q2min   was reset to " <<Q2_min<<" GeV^2\n";
};
var_f = 0.;
end = "";
};

if (arg=="--q2max") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Q2_max = var_f;
cout << "Accoring to '--q2max', Q2max   was reset to " <<Q2_max<<" GeV^2\n";
};
var_f = 0.;
};

if (arg=="--thmin") {
var_f = strtof(argv[i+1], &end);
end_str = end;
if (end_str.size()==0){
Theta_min = abs(var_f);
cout << "Accoring to '--thmin', THmin   was reset to " <<Theta_min<<" deg\n";
};
var_f = 0.;
end = "";
};

if (arg=="--thmax") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Theta_max = var_f;
cout << "Accoring to '--thmax', THmax   was reset to " <<Theta_max<<" deg\n";
};
var_f = 0.;
};

if (arg=="--emin") {
var_f = strtof(argv[i+1], &end);
end_str = end;
if (end_str.size()==0){
E_eprime_min = abs(var_f);
cout << "Accoring to '--emin',  Emin    was reset to " <<E_eprime_min <<" GeV\n";
};
var_f = 0.;
end = "";
};
};

//--------------Target--------------------------------
if (argc > 3) cout << endl;  
for (Short_t i=0;i<argc;i++){
arg = argv[i];

if (arg=="--trad") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Targ_rad = var_f;
cout << "Accoring to '--trad',  Target radius    was reset to " <<Targ_rad <<" cm\n";
};
var_f = 0.;
};
if (arg=="--tlen") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Targ_len = var_f;
cout << "Accoring to '--tlen',  Target length    was reset to " <<Targ_len<<" cm\n";
};
var_f = 0.;
};

if (arg=="--toff") {
var_f = strtof(argv[i+1], &end);
end_str = end;
if (end_str.size()==0){
Targ_off = var_f;
cout << "Accoring to '--toff',  Target offset    was reset to " <<Targ_off<<" cm\n";
};
var_f = 0.;
end = "";
};

if (arg=="--tden") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Targ_dens = var_f;
cout << "Accoring to '--tden',  Target density   was reset to " <<Targ_dens<<" g/cm^3\n";
};
var_f = 0.;
};
if (arg=="--trdln") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Targ_radlen = var_f;
cout << "Accoring to '--trdln', Target radlen    was reset to " <<Targ_radlen<<" cm\n";
};
var_f = 0.;
};
if (arg=="--tz") {
var_sh = abs(atoi(argv[i+1]));
if (!(var_sh == 0)){
Targ_Z = var_sh;
cout << "Accoring to '--tz',    Target Z         was reset to " <<Targ_Z<<"\n";
};
var_sh = 0.;
};
if (arg=="--ta") {
var_sh = abs(atoi(argv[i+1]));
if (!(var_sh == 0)){
Targ_A = var_sh;
cout << "Accoring to '--ta',    Target A         was reset to " <<Targ_A<<"\n";
};
var_sh = 0.;
};
};

//--------Target windows------
if (argc > 3) cout << endl;  
for (Short_t i=0;i<argc;i++){
arg = argv[i];

if (arg=="--twlen") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Twi_thick = var_f;
Twf_thick = var_f;
cout << "Accoring to '--twlen',   Target windows thickness (each)    was reset to " <<Twi_thick<<" um\n";
};
var_f = 0.;
}; 
if (arg=="--twflen") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Twf_thick = var_f;
cout << "Accoring to '--twflen',  Target windows thickness (final)   was reset to " <<Twf_thick<<" um\n";
};
var_f = 0.;
}; 
if (arg=="--twden") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Twi_dens = var_f;
Twf_dens = var_f;
cout << "Accoring to '--twden',   Target windows density   (each)    was reset to " <<Twi_dens<<" g/cm^3\n";
};
var_f = 0.;
}; 
if (arg=="--twfden") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Twf_dens = var_f;
cout << "Accoring to '--twfden',  Target windows density   (final)   was reset to " <<Twf_dens<<" g/cm^3\n";
};
var_f = 0.;
};
if (arg=="--twrdln") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Twi_radlen = var_f;
Twf_radlen = var_f;
cout << "Accoring to '--twrdln',  Target windows radlenth  (each)    was reset to " <<Twi_radlen<<" cm\n";
};
var_f = 0.;
};
if (arg=="--twfrdln") {
var_f = abs(atof(argv[i+1]));
if (!(var_f == 0)){
Twf_radlen = var_f;
cout << "Accoring to '--twfrdln', Target windows radlenth  (final)   was reset to " <<Twf_radlen<<" cm\n";
};
var_f = 0.;
};
if (arg=="--twz") {
var_sh = abs(atoi(argv[i+1]));
if (!(var_sh == 0)){
Twi_Z = var_sh;
Twf_Z = var_sh;
cout << "Accoring to '--twz',     Target windows Z         (each)    was reset to " <<Twi_Z<<"\n";
};
var_sh = 0.;
};
if (arg=="--twfz") {
var_sh = abs(atoi(argv[i+1]));
if (!(var_sh == 0)){
Twf_Z = var_sh;
cout << "Accoring to '--twfz',    Target windows Z         (final)   was reset to " <<Twf_Z<<"\n";
};
var_sh = 0.;
};
if (arg=="--twa") {
var_sh = abs(atoi(argv[i+1]));
if (!(var_sh == 0)){
Twi_A = var_sh;
Twf_A = var_sh;
cout << "Accoring to '--twa',     Target windows A         (each)    was reset to " <<Twi_A<<"\n";
};
var_sh = 0.;
};
if (arg=="--twfa") {
var_sh = abs(atoi(argv[i+1]));
if (!(var_sh == 0)){
Twf_A = var_sh;
cout << "Accoring to '--twfa',    Target windows A         (final)   was reset to " <<Twf_A<<"\n";
};
var_sh = 0.;
};
};

//-----Other (flags and outputs)-------
if (argc > 3) cout << endl;  
for (Short_t i=0;i<argc;i++){
arg = argv[i];

#ifdef BOS
if (arg=="--flagbos") {
var_f = strtof(argv[i+1], &end);
end_str = end;
if ((end_str.size() == 0)&&(abs(var_f) < 3)){
flag_bos = abs(int(var_f));
cout << "Accoring to '--flagbos',    Flag bos         was reset to " <<flag_bos<<"\n";
};
var_f = 0.;
end = "";
}; 
if ((flag_bos == 1)||(flag_bos == 2)){
if (arg=="--bosname") {
var_str = argv[i+1];
out_bos_file = var_str;
cout << "Accoring to '--bosname',    BOS  output      was renamed as " <<out_bos_file<<"\n";
var_str = "";
};
};
#endif

if (arg=="--flaglund") {
var_f = strtof(argv[i+1], &end);
end_str = end;
if ((end_str.size() == 0)&&(abs(var_f) < 2)){
flag_lund = abs(int(var_f));
cout << "Accoring to '--flaglund',   Flag lund        was reset to " <<flag_lund<<"\n";
};
var_f = 0.;
end = "";
}; 

if (arg=="--lundname") {
var_str = argv[i+1];
if (!(flag_lund == 0)){
out_lund_file = var_str;
cout << "Accoring to '--lundname',   LUND output      was renamed as " <<out_lund_file<<"\n";
};
var_str = "";
};

if (arg=="--flagrad") {
var_f = strtof(argv[i+1], &end);
end_str = end;
if ((end_str.size()==0)&&(abs(var_f) < 3)){
flag_radmod = abs(int(var_f));
cout << "Accoring to '--flagrad',    Flag rad mod     was reset to " <<flag_radmod<<"\n";
};
var_f = 0.;
end = "";
}; 
if (arg=="--flagfermi") {
var_f = strtof(argv[i+1], &end);
end_str = end;
if ((end_str.size()==0)&&(abs(var_f) < 2)){
flag_fermi = abs(int(var_f));
cout << "Accoring to '--flagfermi',  Flag Fermi       was reset to " <<flag_fermi<<"\n";
};
var_f = 0.;
end = "";
}; 
if (arg=="--flagflux") {
var_f = strtof(argv[i+1], &end);
end_str = end;
if ((end_str.size()==0)&&(abs(var_f) < 2)){
flag_flux = abs(int(var_f));
cout << "Accoring to '--flagflux',   Flag flux        was reset to " <<flag_flux<<"\n";
};
var_f = 0.;
end = "";
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
};
cout << "Ouput BOS flag  " << flag_bos << "  - no BOS output as 'make nobos' was used (for BOS output compile 'make bos') " << "\n";
#endif     
    
#ifdef BOS      
if ((flag_bos == 1)||(flag_bos == 2)) cout << "BOS output file name " << out_bos_file <<"\n";
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
     
//------------------
void cmdl_help(){
cout << "	option 		default		comment\n";
cout << "GENERAL:\n";
cout <<"	--seed		n/a   	Random seed, by default taken from time(NULL) \n";
cout <<"	--trig		10000 	Number of events to be generated \n";

cout << "KINEMATICS:\n";
cout <<"	--ebeam		10.6    Beam energy (GeV)\n";
cout <<"	--wmin		1.4    	Wmin   (GeV)\n";
cout <<"	--wmax		2.3   	Wmax   (GeV)\n";
cout <<"	--q2min		0.05 	Q2min  (GeV^2)\n";
cout <<"	--q2max		3.5   	Q2max  (GeV^2)\n";
cout <<"	--thmin		1.0   	Theta  min of scattered electron (deg)\n";
cout <<"	--thmax		50.   	Theta  max of scattered electron (deg)\n";
cout <<"	--emin		0.1   	Energy min of scattered electron (GeV)\n";


cout << "VERTEX:\n";
cout <<"	--trad		0.6    	Target radius  (cm)\n";
cout <<"	--tlen		2.0    	Target length  (cm)\n";
cout <<"	--toff		0.0    	Target offset  (cm)\n";

cout << "FOR RAD MODE:\n";
cout <<"	--tden		0.071	Target density   (g/cm^3)\n";
cout <<"	--trdln		888.03	Target radlenth  (cm)\n";
cout <<"	--tz		1    	Target Z	    \n";
cout <<"	--ta		1    	Target A	    \n";

cout <<"	--twlen		15.   	Target windows thickness (each)  (um)	   \n";
cout <<"	--twflen	15.   	Target windows thickness (final) (um) 	   if diff from initial\n";
cout <<"	--twden		2.699	Target windows density   (each)	 (g/cm^3)  \n";
cout <<"	--twfden	2.699	Target windows density   (final) (g/cm^3)  if diff from initial\n";
cout <<"	--twrdln	8.897	Target windows radlenth  (each)  (cm)      \n";
cout <<"	--twfrdln	8.897	Target windows radlenth  (final) (cm)	   if diff from initial\n";
cout <<"	--twz		13    	Target windows Z	 (each)  	  \n";
cout <<"	--twfz		13   	Target windows Z	 (final) 	   if diff from initial\n";
cout <<"	--twa		27    	Target windows A	 (each)  \n";
cout <<"	--twfa		27    	Target windows A	 (final) 	   if diff from initial\n";

cout << "FLAGS:\n";

cout <<"	--flagbos	0/1    	Output BOS file, 0 - no, (for 'make nobos' always 0)\n";
cout <<"	                                      	 1 - MCTK, MCVX banks, (for 'make bos' 1 by default)\n";
cout <<"	                                       	 2 - PART bank\n";
cout <<"	--flaglund	1    	Output LUND file, 0 - no, 1 - yes\n";
cout <<"	--flagrad	0    	Rad mode, 0 - no rad eff,\n";
cout <<"                                   	  1 - rad eff without straggling,\n";
cout <<"                                   	  2 - rad eff with    straggling\n";
cout <<"	--flagfermi	0     	Fermi smearing, 0 - no, 1 - yes\n";
cout <<"	--flagflux	1     	Multiply by Flux, 0 - no  (virtual photoproduction),\n";
cout <<"                                   		  1 - yes (electroproduction)\n";


cout << "OUTPUT:\n";
cout <<"	--bosname 	out.bos		BOS  output name (n/a for 'make nobos')  \n";
cout <<"	--lundname 	twopeg.dat 	LUND output name	    \n";
};     
     
     
     
     
     
     
