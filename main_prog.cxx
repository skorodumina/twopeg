#include "TROOT.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TObject.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <iostream>

#include "input.h"
#include "hist_def.h"
#include "read_xsect_files.h"
#include "read_fit_param_files.h"
#include "out_file_fill.h"
#include "out_file_open.h"
#include "out_file_close.h"
#include "anti_rot.h"
#include "get_xsect_ripani.h"
#include "get_xsect_near_threshold.h"

#include "get_xsect_golovach.h"
#include "get_xsect_fedotov.h"
#include "get_xsect_rip_fed_join.h"
#include "get_xsect_w16_18_lowq2_fit.h"
#include "get_xsect_w18_25.h"
#include "get_xsect_w25_30.h"
#include "get_xsect_w30_45.h"
#include "interpol_int.h"
#include "radcorr.h"
#include "fermi_bonn.h"
#include "fermi_rot.h"
#include "fermi_anti_rot.h"

#include "rot.h"
#include "global.h"
#include <stdlib.h>
#include <time.h>
#include <TLorentzVector.h>
#include <sstream>
#include <TRandom3.h>
#include <fstream>


 using namespace std;      
     

//Byckling function declaration
     Float_t G_BYCKLING(Float_t x, Float_t y, Float_t z, Float_t u, Float_t v, Float_t w) {
     return x*x*y+x*y*y+z*z*u+z*u*u+v*v*w+v*w*w+x*z*w+x*u*v+y*z*v+y*u*w-x*y*(z+u+v+w)-z*u*(x+y+v+w)-v*w*(x+y+z+u);
     };

//This is the function that performs the generation. It is declared here and executes below.   
void generate(Int_t argc, char *argv[]);


int main(int argc, char** argv){
std::string arg;

if (argc < 2) generate(argc,argv);

if (argc >= 2){
arg =   argv[1];
if (arg == "--help") cmdl_help();
if (!(arg == "--help"))  generate(argc,argv);
};

return 0;
};


//This is the function that performs the generation.
void generate(Int_t argc, char *argv[]) {

global();

    Float_t sigma_t_final = 0.;
    Float_t sigma_l_final = 0.; 
    Float_t sigma_c2f_final = 0.;
    Float_t sigma_s2f_final = 0.;
    Float_t sigma_cf_final = 0.;
    Float_t sigma_sf_final = 0.;
      
    Float_t sigma_t_final_1, sigma_l_final_1 , sigma_c2f_final_1, sigma_s2f_final_1,sigma_cf_final_1 ,sigma_sf_final_1 ;
    Float_t sigma_t_final_2, sigma_l_final_2 , sigma_c2f_final_2, sigma_s2f_final_2,sigma_cf_final_2 ,sigma_sf_final_2; 
    Float_t sigma_total_1, sigma_total_2;
    
    Float_t E_beam, E_beam_fermi, theta_rot2,E_E_prime_ferm,Theta_e_prime_ferm,W_ferm;
    Float_t Wnew, Q2new, E_beam_new, W_tmp;
    Float_t W,Q2,phi_e,W_old,Q2nodata,Q2_old, Q2lim1, Q2lim2;
    Float_t nu,E_E_prime,Theta_e_prime,E_E_prime_new;
    Float_t M1,M2,M3;
    Float_t inv_m12,inv_m23,inv_m13,th_hadr,alph_hadr,ph_hadr,s12,s23,s13,en1,en2,en3,mag1,mag2,mag3; 
    Float_t inv_m12_tst,inv_m23_tst,inv_m13_tst,th_hadr_tst,alph_hadr_tst,ph_hadr_tst;  
    Float_t z_EL, x_EL, y_EL,r_vert, phi_vert;
    Float_t e_rad_phot, cr_rad_fact;
    Short_t flag_seed = 0;
    Int_t seed = 0;
    
    
    Float_t V_flux = 0.;
    Float_t alph_const = 1./137.035;
    
    const Float_t phi_e_min = 0;
    const Float_t phi_e_max = 2*M_PI;
    
    TLorentzVector P4_E_prime, P4_PIP,P4_Pfin,P4_PIM;
    TLorentzVector P4_Eini, P4_Eini_new, P4_E_prime_new, P4_Pini;
    TLorentzVector P4_0_miss, P4_pim_miss; 
    TLorentzVector P4_0_miss_2, P4_pim_miss_2; 
    TLorentzVector P4_0_miss_fermi, P4_pim_miss_fermi; 
    TLorentzVector P4_0_miss_fermi_2, P4_pim_miss_fermi_2; 
    TLorentzVector P4_E_prime_boosted;
    TLorentzVector P4_Pini_fermi;
    
    Int_t dummy_int;
    Int_t k=0;
   
    ostringstream qqq;
    
//This is the environment variable that points to the "data" folder with the cross section files.
//By default it points to the current directory (in which the EG executes). 
data_dir = getenv("TWOPEG_DATA_DIR");
data_dir_2pi << data_dir << "/";
cout << "TWOPEG_DATA_DIR is " << data_dir_2pi.str() << endl;
 
 
//This needed for taking masses of the particles from pdg_table located in ROOT_DIR
const char *HOME_ROOT;
const char *HOME_ROOT1;
//HOME_ROOT = getenv("ROOT");
system("root_home=`root-config --etcdir`");
HOME_ROOT = getenv("root_home");
ostringstream ROOT_DIR;
ROOT_DIR << HOME_ROOT << "/pdg_table.txt";

TDatabasePDG *pdg = new TDatabasePDG();
pdg->ReadPDGTable(ROOT_DIR.str().c_str());
TParticlePDG *part1 = new TParticlePDG();
part1 = pdg->GetParticle("proton");
MP= part1->Mass();
part1 = pdg->GetParticle("pi+");
MPIP= part1->Mass();  
part1 = pdg->GetParticle("pi-");
MPIM= part1->Mass();  
part1 = pdg->GetParticle("e-");
Me= part1->Mass();  


//Reading input parameters
if (argc<2) input_stream(E_beam);
if (argc>=2) input_cline(E_beam, argc,argv, flag_seed, seed);
inp_couts(E_beam);

//Reading diff cross section from the tables in .dat files (filling out GLOBAL arrays)
read_xsect_files();
//Reading fit parameterms, which are needed for cross section extrapolation
read_fit_param_files();
     
//Reasonably changing max&min limits of kinematical variables if needed
cout << "___________________________________________\n\n";
    
if (W_min < (1.2375)) {
W_min = 1.2375;
cout << "Minimum W  has been changed to " << W_min << "\n";
}; 
      
if (W_max*W_max > MP*MP +2.*MP*(E_beam - E_eprime_min)) {
W_max = sqrt(MP*MP +2.*MP*(E_beam - E_eprime_min));
};
   
Q2lim1 = 2*E_beam*sin(Theta_min*M_PI/180./2.)*sin(Theta_min*M_PI/180./2.);
Q2lim1 = Q2lim1*(2*E_beam*MP-W_max*W_max+MP*MP);
Q2lim1 = Q2lim1/(MP+2*E_beam*sin(Theta_min*M_PI/180./2.)*sin(Theta_min*M_PI/180./2.));
if (Q2lim1 <= 0.00005) Q2lim1 = 0.00005;

if ((Q2_min < Q2lim1)) {
Q2_min = Q2lim1;
cout << "Minimum Q2 has been changed to " << Q2_min << "\n";
};
  
Q2lim2 = 2*E_beam*sin(Theta_max*M_PI/180./2.)*sin(Theta_max*M_PI/180./2.);
Q2lim2 = Q2lim2*(2*E_beam*MP-W_min*W_min+MP*MP);
Q2lim2 = Q2lim2/(MP+2*E_beam*sin(Theta_max*M_PI/180./2.)*sin(Theta_max*M_PI/180./2.));
   
if (Q2_max > Q2lim2) {
Q2_max = Q2lim2;
cout << "Maximum Q2 has been changed to " << Q2_max << "\n";
};
   
if ((W_max*W_max > MP*MP +2.*MP*(E_beam - E_eprime_min) -Q2_min)&&(W_max <= 4.5375)) {
W_max = sqrt(MP*MP +2.*MP*(E_beam - E_eprime_min) -Q2_min);
cout << "Maximum W  has been changed to " << W_max << "\n";
};

if (W_max > 4.5375) {
W_max = 4.5375;
cout << "Maximum W  has been changed to " << W_max << "\n";
};

//Defining some histograms
hist_def(E_beam); 
  
//Defining the 4-vectors of the initial particles (e and p)
P4_Pini.SetXYZT(0.,0.,0.,MP); 
P4_Eini.SetXYZT(0.,0.,E_beam,E_beam);
    
     
//Chosing the variable set  
//Second set of variables. FOR 1-PIM, 2-PIP, 3-P
     M1 = MPIM;
     M2 = MPIP;
     M3 = MP; 
//First set of variables. FOR 1-P, 2-PIP, 3-PIM
//   M1 = MP;
//   M2 = MPIP;
//   M3 = MPIM;
//Third set of variables. FOR 1-PIP, 2-PIM, 3-P
//   M1 = MPIP;
//   M2 = MPIM;
//   M3 = MP;
    
//open input file
out_file_open();


cout <<"\n";
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%!I. GENERATION STARTS!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

 cout << "____________EVENT GENERATION:______________\n\n";
 
if (flag_seed == 0) seed = time(NULL);
srand (seed);
cout << "RandomSeedActuallyUsed: " << seed << endl;       
       
 TRandom3 ph_e_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 th_hadr_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 W_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 Q2_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 z_EL_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 alph_hadr_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.)); 
 TRandom3 s12_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 s23_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));	  
 TRandom3 ph_hadr_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));	   
 TRandom3 r_vert_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 phi_vert_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
//This is for Fermi motion only 
 TRandom3 fermi_R_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 fermi_R1_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 fermi_theta_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 fermi_phi_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.)); 
//This is for radeff 
 TRandom3 phot_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.)); 
 TRandom3 hardini_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.));
 TRandom3 hardfin_rndm(UInt_t(((float) rand() / (float)(RAND_MAX))*4000000000.)); 
 
// Start to generate electrons    
for (Int_t i=1; i<=Nevents; i++) {
    
phi_e =ph_e_rndm.Uniform(phi_e_min,phi_e_max);
z_EL = z_EL_rndm.Uniform(Targ_off  - Targ_len/2.,Targ_off  + Targ_len/2.);

r_vert = r_vert_rndm.Uniform(0.,Targ_rad);
phi_vert = phi_vert_rndm.Uniform(0.,6.28318);
x_EL = r_vert*cos(phi_vert);
y_EL = r_vert*sin(phi_vert);

alph_hadr = alph_hadr_rndm.Uniform(0.,6.28318);
ph_hadr = ph_hadr_rndm.Uniform(0.,6.28318);    

do {
    
k++;
if ((k % 1000) == 0) cout  << "N generated = " << k << "; N accepted = " << i << "\n";

W = W_rndm.Uniform(W_min,W_max);
Q2 = Q2_rndm.Uniform(Q2_min,Q2_max);
  
nu=(W*W+Q2-MP*MP)/2./MP;
E_E_prime=E_beam-nu;

    Q2lim1 = 2*E_beam*sin(Theta_min*M_PI/180./2.)*sin(Theta_min*M_PI/180./2.);
    Q2lim1 = Q2lim1*(2*E_beam*MP-W*W+MP*MP);
    Q2lim1 = Q2lim1/(MP+2*E_beam*sin(Theta_min*M_PI/180./2.)*sin(Theta_min*M_PI/180./2.));

    Q2lim2 = 2*E_beam*sin(Theta_max*M_PI/180./2.)*sin(Theta_max*M_PI/180./2.);
    Q2lim2 = Q2lim2*(2*E_beam*MP-W*W+MP*MP);
    Q2lim2 = Q2lim2/(MP+2*E_beam*sin(Theta_max*M_PI/180./2.)*sin(Theta_max*M_PI/180./2.));
    
} while ((nu>E_beam-E_eprime_min)||(Q2>Q2lim2)||(Q2<Q2lim1)||(isnan(acos((Q2+2.*E_beam*nu)/2./E_beam/(sqrt(Q2+nu*nu))))));

//isnan(acos(1.-Q2/E_beam/E_E_prime/2.))
E_E_prime=E_beam-nu;
Theta_e_prime = acos(1.-Q2/E_beam/E_E_prime/2.);
    
//The four-momentum of the scattered electron in the Lab frame 
P4_E_prime.SetXYZT(E_E_prime*cos(phi_e)*sin(Theta_e_prime),E_E_prime*sin(phi_e)*sin(Theta_e_prime),E_E_prime*cos(Theta_e_prime),E_E_prime);
    
E_beam_new = E_beam;

W_old = W;
Q2_old = Q2;

Wnew = W;
W_ferm = W;
Q2new = Q2;

//Radiative effects
if ((flag_radmod == 1)||(flag_radmod == 2)){

radcorr(phot_rndm.Uniform(0.,1.),hardini_rndm.Uniform(0.,1.),hardfin_rndm.Uniform(0.,1.), E_beam,Q2,W,Wnew,Q2new,E_beam_new,e_rad_phot,cr_rad_fact);
h_eradgam->Fill(e_rad_phot,1.);
W = Wnew;
Q2 = Q2new;

};//end if rad mode
   
//for rad_eff !!!
P4_Eini_new.SetXYZT(0.,0.,E_beam_new, E_beam_new);

E_E_prime_new = E_beam_new - (Wnew*Wnew+Q2new-MP*MP)/2./MP;

//The four-momentum of the scattered electron in the Lab frame if rad eff happened   
P4_E_prime_new.SetXYZT(E_E_prime_new*cos(phi_e)*sin(Theta_e_prime),E_E_prime_new*sin(phi_e)*sin(Theta_e_prime),E_E_prime_new*cos(Theta_e_prime),E_E_prime_new);

E_beam_fermi = E_beam_new;

//for fermi motion    
if (flag_fermi == 1) {
do {

fermi_bonn(fermi_R_rndm.Uniform(0.,0.9999679),fermi_R1_rndm.Uniform(0.,1.),fermi_theta_rndm.Uniform(-1.,1.),fermi_phi_rndm.Uniform(0.,2.*M_PI));

//The four-momentum of the moving initial proton in the Lab frame
P4_Pini_fermi.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 
W_ferm = (P4_Pini_fermi + P4_Eini_new - P4_E_prime_new).Mag();
  
} while (W_ferm < 1.2375 );

//for fermi motion
fermi_rot(E_beam_fermi,theta_rot2,E_beam_new,P4_E_prime_new,P4_E_prime_boosted); 
phi_e = P4_E_prime_boosted.Phi();
W = W_ferm;
h_e_beam_eff->Fill(W_ferm,1.);
};
  
//  W_tmp = W;
// W = 1.6125;
dummy_int = 0;

do {

s12 = s12_rndm.Uniform((M1+M2)*(M1+M2),(W-M3)*(W-M3));
s23 = s23_rndm.Uniform((M2+M3)*(M2+M3),(W-M1)*(W-M1));
inv_m12 = sqrt(s12);
inv_m23 = sqrt(s23); 


//this variables are calculated here for the check of correct near-boundaries generation - conditions in while (...)
s13 = W*W+M1*M1+M2*M2+M3*M3-s12-s23;
en1 = (W*W+M1*M1-s23)/2./W;
en2 = (W*W+M2*M2-s13)/2./W;
en3 = (W*W+M3*M3-s12)/2./W;
mag1 = sqrt(en1*en1 - M1*M1);
mag2 = sqrt(en2*en2 - M2*M2);
mag3 = sqrt(en3*en3 - M3*M3);

dummy_int++;

} while ((G_BYCKLING(inv_m12*inv_m12,inv_m23*inv_m23,W*W,M2*M2,M1*M1,M3*M3) > 0.)||(isnan(acos((M1*M1+M2*M2+2*en1*en2-s12)/2./mag1/mag2)))||(isnan(acos((M1*M1+M3*M3+2*en1*en3-s13)/2./mag1/mag3)))||(en1 < M1)||(en2 < M2)||(en3 < M3)||(sqrt(s13) < M1+M3)||(sqrt(s13)>W-M2));
  

th_hadr = acos(th_hadr_rndm.Uniform(-1.,1.));


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%!II. WEIGHTS CALCULATION!%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Note that weights are calculated in assumption of the 2nd set of variables,i.e.
//s12 -> s_pim_pip and s23 -> s_pip_pr

Q2nodata = Q2;

if (Q2 > 1.299)Q2 = 1.299;
if (Q2 < 0.0005)Q2 = 0.0005;

sigma_t_final = 0.;
sigma_l_final = 0.;
sigma_c2f_final = 0.;
sigma_s2f_final = 0.;
sigma_cf_final = 0.;
sigma_sf_final = 0.;


//Getting cross section in given generated (W, Q2, s12, s23, theta, alpha)-point
sigma_total = 0.;

if ((W>=1.4125)&&(W<=1.8125)&&(Q2>=0.65)&&(Q2<=1.3)) {
get_xsect_ripani(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);
};

if ((W>=1.6125)&&(W<=1.8125)&&(Q2>0.000001)&&(Q2<0.65)) {
get_xsect_w16_18_lowq2_fit(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);
};

if (((W>=1.3125)&&(W<=1.4375)&&(Q2>0.275)&&(Q2<0.575))||((W>=1.4125)&&(W<=1.4875)&&(Q2>0.275)&&(Q2<0.525))||((W>=1.4125)&&(W<=1.5125)&&(Q2>0.275)&&(Q2<0.425))||((W>=1.4125)&&(W<=1.5375)&&(Q2>0.225)&&(Q2<0.375))||((W>=1.4125)&&(W<=1.5625)&&(Q2>0.225)&&(Q2<0.275))){
 get_xsect_fedotov(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr, sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final);
};

 if (((W>=1.4125)&&(W<=1.4375)&&(Q2>=0.575)&&(Q2<=0.65))||((W>=1.4375)&&(W<=1.4875)&&(Q2>=0.525)&&(Q2<=0.65))||((W>=1.4875)&&(W<=1.5125)&&(Q2>=0.425)&&(Q2<=0.65))||((W>=1.5125)&&(W<=1.5375)&&(Q2>=0.325)&&(Q2<=0.65))||((W>=1.5375)&&(W<=1.5625)&&(Q2>=0.275)&&(Q2<=0.65))||((W>=1.5625)&&(W<=1.5875)&&(Q2>=0.225)&&(Q2<=0.65))||((W>=1.3125)&&(W<1.4125)&&(Q2>=0.575)&&(Q2<=1.3))||((W>=1.3125)&&(W<1.5125)&&(Q2>=0.0002)&&(Q2<=0.275))||((W>=1.5125)&&(W<=1.5875)&&(Q2>=0.0002)&&(Q2<=0.225))){
get_xsect_rip_fed_join(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final); 
};

if ((W>=1.2375)&&(W<1.3125)&&(Q2>0.00002)&&(Q2<1.3)) {
get_xsect_near_threshold(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );
};

if ((W>=1.8375)&&(W<=2.5375)&&(Q2>0.00002)&&(Q2<1.3)) {
get_xsect_w18_25(E_beam,Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );
};

if ((W>=2.5875)&&(W<=3.0375)&&(Q2>0.00002)&&(Q2<1.3)) {
get_xsect_w25_30(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );
};


if ((W>=1.5875)&&(W<=1.6125)&&(Q2>0.00001)&&(Q2<0.65)) {
get_xsect_w16_18_lowq2_fit(Q2, 1.6125, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

get_xsect_rip_fed_join(Q2,1.5875, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1); 
 
sigma_t_final = 1./0.025;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(1.5875-W)+sigma_t_final_1*fabs(1.6125-W));

sigma_l_final = 1./0.025;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(1.5875-W)+sigma_l_final_1*fabs(1.6125-W));

sigma_c2f_final = 1./0.025;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_2*fabs(1.5875-W)+sigma_c2f_final_1*fabs(1.6125-W));

sigma_s2f_final = 1./0.025;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_2*fabs(1.5875-W)+sigma_s2f_final_1*fabs(1.6125-W));

sigma_cf_final = 1./0.025;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_2*fabs(1.5875-W)+sigma_cf_final_1*fabs(1.6125-W));

sigma_sf_final = 1./0.025;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_2*fabs(1.5875-W)+sigma_sf_final_1*fabs(1.6125-W));
};


if ((W>=1.8125)&&(W<=1.8375)&&(Q2>0.00002)&&(Q2<=1.3)) {

if ((Q2>0.00002)&&(Q2<=0.65)) get_xsect_w16_18_lowq2_fit(Q2, 1.8125, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

if ((Q2>0.65)&&(Q2<=1.3)) get_xsect_ripani(Q2, 1.8125, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

get_xsect_w18_25(E_beam,Q2, 1.8375, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1); 

sigma_t_final = 1./0.025;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(1.8375-W)+sigma_t_final_1*fabs(1.8125-W));

sigma_l_final = 1./0.025;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(1.8375-W)+sigma_l_final_1*fabs(1.8125-W));

sigma_c2f_final = 1./0.025;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_2*fabs(1.8375-W)+sigma_c2f_final_1*fabs(1.8125-W));

sigma_s2f_final = 1./0.025;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_2*fabs(1.8375-W)+sigma_s2f_final_1*fabs(1.8125-W));

sigma_cf_final = 1./0.025;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_2*fabs(1.8375-W)+sigma_cf_final_1*fabs(1.8125-W));

sigma_sf_final = 1./0.025;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_2*fabs(1.8375-W)+sigma_sf_final_1*fabs(1.8125-W));
};



if ((W>=2.5375)&&(W<=2.5875)&&(Q2>0.00001)&&(Q2<1.3)) {
get_xsect_w18_25(E_beam,Q2, 2.5375, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2,sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

get_xsect_w25_30(Q2,2.5875, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1); 
 
sigma_t_final = 1./0.05;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(2.5875-W)+sigma_t_final_1*fabs(2.5375-W));

sigma_l_final = 1./0.05;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(2.5875-W)+sigma_l_final_1*fabs(2.5375-W));

sigma_s2f_final = 1./0.05;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_2*fabs(2.5875-W)+sigma_s2f_final_1*fabs(2.5375-W));

sigma_c2f_final = 1./0.05;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_2*fabs(2.5875-W)+sigma_c2f_final_1*fabs(2.5375-W));

sigma_cf_final = 1./0.05;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_2*fabs(2.5875-W)+sigma_cf_final_1*fabs(2.5375-W));

sigma_sf_final = 1./0.05;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_2*fabs(2.5875-W)+sigma_sf_final_1*fabs(2.5375-W));
};


 if ((W>=3.1375)&&(W<=4.5375)&&(Q2>0.00001)&&(Q2<1.3)) {
get_xsect_w30_45(Q2, W, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final, sigma_l_final, sigma_c2f_final,sigma_s2f_final,sigma_cf_final,sigma_sf_final );
};

if ((W>=3.0375)&&(W<=3.1375)&&(Q2>0.00001)&&(Q2<1.3)) {

get_xsect_w25_30(Q2,3.0375, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_2, sigma_l_final_2, sigma_c2f_final_2,sigma_s2f_final_2,sigma_cf_final_2,sigma_sf_final_2);

get_xsect_w30_45(Q2, 3.1375, s12,s23, th_hadr, alph_hadr, ph_hadr,sigma_t_final_1, sigma_l_final_1, sigma_c2f_final_1,sigma_s2f_final_1,sigma_cf_final_1,sigma_sf_final_1 );


sigma_t_final = 1./0.1;
sigma_t_final = sigma_t_final*(sigma_t_final_2*fabs(3.1375-W)+sigma_t_final_1*fabs(3.0375-W));

sigma_l_final = 1./0.1;
sigma_l_final = sigma_l_final*(sigma_l_final_2*fabs(3.1375-W)+sigma_l_final_1*fabs(3.0375-W));

sigma_c2f_final = 1./0.1;
sigma_c2f_final = sigma_c2f_final*(sigma_c2f_final_2*fabs(3.1375-W)+sigma_c2f_final_1*fabs(3.0375-W));

sigma_s2f_final = 1./0.1;
sigma_s2f_final = sigma_s2f_final*(sigma_s2f_final_2*fabs(3.1375-W)+sigma_s2f_final_1*fabs(3.0375-W));

sigma_cf_final = 1./0.1;
sigma_cf_final = sigma_cf_final*(sigma_cf_final_2*fabs(3.1375-W)+sigma_cf_final_1*fabs(3.0375-W));

sigma_sf_final = 1./0.1;
sigma_sf_final = sigma_sf_final*(sigma_sf_final_2*fabs(3.1375-W)+sigma_sf_final_1*fabs(3.0375-W));
};

//Extending to the region Q2 > 1.3 GeV2
if (Q2nodata > 1.299){

sigma_t_final = sigma_t_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_l_final = sigma_l_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_c2f_final = sigma_c2f_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_s2f_final = sigma_s2f_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_cf_final = sigma_cf_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);
sigma_sf_final = sigma_sf_final*Func_q2_dep(Q2nodata)/Func_q2_dep(1.299);

};

Q2 = Q2nodata;

if ((W < 1.2375)||(W > 4.5375)){
sigma_t_final = 0.;
sigma_l_final = 0.;
sigma_c2f_final = 0.;
sigma_s2f_final = 0.;
sigma_cf_final = 0.;
sigma_sf_final = 0.;
};


//calculating sigma_total from different sigmas, eps_l and eps_t
Float_t eps_l,eps_t,nu_g,theta_el;

nu_g = (W*W + Q2 - MP*MP)/2./MP;


theta_el = acos(1.- Q2/E_beam_fermi/(E_beam_fermi - nu_g)/2.);
eps_t = 1./(1.+ 2.*(1. + nu_g*nu_g/Q2)*tan(theta_el/2.)*tan(theta_el/2.));
eps_l = Q2*eps_t/nu_g/nu_g;

sigma_total =0.;

if ((isnan(eps_l))||(isnan(eps_t))) cout << eps_l<< " "<< eps_t<<" eps in nan!\n";

if  (!(eps_l>0.)&&!(eps_l<0)) eps_l = 0.;
if  (!(eps_t>0.)&&!(eps_t<0)) eps_t = 0.;

//combining structure function into the full cross section
sigma_total = sigma_t_final;
sigma_total = sigma_total + eps_l*sigma_l_final;
sigma_total = sigma_total + eps_t*(sigma_c2f_final*cos(2.*ph_hadr) + sigma_s2f_final*sin(2.*ph_hadr));
sigma_total = sigma_total + sqrt(2.*eps_l*(eps_t+1))*(sigma_cf_final*cos(ph_hadr) + sigma_sf_final*sin(ph_hadr));

//Adding additional rad corr weight factor, if needed
if ((flag_radmod == 1)||(flag_radmod == 2)) {
sigma_total = sigma_total*cr_rad_fact;
if (isnan(cr_rad_fact)) cout<< sigma_total<<" "<<cr_rad_fact<<" rad corr factor is nan!\n";
};


//multiply sigma_total by virtual photon flux
if ((flag_flux == 1)&&(flag_fermi == 0)){
nu_g = (W_old*W_old + Q2_old - MP*MP)/2./MP;
theta_el = acos(1.- Q2_old/E_beam/(E_beam - nu_g)/2.);
eps_t = 1./(1.+ 2.*(1. + nu_g*nu_g/Q2_old)*tan(theta_el/2.)*tan(theta_el/2.));

V_flux = alph_const/4./M_PI;
V_flux = V_flux/E_beam/E_beam/MP/MP;
V_flux = V_flux/(1.-eps_t)/Q2_old;
V_flux = V_flux*W_old*(W_old*W_old-MP*MP);

sigma_total = sigma_total*V_flux;    
};


if ((flag_flux == 1)&&(flag_fermi == 1)){
nu_g = (W*W + Q2 - MP*MP)/2./MP;
theta_el = acos(1.- Q2/E_beam_fermi/(E_beam_fermi - nu_g)/2.);
eps_t = 1./(1.+ 2.*(1. + nu_g*nu_g/Q2)*tan(theta_el/2.)*tan(theta_el/2.));
//cout << eps_t << " t\n";
V_flux = alph_const/4./M_PI;
V_flux = V_flux/E_beam_fermi/E_beam_fermi/MP/MP;
V_flux = V_flux/(1.-eps_t)/Q2;
V_flux = V_flux*W*(W*W-MP*MP);

sigma_total = sigma_total*V_flux;    
};

if ((isnan(sigma_total))||(isnan(V_flux))) cout<<W_old<< " "<<W<<" "<<Q2_old<< " "<< Q2<<" "<< sigma_total<<" "<<sigma_t_final<< " "<< sigma_l_final<<" "<< sigma_c2f_final<< " "<< eps_l<<" weight is nan!\n";

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%!III. OBTAINING THE FINAL PARTICLE FOUR-MOMENTA IN THE LAB FRAME!%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//   cout <<"qqq1  "<<th_hadr <<"  "<< ph_hadr <<"  "<<alph_hadr<<  "\n";
 
    
//Second set of varibles. FOR 1-PIM, 2-PIP, 3-P
anti_rot(W, Q2, phi_e, E_beam_fermi, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr,  MPIM, MPIP, MP, P4_PIM, P4_PIP,  P4_Pfin);
  
//First set of varibles. FOR 1-P, 2-PIP, 3-PIM
//anti_rot(W, Q2, phi_e, E_beam_fermi, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MP,MPIP, MPIM, P4_Pfin,P4_PIP, P4_PIM);
   
// Third set of variables. FOR 1-PIP, 2-PIM, 3-P 
//anti_rot(W, Q2, phi_e, E_beam_fermi, inv_m12, inv_m23, th_hadr, alph_hadr,  ph_hadr, MPIP,MPIM, MP, P4_PIP,P4_PIM, P4_Pfin);
	
//-------------------------------------------------------------------------------------  
    
//This is transformation back-to-cms, which is needed here only for testing purposes.

//Second set of varibles. FOR 1-PIM, 2-PIP, 3-P
if (flag_fermi ==0) rot(Q2,E_beam_fermi,P4_E_prime_new,P4_PIM,P4_PIP,P4_Pfin,inv_m12_tst, inv_m23_tst,th_hadr_tst,alph_hadr_tst,ph_hadr_tst);
if (flag_fermi ==1) rot(Q2,E_beam_fermi,P4_E_prime_boosted,P4_PIM,P4_PIP,P4_Pfin,inv_m12_tst, inv_m23_tst,th_hadr_tst,alph_hadr_tst,ph_hadr_tst);

//First set of varibles. FOR 1-P, 2-PIP, 3-PIM
/*if (flag_fermi ==0) rot(Q2,E_beam_fermi,P4_E_prime_new,P4_Pfin,P4_PIP,P4_PIM,inv_m12_tst,inv_m23_tst,th_hadr_tst,alph_hadr_tst,ph_hadr_tst);
if (flag_fermi ==1) rot(Q2,E_beam_fermi,P4_E_prime_boosted,P4_Pfin,P4_PIP,P4_PIM,inv_m12_tst,inv_m23_tst,th_hadr_tst,alph_hadr_tst,ph_hadr_tst);*/
  
//Third set of variables. FOR 1-PIP, 2-PIM, 3-P
/*if (flag_fermi ==0) rot(Q2,E_beam_fermi,P4_E_prime_new,P4_PIP,P4_PIM,P4_Pfin,inv_m12_tst,inv_m23_tst,th_hadr_tst,alph_hadr_tst,ph_hadr_tst);
if (flag_fermi ==1) rot(Q2,E_beam_fermi,P4_E_prime_boosted,P4_PIP,P4_PIM,P4_Pfin,inv_m12_tst,inv_m23_tst,th_hadr_tst,alph_hadr_tst,ph_hadr_tst);*/

//test
//if ((fabs(th_hadr-th_hadr_tst)>0.001)||(fabs(ph_hadr-ph_hadr_tst)>0.001)||(fabs(alph_hadr-alph_hadr_tst)>0.001)) cout <<"ALARM!  Four-momenta of the final particles are probably WRONG!\n";     
if (fabs(th_hadr-th_hadr_tst)>0.001) cout <<"ALARM!  Four-momenta of the final particles are probably WRONG!\n";
 
// cout <<"qqq2  "<<th_hadr_tst <<"  "<< ph_hadr_tst <<"  "<<alph_hadr_tst<<  "\n";
//------------------------------------------------------------

//In the Fermi mode - transformation of the momenta from quasiLab to Lab
if (flag_fermi ==1) fermi_anti_rot(W,Q2,E_beam_new,E_beam_fermi, phi_e,theta_rot2,P4_PIM,P4_PIP, P4_Pfin,P4_E_prime_new);

//cout << phi_e<<" ppp\n";	
//do not remove. this is the nan-check before writing all momenta to the output file
if ((isnan(P4_E_prime[0]))||(isnan(P4_E_prime[1]))||(isnan(P4_E_prime[2]))||(isnan(P4_E_prime[3]))) cout << P4_E_prime[0]<< " "<< P4_E_prime[1]<< " "<<P4_E_prime[2]<< " "<<P4_E_prime[3]<<" "<<W<<" "<<Q2<< " final electron is nan \n";
	
if ((isnan(P4_Pfin[0]))||(isnan(P4_Pfin[1]))||(isnan(P4_Pfin[2]))||(isnan(P4_Pfin[3]))) cout << P4_Pfin[0]<< " "<< P4_Pfin[1]<< " "<<P4_Pfin[2]<< " "<<P4_Pfin[3]<<" "<<W<<" "<<Q2<< " final proton is nan \n";
	
if ((isnan(P4_PIP[0]))||(isnan(P4_PIP[1]))||(isnan(P4_PIP[2]))||(isnan(P4_PIP[3]))) cout << P4_PIP[0]<< " "<< P4_PIP[1]<< " "<<P4_PIP[2]<< " "<<P4_PIP[3]<<" "<<W<<" "<<Q2<< " pip is nan \n";
		
if ((isnan(P4_PIM[0]))||(isnan(P4_PIM[1]))||(isnan(P4_PIM[2]))||(isnan(P4_PIM[3]))) cout << P4_PIM[0]<< " "<< P4_PIM[1]<< " "<<P4_PIM[2]<< " "<<P4_PIM[3]<<" "<<W<<" "<<Q2<< " pim is nan \n";
	
//cout << P4_PIP.Phi()<<" oo\n";	

//writing generated events into the desired input file+filling the root-tree with weights
out_file_fill(i,sigma_total, W, Q2,  P4_E_prime,P4_Pfin, P4_PIP,P4_PIM,z_EL,x_EL,y_EL);

//filling out histograms   
hist_fill(E_beam,W,W_old,Q2,Q2_old,phi_e,z_EL,s12,s23, th_hadr,alph_hadr,ph_hadr,sigma_t_final,sigma_l_final,eps_l,sigma_total,P4_E_prime,P4_Pfin, P4_PIP,P4_PIM,P4_Eini_new,P4_E_prime_new);
  
};//end event-loop
    
    
//closing output file including root-tree with weights
out_file_close();

//Writing histograms into the root file 
hist_write();
    
  
/*for (Short_t j=0;j<=26;j++){
 qqq.str("");
  //  qqq << "Q2_" << 0.0005+0.05*j <<  ".dat";
    qqq << "intsec_q2_" <<j <<  ".dat";
     std::ofstream ofs (qqq.str().c_str(), std::ofstream::out);
     qqq.str(""); 
 for (Short_t i=1;i<=71;i++){
  ofs  << 1.2375 +0.025*i<< "\n";    
  ofs <<h_int_crsect_t[j]->GetBinContent(i) << "\n";  
  ofs <<h_int_crsect_l[j]->GetBinContent(i) << "\n"; 
 };
  }; */


};
