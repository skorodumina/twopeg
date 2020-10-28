#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"
#include "interpol_int.h"
#include "get_xsect_14_18_lowq2_fit.h"
#include "TH1.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "TMath.h"

    //Min Energy of radiated photon for radcorr
    Float_t delta = 0.01;
    Float_t alpha = 1./137.035;
    TH1F *h;
    const  Int_t Npoints = 801;

    Float_t ARR_e_radhardini[Npoints];
    Float_t ARR_r_radhardini[Npoints];
    
    
     Float_t ARR_e_radhardfin[Npoints];
     Float_t ARR_r_radhardfin[Npoints];   
    
    
    
    //Calculating lengthes of the target and i/f windows using input-file information
    Float_t T_targ ;
    Float_t T_wi;
    Float_t T_wf;
    
    

using namespace std;


//Virtual photon flux-----------------------------------------------------------------------
Float_t Gamma_V_dE_dOm (Float_t E_beam, Float_t Wgen, Float_t Q2gen){

//Float_t alpha = 1./137.035;
Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Float_t theta_el = acos(1.- Q2gen/E_beam/(E_beam - E_gamma_lab)/2.);
//cout << theta_el << "  radcorr \n";
Float_t eps_t = 1./(1.+ 2.*(1. + E_gamma_lab*E_gamma_lab/Q2gen)*tan(theta_el/2.)*tan(theta_el/2.));
Float_t GV = alpha/4./M_PI/M_PI;

GV = GV/E_beam/Q2gen/(1.-eps_t)/MP;
GV = GV*(E_beam - E_gamma_lab)*(Wgen*Wgen - MP*MP);

return GV;
};


//Auxiliary B-function according to Mo&Tsai Appendix A, furmalae A.4 and A.5 page 226------
Float_t BZFUN(Float_t Z){
Float_t xi = (log(1440.)-2.*log(Z)/3.)/(log(183.)-log(Z)/3.);

return 4./3.*(1. + ((Z+1.)/(Z+xi))/(log(183.) - log(Z)/3.)/9.);

};


//Spence-function-----------------------------------------------------------------------
Float_t Spence(Float_t x){

Float_t value;
Float_t null = 1E-8;

if ((x>=null)&&(x<=1.)) {
   TF1 f("Spence Func", "-log(abs(1.-x))/x", null, x);
   ROOT::Math::WrappedTF1 wf1(f);
   ROOT::Math::GaussIntegrator ig;
   ig.SetFunction(wf1);
   ig.SetRelTolerance(1E-8);
   value = ig.Integral(null, x);
};

if (x>=1.) {
   TF1 f("Spence Func", "-log(abs(1.-x))/x", 1., x);
   ROOT::Math::WrappedTF1 wf1(f);
   ROOT::Math::GaussIntegrator ig;
   ig.SetFunction(wf1);
   ig.SetRelTolerance(1E-8);
   value = M_PI*M_PI/6. + ig.Integral(1., x);
  };
 
 
 if (x<=-1.*null) {
   TF1 f("Spence Func", "-log(abs(1.-x))/x",x, -1.*null);
   ROOT::Math::WrappedTF1 wf1(f);
   ROOT::Math::GaussIntegrator ig;
   ig.SetFunction(wf1);
   ig.SetRelTolerance(1E-8);
   value = -1.*ig.Integral(x, -1*null);
  };
  
 if ((x>=-1.*null)&&(x<=null)) value = 0.;

return value; 
 };



//This subroutine calculates the first term of the formula (IV.1) Mo&Tsai  (page 213)
//It includes:
//* calculation of two factors Delta_r and Delta_t (see section IV.A Mo and Tsai pages 213-214)
//* calculation of 2dim linearly interpolated intergated 2pi cross-section (in W-Q2 variables) including calculation sigma_t and sigma_l separately and their combonation with eps_L
//*  multiplying by the Gamma_V_dE_dOm - virtual photon flux for the case (E_f, Omega variables)
Float_t s1_radsoft(Float_t E_beam, Float_t Q2gen, Float_t Wgen){


Float_t sec_test_t,sec_test_l,sec_total;

Float_t E_gamma_lab, Delta_t, Delta_r, factor_radsoft, s1;


Float_t Btarg = BZFUN(Targ_Z);
Float_t Bwi = BZFUN(Twi_Z);
Float_t Bwf = BZFUN(Twf_Z);


E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 

Delta_t = (Bwi*T_wi + 0.5*Btarg*T_targ)*log(E_beam/delta)+(Bwf*T_wf + 0.5*Btarg*T_targ)*log((E_beam-E_gamma_lab)/delta);
Delta_t = -1.*Delta_t;

Delta_r = 28./9. - (13./6.)*log(Q2gen/Me/Me)+(log(E_beam/delta)+log((E_beam-E_gamma_lab)/delta))*(log(Q2gen/Me/Me)-1.) - Spence(-E_gamma_lab/(E_beam-E_gamma_lab))-Spence(E_gamma_lab/E_beam);

Delta_r = -1.*Delta_r*alpha/M_PI;

factor_radsoft = exp(Delta_t + Delta_r);

//Cross-section interpolation
interpol_int(Q2gen,Wgen,sec_test_t, sec_test_l);

sec_total = sec_test_t + getEpsL(E_beam,Wgen,Q2gen)*sec_test_l;

sec_total = sec_total*Gamma_V_dE_dOm(E_beam,Wgen,Q2gen);

s1 = factor_radsoft*sec_total;

return s1;
};


//-------------------------------------------
Float_t Omega_max_ini(Float_t E_beam, Float_t Q2gen, Float_t Wgen){

Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Float_t theta_el = acos(1.- Q2gen/E_beam/(E_beam - E_gamma_lab)/2.);
Float_t Ep = E_beam - E_gamma_lab;

//For the case of 1pi-treshold (see Mo&Tsai formula (A.17)&(A.18))
//Float_t E_beam_min = (MPIM*MPIM + 2.*MP*MPIM + 2.*MP*Ep)/(2.*MP - 2.*Ep*(1. - cos(theta_el)));

//For the case of 2pi-treshold (see Mo&Tsai formula (A.17)&(A.18))
Float_t E_beam_min = (2.*MPIM*MPIM + 2.*MP*MPIM + MP*Ep)/(MP - Ep*(1. - cos(theta_el)));

return E_beam - E_beam_min;

}; 


Float_t Rho_s_radhard(Float_t E_beam, Float_t Q2gen, Float_t Wgen, Float_t omega){


Float_t Btarg = BZFUN(Targ_Z);
Float_t Bwi = BZFUN(Twi_Z);
Float_t Bwf = BZFUN(Twf_Z);


Float_t E_beam_new = E_beam - omega;

Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Float_t Ep = E_beam - E_gamma_lab;


Float_t cs = (Bwi*T_wi + 0.5*Btarg*T_targ);
Float_t cp = (Bwf*T_wf + 0.5*Btarg*T_targ);
Float_t tr = alpha/M_PI/Btarg*(log(Q2gen/Me/Me) - 1.);
Float_t fs = Btarg*tr + cs;
Float_t fp = Btarg*tr + cp;
Float_t xs = E_beam_new/E_beam;

Float_t ts = alpha/M_PI*(0.5*(1. + xs*xs)*log(Q2gen/Me/Me) - xs);



Float_t Rho = 1./(E_beam - E_beam_new);
Rho = Rho*(ts+cs*(xs + 0.75*(1. - xs)*(1. - xs)));
Rho = Rho*pow(log(1./xs),fs);
Rho = Rho*pow(delta/Ep,0.5*fp);

return Rho; 
};

//This is for performing integration by Simpson method
Float_t DSIMPS_s2(){
Float_t Integral, h;
Integral = 0.;
for(Int_t i=2;i<Npoints;i += 2){

h = ARR_e_radhardini[i] - ARR_e_radhardini[i-2];
Integral += (ARR_r_radhardini[i-2] + 4.*ARR_r_radhardini[i-1] + ARR_r_radhardini[i])*h; 
};

Integral /= 6.0;

return Integral;
};


//This subroutine calculates the second term of the formula (IV.1) Mo&Tsai (page 213)
//It corresponds to the case when initial electron emits photon, so the beam energy changes
//The final electron is asumed unchanged!
//It includes:
//* determination of the limits of integration. Intergation is performed over the possible photon energy,
//so the lower limit is delta and the upper is calculated by subroutine Omega_max_ini according to the 
//reaction treshold - see Mo&Tsai formula (A.17)&(A.18)
//* the loop over the possible photon energies from min to max 
//------for each loop iteration we are doing the following:
//-------*calculating "new" beam energy 
//-------*calculating "new" W and Q2 values assuming the final electron unchanged!
//-------*picking the cross-section for "new" W and Q2  
//-------*calculating the rest of the integrand of the second term of the formula (IV.1) Mo&Tsai 
//--------by using the subroutin Rho_s_radhard
//*performing integration of the function Cross-section*Rho_s_radhard defined on the grid given by loop
Float_t s2_radhardini(Float_t E_beam,Float_t Q2gen, Float_t Wgen){

Float_t Omega_maxini =  Omega_max_ini(E_beam, Q2gen,Wgen);

if (Omega_maxini <= delta){
return 0;
};

Float_t omega_current,s2;
Float_t sec_test_t,sec_test_l,sec_total;
Float_t E_beam_new,w_new,q2_new;

Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Float_t theta_el = acos(1.- Q2gen/E_beam/(E_beam - E_gamma_lab)/2.);
Float_t Ep = E_beam - E_gamma_lab;


//loop
for (Int_t i=1;i<=Npoints;i++){

omega_current = delta + (Omega_maxini - delta)/(Npoints-1)*(i-1);
ARR_e_radhardini[i-1] = omega_current;


E_beam_new = E_beam - omega_current;

q2_new = 2*E_beam_new*Ep*(1. - cos(theta_el));
w_new = sqrt(MP*MP - q2_new + 2.*(E_beam_new - Ep)*MP);

sec_test_t = 0.;
sec_test_l =0.;


if ((w_new < 1.2375)||(w_new > 4.5375)) sec_total = 0.;



if ((w_new > 1.2375)&&(w_new < 4.5375)) {

//Cross-section interpolation
interpol_int(q2_new,w_new,sec_test_t, sec_test_l);
sec_total = sec_test_t + getEpsL(E_beam-omega_current,w_new,q2_new)*sec_test_l;
};

sec_total = sec_total*Gamma_V_dE_dOm(E_beam-omega_current,w_new,q2_new);
ARR_r_radhardini[i-1] = Rho_s_radhard(E_beam, Q2gen, Wgen,omega_current)*sec_total;

};

s2 = DSIMPS_s2();

return s2;
};


//---------------------------------------------------------------------
Float_t Omega_max_fin(Float_t E_beam, Float_t Q2gen, Float_t Wgen){

Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Float_t theta_el = acos(1.- Q2gen/E_beam/(E_beam - E_gamma_lab)/2.);
Float_t Ep = E_beam - E_gamma_lab;



//For the case of 1pi-treshold (see Mo&Tsai formula (A.17)&(A.19))
//Float_t E_p_max = (2.*MP*E_beam - 2.*MP*MPIM -MPIM*MPIM)/(2.*MP + 2.*E_beam*(1. - cos(theta_el)));

//For the case of 2pi-treshold (see Mo&Tsai formula (A.17)&(A.19))
Float_t E_p_max = (MP*E_beam - 2.*MP*MPIM - 2.*MPIM*MPIM)/(MP + E_beam*(1. - cos(theta_el)));

return E_p_max - Ep;
}; 
//----------------------------------------------------------------------

//-------------------------------------------------------------------------
Float_t Rho_p_radhard(Float_t E_beam, Float_t Q2gen, Float_t Wgen, Float_t omega){


Float_t Btarg = BZFUN(Targ_Z);
Float_t Bwi = BZFUN(Twi_Z);
Float_t Bwf = BZFUN(Twf_Z);

Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Float_t Ep = E_beam - E_gamma_lab;
Float_t E_p_new = Ep + omega;

    

Float_t cs = (Bwi*T_wi + 0.5*Btarg*T_targ);
Float_t cp = (Bwf*T_wf + 0.5*Btarg*T_targ);


Float_t tr = alpha/M_PI/Btarg*(log(Q2gen/Me/Me) - 1.);
Float_t fs = Btarg*tr + cs;
Float_t fp = Btarg*tr + cp;
Float_t xp = Ep/E_p_new;
Float_t tp = alpha/M_PI*(0.5*(1. + xp*xp)*log(Q2gen/Me/Me) - xp);



Float_t Rho = 1./(E_p_new - Ep);
Rho = Rho*(tp+cp*(xp + 0.75*(1. - xp)*(1. - xp)));
Rho = Rho*pow(log(1./xp),fp);
Rho = Rho*pow(delta/E_beam,0.5*fs);

return Rho; 
};
//-------------------------------------------------------------------------


//This is for performing integration by Simpson method
Float_t DSIMPS_s3(){
Float_t Integral, h;
Integral = 0.;
for(Int_t i=2;i<Npoints;i += 2){

h = ARR_e_radhardfin[i] - ARR_e_radhardfin[i-2];
Integral += (ARR_r_radhardfin[i-2] + 4.*ARR_r_radhardfin[i-1] + ARR_r_radhardfin[i])*h; 
};

Integral /= 6.0;

return Integral;
};



//This subroutine calculates the third term of the formula (IV.1) Mo&Tsai (page 213) 
//It corresponds to the case when final electron emits photon, so the energy of the final electron changes
//The initial electron is asumed unchanged!
//It includes:
//* determination of the limits of integration. Intergation is performed over the possible photon energy,
//so the lower limit is delta and the upper is calculated by subroutine Omega_max_fin according to the 
//reaction treshold - see Mo&Tsai formula (A.17)&(A.18)
//* the loop over the possible photon energies from min to max 
//------for each loop iteration we are doing the following:
//-------*calculating "new" energy of the final electron 
//-------*calculating "new" W and Q2 values assuming the initial electron unchanged!
//-------*picking the cross-section for "new" W and Q2  
//-------*calculating the rest of the integrand of the third term of the formula (IV.1) Mo&Tsai 
//--------by using the subroutin Rho_p_radhard
//*performing integration of the function Cross-section*Rho_p_radhard defined on the grid given by loop
Float_t s3_radhardfin(Float_t E_beam,Float_t Q2gen, Float_t Wgen){

Float_t Omega_maxfin =  Omega_max_fin(E_beam, Q2gen,Wgen);

if (Omega_maxfin <= delta){
return 0;
};

Float_t omega_current,s3;
Float_t sec_test_t,sec_test_l,sec_total;
Float_t Ep,E_p_new,w_new,q2_new;

Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Float_t theta_el = acos(1.- Q2gen/E_beam/(E_beam - E_gamma_lab)/2.);
Ep = E_beam - E_gamma_lab;


for (Int_t i=1;i<=Npoints;i++){

omega_current = delta + (Omega_maxfin - delta)/(Npoints-1)*(i-1);
ARR_e_radhardfin[i-1] = omega_current;

E_p_new = Ep + omega_current;

q2_new = 2.*E_beam*E_p_new*(1. - cos(theta_el));
w_new = sqrt(MP*MP - q2_new + 2.*(E_beam - E_p_new)*MP);

sec_test_t = 0.;
sec_test_l =0.;

if ((w_new < 1.2375)||(w_new >4.5375)) {
sec_total = 0.;
};


if ((w_new > 1.2375)&&(w_new < 4.5375)) {

//Cross-section interpolation
interpol_int(q2_new,w_new,sec_test_t, sec_test_l);
sec_total = sec_test_t + getEpsL(E_beam,w_new,q2_new)*sec_test_l;

};

sec_total = sec_total*Gamma_V_dE_dOm(E_beam,w_new,q2_new);
ARR_r_radhardfin[i-1] = Rho_p_radhard(E_beam, Q2gen, Wgen,omega_current)*sec_total;


};

s3 = DSIMPS_s3();

return s3;
};

void radcorr(Float_t R, Double_t R_ini, Double_t R_fin, Float_t E_beam, Float_t Q2gen, Float_t Wgen, Float_t &Wnew, Float_t &Q2new, Float_t &E_beam_new, Float_t &e_radgam, Float_t &cr_rad_fact){

//Calculating lengthes of the target and i/f windows using input-file information
Float_t sec_test_t,sec_test_l,sec_total;

if (flag_radmod == 1) {
T_targ  = 0.;
T_wi  = 0.;
T_wf  = 0.;
};


if (flag_radmod == 2) {
T_targ = Targ_len/Targ_radlen;
T_wi = Twi_thick/Twi_radlen*1E-4;
T_wf = Twf_thick/Twf_radlen*1E-4;
};
 
Double_t R1[1];
Double_t eran_hardini[1];

Double_t R2[1];
Double_t eran_hardfin[1];


Float_t s1 = s1_radsoft(E_beam, Q2gen,Wgen);

Float_t s2 = s2_radhardini(E_beam,Q2gen,Wgen);

Float_t s3 = s3_radhardfin(E_beam,Q2gen,Wgen);

//cout << "s1 = "<< s1_radsoft(E_beam, 0.9,1.35)<< ", s2 = "<< s2_radhardini(E_beam,0.9,1.35)<< ", s3 = "<<s3_radhardfin(E_beam,0.9,1.35)<<"\n";

//cout << s1<< " "<< s2<< " "<< s3<< "\n";

//E_beam_new = E_beam;
Float_t Ep_new;


//R is a random number from 0 to 1. It is an input parameter.
if (R < s1/(s1+s2+s3)){
Wnew = Wgen;
Q2new = Q2gen;
E_beam_new = E_beam;


Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Ep_new = E_beam - E_gamma_lab;
e_radgam = 0.;

};

if ((R >= s1/(s1+s2+s3))&&(R < (s1+s2)/(s1+s2+s3))){

TH1F*h_radhardini = new TH1F("h_radhardini","h_radhardini",800,ARR_e_radhardini[0],ARR_e_radhardini[Npoints-1]);

for (Int_t i=0;i<Npoints-1;i++){
h_radhardini->SetBinContent(i+1,(ARR_r_radhardini[i]+ARR_r_radhardini[i+1])/2.);
};
//R_ini is a random number from 0 to 1. It is an input parameter.
R1[0] = R_ini;

h_radhardini->GetQuantiles(1,eran_hardini,R1);
//cout << eran_hardini[0]<<" "<< ARR_e_radhardini[0]<<" "<<ARR_e_radhardini[Npoints-1] << "\n";


E_beam_new = E_beam - Float_t(eran_hardini[0]);

Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Float_t Ep = E_beam - E_gamma_lab;
Float_t theta_el = acos(1.- Q2gen/E_beam/(E_beam - E_gamma_lab)/2.);
Q2new = 2*E_beam_new*Ep*(1. - cos(theta_el));
Wnew = sqrt(MP*MP - Q2new + 2.*(E_beam_new - Ep)*MP);
Ep_new = Ep;

e_radgam = Float_t(eran_hardini[0]);

h_radhardini->Delete();
};



if ((R >= (s1+s2)/(s1+s2+s3))&&(R < 1.)){
TH1F*h_radhardfin = new TH1F("h_radhardfin","h_radhardfin",800,ARR_e_radhardfin[0],ARR_e_radhardfin[Npoints-1]);


for (Int_t i=0;i<Npoints-1;i++){
h_radhardfin->SetBinContent(i+1,(ARR_r_radhardfin[i]+ARR_r_radhardfin[i+1])/2.);
};
//R_fin is a random number from 0 to 1. It is an input parameter.
R2[0] = R_fin;

h_radhardfin->GetQuantiles(1,eran_hardfin,R2);

Float_t E_gamma_lab = (Wgen*Wgen+Q2gen-MP*MP)/2./MP; 
Float_t Ep = E_beam - E_gamma_lab;
Float_t theta_el = acos(1.- Q2gen/E_beam/(E_beam - E_gamma_lab)/2.);

Ep_new = Ep + Float_t(eran_hardfin[0]);

E_beam_new = E_beam;


Q2new = 2.*E_beam*Ep_new*(1. - cos(theta_el));
Wnew = sqrt(MP*MP - Q2new + 2.*(E_beam - Ep_new)*MP);


e_radgam = Float_t(eran_hardfin[0]);

h_radhardfin->Delete();
};
//cout <<Wgen<< " "<< s1/(s1+s2+s3)<<" "<<(s1+s2)/(s1+s2+s3)<<" R\n";
sec_test_t = 0.;
sec_test_l =0.;

if ((Wnew < 1.2375)||(Wnew > 4.5375)) {
sec_total = 0.;
cr_rad_fact = 0.;
};



if ((Wnew > 1.2375)&&(Wnew < 4.5375)) {
//Cross-section interpolation
interpol_int(Q2new,Wnew,sec_test_t, sec_test_l);

sec_total = sec_test_t + getEpsL(E_beam_new,Wnew,Q2new)*sec_test_l;
sec_total = sec_total*Gamma_V_dE_dOm(E_beam,Wgen,Q2gen);
cr_rad_fact = (s1+s2+s3)/sec_total;

};





};
