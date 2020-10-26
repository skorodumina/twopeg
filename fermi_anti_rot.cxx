#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"


using namespace std;

//----------------------------------------------------
//This subroutine performes the quasiLab-->Lab transformation, where quasiLab is the system, where the initial proton is at rest, while the initial electron moves along the Z-axis. The transformation is performed via the sequence of auxiliary systems: quasiLab(System 3)-->System 2-->System 1-->Lab.
//As an output it gives the four-momenta of the final hadrons (P4_1, P4_2, and P4_3) in the Lab frame.

 void fermi_anti_rot(Float_t W_ferm, Float_t Q2,Float_t E_beam, Float_t E_beam_fermi,Float_t phi_e, Float_t theta_rot2, TLorentzVector &P4_1, TLorentzVector &P4_2, TLorentzVector &P4_3, TLorentzVector P4_E_prime_Lab) {
 
//For the testing purposes let's construct the four-momenta of the incoming and scattered electons in the quasiLab frame 
TLorentzVector P4_Eini_qLab, P4_E_prime_qLab;
Float_t E_E_prime_new2,Theta_e_prime2;

//The four-momentum of the initial electron in the quasiLab frame 
P4_Eini_qLab.SetXYZT(0.,0.,E_beam_fermi,E_beam_fermi); 


E_E_prime_new2 = E_beam_fermi - (W_ferm*W_ferm+Q2-MP*MP)/2./MP; 
Theta_e_prime2 = acos(1.-Q2/E_beam_fermi/E_E_prime_new2/2.);

//The four-momentum of the scattered electron in the quasiLab frame 
P4_E_prime_qLab.SetXYZT(E_E_prime_new2*cos(phi_e)*sin(Theta_e_prime2),E_E_prime_new2*sin(phi_e)*sin(Theta_e_prime2),E_E_prime_new2*cos(Theta_e_prime2),E_E_prime_new2);
  
//Defining of the spatial angles of the Fermi momentum   
Float_t theta_fermi, phi_fermi;

theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;


//Defininh beta of the boost  
Float_t beta;

 beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);

 //--------------------------------------------
 //Transformation quasiLab-->System 2
P4_Eini_qLab.RotateY(-theta_rot2); 
P4_E_prime_qLab.RotateY(-theta_rot2);
P4_1.RotateY(-theta_rot2);
P4_2.RotateY(-theta_rot2);
P4_3.RotateY(-theta_rot2);
 

//Transformation System 2--> System 1  
P4_Eini_qLab.Boost(0,0,beta); 
P4_E_prime_qLab.Boost(0,0,beta);
P4_1.Boost(0,0,beta);
P4_2.Boost(0,0,beta);
P4_3.Boost(0,0,beta);

//Transformation System 1--> Lab
P4_Eini_qLab.RotateY(theta_fermi);
P4_Eini_qLab.RotateZ(phi_fermi); 
  
P4_E_prime_qLab.RotateY(theta_fermi);
P4_E_prime_qLab.RotateZ(phi_fermi);
    
P4_1.RotateY(theta_fermi);
P4_1.RotateZ(phi_fermi);

P4_2.RotateY(theta_fermi);
P4_2.RotateZ(phi_fermi);

P4_3.RotateY(theta_fermi);
P4_3.RotateZ(phi_fermi);

//Now, the four momenta of all final hadrons as well as the initial (P4_Eini_qLab) and scatteres electrons (P4_E_prime_qLab) are written in the usual Lab frame!
 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
//For the testing purposes let's compare the obtained four-momenta of the initial (P4_Eini_qLab) and scatteres electrons (P4_E_prime_qLab) with those already exsisted in the Lab frame. 

//The four-momentum of the initial electron in the Lab frame
TLorentzVector P4_Eini_Lab;
P4_Eini_Lab.SetXYZT(0.,0.,E_beam,E_beam); 

 
 if ((fabs(P4_E_prime_qLab[0]-P4_E_prime_Lab[0])>0.0001)||(fabs(P4_E_prime_qLab[1]-P4_E_prime_Lab[1])>0.0001)||(fabs(P4_E_prime_qLab[2]-P4_E_prime_Lab[2])>0.0001)||(fabs(P4_E_prime_qLab[2]-P4_E_prime_Lab[2])>0.0001))  cout << "ALARM! Wrong Lab-->quasiLab-->Lab momenta transformation! \n";

 if ((fabs(P4_Eini_qLab[0]-P4_Eini_Lab[0])>0.0001)||(fabs(P4_Eini_qLab[1]-P4_Eini_Lab[1])>0.0001)||(fabs(P4_Eini_qLab[2]-P4_Eini_Lab[2])>0.0001)||(fabs(P4_Eini_qLab[2]-P4_Eini_Lab[2])>0.0001)) cout << "ALARM! Wrong Lab-->quasiLab-->Lab momenta transformation! \n"; 
 
 
  
 };
 
