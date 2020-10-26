#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"


using namespace std;

//----------------------------------------------------
//This subroutine performs the Lab-->quasiLab transformation, where quasiLab is the system, where the initial proton is at rest, while the initial electron moves along the Z-axis. The transformation is performed via the sequence of auxiliary systems: Lab-->System 1-->System 2-->quasiLab(System 3). 
//As an output it gives the energy of the incoming electron in the quasiLab (E_beam_fermi), the four-momentum of the scattered electron in the quasiLab (P4_E_prime_boosted), and the rotation angle of the System 2-->quasiLab transformation (theta_rot2).

 void fermi_rot(Float_t &E_beam_fermi, Float_t &theta_rot2,Float_t E_beam,TLorentzVector  P4_E_prime, TLorentzVector  &P4_E_prime_boosted) {

 Float_t theta_fermi, phi_fermi, beta;
 TLorentzVector P4_EL, P4_in_Prot;


//The four-momenta of the incoming electron and the target proton in the Lab frame. 
 P4_EL.SetXYZT(0.,0.,E_beam,E_beam);
 P4_in_Prot.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 

//Transformation Lab --> System 1

theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;
 

P4_EL.RotateZ(-phi_fermi);
P4_EL.RotateY(-theta_fermi);

P4_E_prime.RotateZ(-phi_fermi);
P4_E_prime.RotateY(-theta_fermi);
 
P4_in_Prot.RotateZ(-phi_fermi);
P4_in_Prot.RotateY(-theta_fermi);
 

//Transformation System 1 --> System 2
 beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);

P4_in_Prot.Boost(0,0,-beta);
P4_EL.Boost(0,0,-beta);
P4_E_prime.Boost(0,0,-beta);


//Transformation System 2 --> quasiLab (Syatem 3)
theta_rot2 = acos(P4_EL[2]/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));


P4_EL.RotateY(theta_rot2);
P4_E_prime.RotateY(theta_rot2);
P4_in_Prot.RotateY(theta_rot2);

E_beam_fermi = P4_EL[3];
P4_E_prime_boosted=P4_E_prime;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//Test whether the target proton is at rest in the quasiLab.  If not --> smth is wrong!
 if ((fabs(P4_in_Prot[0])>0.001)||(fabs(P4_in_Prot[1])>0.001)||(fabs(P4_in_Prot[2])>0.001)||(fabs(P4_in_Prot[3]-MP)>0.0001))  cout << "ALARM! Wrong Lab-->quasiLab transformation! Proton is not at rest in the quasiLab! \n";

//Test whether the incoming electron moves along Z-axis in the quasiLab.  If not --> smth is wrong! 
 if ((fabs(P4_EL[0])>0.001)||(fabs(P4_EL[1])>0.001)||(fabs(P4_EL[2]-P4_EL[3])>0.001))  cout << "ALARM! Wrong Lab-->quasiLab transformation! Electron does not move along Z-axis in the quasiLab! \n";


};
