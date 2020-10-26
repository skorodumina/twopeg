#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"


using namespace std;
//Byckling function declaration
     Float_t G_BYCKLING1(Float_t x, Float_t y, Float_t z, Float_t u, Float_t v, Float_t w) {
     return x*x*y+x*y*y+z*z*u+z*u*u+v*v*w+v*w*w+x*z*w+x*u*v+y*z*v+y*u*w-x*y*(z+u+v+w)-z*u*(x+y+v+w)-v*w*(x+y+z+u);
     };
//----------------------------------------------------
//WE ASSUMED THAT IN REACTION A + B = 1 + 2 + 3
//A - VIRTUAL PHOTON
//B - INITIAL PROTON
//1 - PARTICLE WITH MASS m1
//2 - PARTICLE WITH MASS m2
//3 - PARTICLE WITH MASS m3

//We assume as global variables: MP - initial proton mass (is taken from pdg-table), E_beam - beam energy (is taken from input file).

//Here all derivations are performed via special root functions. If you need further level of explanation please see anti_rot_explanation, where all calculations are shown in more details and rotations are performed via matrices.

 void anti_rot(Float_t W, Float_t Q2, Float_t phi_el, Float_t E_beam, Float_t M12, Float_t M23, Float_t theta_hadr, Float_t alpha_hadr, Float_t phi_hadr, Float_t m1, Float_t m2, Float_t m3, TLorentzVector &P4_1, TLorentzVector &P4_2, TLorentzVector &P4_3) {




//Energy of virtual photon is CMS-frame
Float_t E_gamma = (W*W-Q2-MP*MP)/2./W;
//cout << MP << "\n";

 //4-momenta of initial particles in CMS-frame
TLorentzVector P4_gamma,P4_P_ini;
P4_gamma.SetXYZT(0,0,sqrt(E_gamma*E_gamma+Q2),E_gamma); 
P4_P_ini.SetXYZT(0,0,-P4_gamma[2],sqrt(MP*MP+P4_gamma[2]*P4_gamma[2]));


// One of the invariant masses should be derived from the other two (according to Byckling-Kajantie)
Float_t M13 = sqrt(W*W+m1*m1+m2*m2+m3*m3-M12*M12-M23*M23);

//cout << m1<< "  "<< m2 << " "<< m3<< "\n";
//Energies and magnitudes of 3-momenta of final particles in CMS-system
// Float_t en_1,en_2,en_3,mag_1, mag_2, mag_3, th_1, th_2, th_3, ph_1, ph_2, ph_3;
//Set energy components with arbitrary spatial components

 P4_1.SetXYZT(1.,1.,1.,(W*W+m1*m1-M23*M23)/2./W);
 P4_2.SetXYZT(1.,1.,1.,(W*W+m2*m2-M13*M13)/2./W);
 P4_3.SetXYZT(1.,1.,1.,(W*W+m3*m3-M12*M12)/2./W);
// en_1 = (W*W+m1*m1-M23*M23)/2./W;
// en_2 = (W*W+m2*m2-M13*M13)/2./W;
// en_3 = (W*W+m3*m3-M12*M12)/2./W;
 
 
// cout << en_3 << "qqq \n";
// cout << E_beam + MP - en_1 -en_2<< "qqq1 \n";
 
 //Set spatial magnitudes = magnitudes of 3-momenta

 P4_1.SetRho(sqrt(P4_1[3]*P4_1[3]-m1*m1));
 P4_2.SetRho(sqrt(P4_2[3]*P4_2[3]-m2*m2));
 P4_3.SetRho(sqrt(P4_3[3]*P4_3[3]-m3*m3));
 
 if ((P4_1[3]*P4_1[3]-m1*m1<0.)||(P4_2[3]*P4_2[3]-m2*m2<0.)||(P4_3[3]*P4_3[3]-m3*m3<0.)) cout << " mag \n";
 
//mag_1 = sqrt(en_1*en_1-m1*m1);
//mag_2 = sqrt(en_2*en_2-m2*m2);
//mag_3 = sqrt(en_3*en_3-m3*m3);

 
 //At this stage only magnitudes of 3-momenta and energy components of 4-momenta of final hadrons are defined (they are the same for ANY axes orientations). To define 4-momenta completely we need to calculate their spatial angles for SOME axis orientation.
 
 //Angles calculation in the system where Z-axis along particle #1, X perpendicular to Z and situated in a-b-1 plane (plane A). 
 //Two remaining final hadrons (#2 & #3) are located in plane B.
 //We choose such a frame, because in this frame azimuthal angles of remaining final hadrons (phi-angles) are equal (or derived via) given alpha-angle that is defined as angle between planes A and B.
 //Here we simultaneously set defined angles (theta and phi) as angles of 4-momenta of final hadrons that means we completely define them.
  
//For final hadron #1


 P4_1.SetTheta(0.);
 P4_1.SetPhi(0.);

//th_1 = 0.;
//ph_1 = 0.; 

//For final hadron #2

 if (!(isnan(acos((m1*m1+m2*m2+2*P4_1[3]*P4_2[3]-M12*M12)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_2[3]*P4_2[3]-m2*m2)))))) P4_2.SetTheta(acos((m1*m1+m2*m2+2*P4_1[3]*P4_2[3]-M12*M12)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_2[3]*P4_2[3]-m2*m2))));
  
 if ((isnan(acos((m1*m1+m2*m2+2*P4_1[3]*P4_2[3]-M12*M12)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_2[3]*P4_2[3]-m2*m2)))))&&((m1*m1+m2*m2+2*P4_1[3]*P4_2[3]-M12*M12)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_2[3]*P4_2[3]-m2*m2))>0.)) P4_2.SetTheta(0.);
 if ((isnan(acos((m1*m1+m2*m2+2*P4_1[3]*P4_2[3]-M12*M12)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_2[3]*P4_2[3]-m2*m2)))))&&((m1*m1+m2*m2+2*P4_1[3]*P4_2[3]-M12*M12)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_2[3]*P4_2[3]-m2*m2))<0.)) P4_2.SetTheta(3.14159); 
// cout << acos(1.)<< " " << " " << acos(-1.)<<"\n";
 

 
 P4_2.SetPhi(alpha_hadr);
 
//  th_2 = acos((m1*m1+m2*m2+2*en_1*en_2-M12*M12)/2./mag_1/mag_2);
//ph_2 = alpha_hadr;

//if (isnan(th_2)&&((m1*m1+m2*m2+2*P4_1[3]*P4_2[3]-M12*M12)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_2[3]*P4_2[3]-m2*m2))>0.)) th_2 = 0.;
// if (isnan(th_2)&&((m1*m1+m2*m2+2*P4_1[3]*P4_2[3]-M12*M12)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_2[3]*P4_2[3]-m2*m2))<0.)) th_2 = M_PI; 
 
//For final hadron #3 
  
 if (!(isnan(acos((m1*m1+m3*m3+2*P4_1[3]*P4_3[3]-M13*M13)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_3[3]*P4_3[3]-m3*m3)))))) P4_3.SetTheta(acos((m1*m1+m3*m3+2*P4_1[3]*P4_3[3]-M13*M13)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_3[3]*P4_3[3]-m3*m3))));
 
 if ((isnan(acos((m1*m1+m3*m3+2*P4_1[3]*P4_3[3]-M13*M13)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_3[3]*P4_3[3]-m3*m3)))))&&((m1*m1+m3*m3+2*P4_1[3]*P4_3[3]-M13*M13)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_3[3]*P4_3[3]-m3*m3))>0.)) P4_3.SetTheta(0.);
 if ((isnan(acos((m1*m1+m3*m3+2*P4_1[3]*P4_3[3]-M13*M13)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_3[3]*P4_3[3]-m3*m3)))))&&((m1*m1+m3*m3+2*P4_1[3]*P4_3[3]-M13*M13)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_3[3]*P4_3[3]-m3*m3))<0.)) P4_3.SetTheta(3.14159); 
  

 if (alpha_hadr > M_PI) P4_3.SetPhi(alpha_hadr - M_PI); 
 if (alpha_hadr <= M_PI) P4_3.SetPhi(alpha_hadr + M_PI);
 
 
/* th_3 = acos((m1*m1+m3*m3+2*en_1*en_3-M13*M13)/2./mag_1/mag_3);
 
 if ((isnan(th_3))&&((m1*m1+m3*m3+2*P4_1[3]*P4_3[3]-M13*M13)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_3[3]*P4_3[3]-m3*m3)))>0.)  th_3 = 0.; 
  if  ((isnan(th_3))&&((m1*m1+m3*m3+2*P4_1[3]*P4_3[3]-M13*M13)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_3[3]*P4_3[3]-m3*m3)))<0.) th_3=M_PI;
 
if (alpha_hadr > M_PI) ph_3 = alpha_hadr - M_PI;
if (alpha_hadr <= M_PI) ph_3 = alpha_hadr + M_PI;

P4_1.SetXYZT(mag_1*cos(ph_1)*sin(th_1),mag_1*sin(ph_1)*sin(th_1),mag_1*cos(th_1),en_1);
P4_2.SetXYZT(mag_2*cos(ph_2)*sin(th_2),mag_2*sin(ph_2)*sin(th_2),mag_2*cos(th_2),en_2);
P4_3.SetXYZT(mag_3*cos(ph_3)*sin(th_3),mag_3*sin(ph_3)*sin(th_3),mag_3*cos(th_3),en_3);
*/

//cout << P4_3 [0] << "\n";
//Now we need to rotate our auxiliary coordinate system to the usual CMS-axes-orientation: Z-axis along particle a (virtual photon) and X-axis is in ee'-plane
//For that purpose we need to do two rotations of coordinate axes:

//1 - rotate Z-axis with the angle theta_hadr in the XZ-plane (around Y-axis) to set  Z-axis along particle a (virtual photon). After this rotation X-axis is still in hadronic plane.
//2 - rotate X-axis with the angle phi_hardr in XY-plane (around Z-axis) to set X-axis is in ee'-plane



if ((isnan(P4_1.Mag()))||(isnan(P4_2.Mag()))|| (isnan(P4_3.Mag()))) {

cout << P4_1[0] << " "<< P4_1[1] << " "<< P4_1[2] << " "<< P4_1[3] << " t\n";
cout << P4_2[0] << " "<< P4_2[1] << " "<< P4_2[2] << " "<< P4_2[3] << " t\n";
cout << P4_3[0] << " "<< P4_3[1] << " "<< P4_3[2] << " "<< P4_3[3] << " t\n";


cout <<" "<<P4_1.Mag()<<" "<< P4_2.Mag()<<" " << P4_3.Mag()<<" " <<"   "<<G_BYCKLING1(M12*M12,M23*M23,W*W,m2*m2,m1*m1,m3*m3)<<" "<< (m1*m1+m2*m2+2*P4_1[3]*P4_2[3]-M12*M12)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_2[3]*P4_2[3]-m2*m2)) <<" "<< (m1*m1+m3*m3+2*P4_1[3]*P4_3[3]-M13*M13)/2./(sqrt(P4_1[3]*P4_1[3]-m1*m1))/(sqrt(P4_3[3]*P4_3[3]-m3*m3)) <<" \n"; 
};
//Here we transfom our final hadrons 4-momenta into the system described above
  P4_1.RotateY(theta_hadr);
  P4_1.RotateZ(phi_hadr); 
  P4_2.RotateY(theta_hadr);
  P4_2.RotateZ(phi_hadr); 
  P4_3.SetVect(-1.*P4_1.Vect()-P4_2.Vect());
  
//After that our final hadrons 4-momenta are defined in the usual CMS-system that that will be in the LAB frame after the boost
 
//cout <<"qqq1 "<< P4_1[0]<< " " << P4_1[1]<< " " <<P4_1[2]<< " " <<P4_1[3]<< "\n "; 
  
//Now we need to perform boost to find ourselves in the LAB-frame.
//beta here is beta of transformation from LAB-system to CMS and since root function requires to use -beta we have no minus sign in the formulae below

//Energy of virtual photon is LAB-frame needed to calculate beta.
Float_t E_gamma_lab = (W*W+Q2-MP*MP)/2./MP; 

P4_1.Boost(0,0,sqrt(E_gamma_lab*E_gamma_lab+Q2)/(E_gamma_lab+MP));
P4_2.Boost(0,0,sqrt(E_gamma_lab*E_gamma_lab+Q2)/(E_gamma_lab+MP));
P4_3.Boost(0,0,sqrt(E_gamma_lab*E_gamma_lab+Q2)/(E_gamma_lab+MP));


//Now we are in the LAB-frame (let's mark it frame 2) with Z2 - along the virtual photon, X2 - in the ee'-plane and Y2 - perpendicular to the e'-plane.
//The only thing left is to rotate our coordinate axes into usual LAB-frame-axes-orientation (lets call it frame 1): Z1-axis along the beam, Y1 - up, X1 - to the south.
//For that purpose we will find coordinates of unit vectors along the axes of frame 2 in frame 1.

  
 TVector3 P3_EL,P3_G,P3_X,P3_Y;
 
 //We define 3-vector of scattered electron in LAB-frame 1 (theta is from the formula for Q2 in the LAB frame: Q2 = 2E_beamE_e'(1-cos(theta_e')))
 P3_EL.SetXYZ(0.,0.,1.);
 P3_EL.SetTheta(acos(1.- Q2/E_beam/(E_beam-E_gamma_lab)/2.));
// cout << acos(1.- Q2/E_beam/(E_beam-E_gamma_lab)/2.) <<" rot\n";
 //cout << (E_beam-E_gamma_lab) << " qqq \n";
 P3_EL.SetPhi(phi_el); 
 
 //And we define 3-vector of virtual photon in LAB-frame 1
 P3_G.SetXYZ(0.,0.,1.);
P3_G.SetTheta(acos((Q2+2.*E_beam*E_gamma_lab)/2./E_beam/(sqrt(Q2+E_gamma_lab*E_gamma_lab))));
 if (phi_el < M_PI) P3_G.SetPhi(phi_el+ M_PI);
 if (phi_el >= M_PI) P3_G.SetPhi(phi_el- M_PI);
 
 
 //P3_G is along Z2 and now it is written in the frame 1
 
 //We define P3_Y that is along Y2 (perpendicular to the scattering plane)
 P3_Y = (P3_G.Cross(P3_EL)).Unit();
 //And P3_X that is along X2 (in scattering plane)
 P3_X = (P3_Y.Cross(P3_G)).Unit(); 
 
 //Now unit vectors of frame 2 are written in the frame 1. We define rotation of frame 2 -> frame 1 via the Euler Angles
  


  TRotation vrot;
  vrot.SetXEulerAngles(atan2(P3_X[2],P3_Y[2]),acos(P3_G[2]),atan2(P3_G[0],-1.*P3_G[1]));
  
  //And transfom our 4-momenta of final hadrons according to axes transformation
 P4_1.Transform(vrot);
 P4_2.Transform(vrot);
 P4_3.Transform(vrot);
 
 //cout <<"qqq1 "<< P4_3[0]<< " " << P4_3[1]<< " " <<P4_3[2]<< " " <<P4_3[3]<< "\n ";
 

    
 return;
 
 };
