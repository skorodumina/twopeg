#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"


using namespace std; 

void rot(Float_t Q2,Float_t E_beam, TLorentzVector P4_E_fin, TLorentzVector P4_1, TLorentzVector P4_2, TLorentzVector P4_3, Float_t &M12, Float_t &M23, Float_t &theta_hadr, Float_t &alpha_hadr, Float_t &phi_hadr){
//, Float_t &theta_hadr, Float_t &alpha_hadr, Float_t &phi_hadr 
//cout <<"qqq2 "<< P4_3[0]<< " " << P4_3[1]<< " " <<P4_3[2]<< " " <<P4_3[3]<< "\n ";
 Float_t a_gamma, b_gamma, a_beta,b_beta, beta;
// cout << E_beam<< "  " << MP<< "\n";
 // cout << "qqq2"<<P4_E_fin [0] << " "<< P4_E_fin [1] << " "<< P4_E_fin [2] <<" "<< P4_E_fin [3] << "\n";
TVector3 Vect3_gamma, Vect3_beta,V3_anti_z(0,0,-1);
TLorentzVector P4_E_ini, P4_gamma; 
 TRotation rot;

M12 = sqrt((P4_1+P4_2)*(P4_1+P4_2));
M23 = sqrt((P4_3+P4_2)*(P4_3+P4_2));
//inv_m31 = sqrt((P4_3+P4_1)*(P4_3+P4_1));


P4_E_ini.SetXYZT(0.,0.,E_beam, E_beam);

P4_gamma = P4_E_ini - P4_E_fin;

 TVector3 uz = P4_gamma.Vect().Unit();
 TVector3 ux = (P4_E_ini.Vect().Cross(P4_E_fin.Vect())).Unit();
 ux.Rotate(3.*M_PI/2,uz);
 rot.SetZAxis(uz,ux).Invert();
 P4_1.Transform(rot);
 P4_2.Transform(rot);
P4_3.Transform(rot);
 P4_gamma.Transform(rot);


 beta = sqrt(P4_gamma[3]*P4_gamma[3]+Q2)/(P4_gamma[3]+MP);
 
 P4_1.Boost(0,0,-beta);
P4_2.Boost(0,0,-beta);
P4_3.Boost(0,0,-beta);
P4_gamma.Boost(0,0,-beta);

theta_hadr = P4_1.Theta();
//theta_hadr = P4_2.Theta();
//theta_hadr = P4_3.Theta();
//cout <<"qqq2 "<< P4_1[0]<< " " << P4_1[1]<< " " <<P4_1[2]<< " " <<P4_1[3]<< "\n ";

if (P4_1.Phi()>0) phi_hadr = P4_1.Phi();
if (P4_1.Phi()<0) phi_hadr = P4_1.Phi()+2.*M_PI;

//if (P4_2.Phi()>0) phi_hadr = P4_2.Phi();
//if (P4_2.Phi()<0) phi_hadr = P4_2.Phi()+2.*M_PI;

//if (P4_3.Phi()>0) phi_hadr = P4_3.Phi();
//if (P4_3.Phi()<0) phi_hadr = P4_3.Phi()+2.*M_PI;

//a_gamma = sqrt(1./(1.-double(pow((P4_1.Vect().Unit() * V3_anti_z)),2)));
a_gamma = sqrt(1000000000000./(1000000000000.-(1000000.*(P4_1.Vect().Unit() * V3_anti_z))*(1000000.*(P4_1.Vect().Unit() * V3_anti_z))));
//cout << a_gamma<<" "<<(1000000.*(P4_1.Vect().Unit() * V3_anti_z))*(1000000.*(P4_1.Vect().Unit() * V3_anti_z))<<" p1 \n";
b_gamma = -(P4_1.Vect().Unit() * V3_anti_z)*a_gamma;

Vect3_gamma = a_gamma*V3_anti_z +b_gamma*P4_1.Vect().Unit();

a_beta = sqrt(1000000000000./(1000000000000.-(1000000.*(P4_1.Vect().Unit() * P4_2.Vect().Unit()))*(1000000.*(P4_1.Vect().Unit() * P4_2.Vect().Unit()))));

//cout << a_beta<<" "<<(1000000.*(P4_1.Vect().Unit() * P4_2.Vect().Unit()))*(1000000.*(P4_1.Vect().Unit() * P4_2.Vect().Unit())) <<" p3 \n";

b_beta = -(P4_1.Vect().Unit() * P4_2.Vect().Unit())*a_beta;

Vect3_beta = a_beta*P4_2.Vect().Unit() + b_beta*P4_1.Vect().Unit();

alpha_hadr = acos(Vect3_gamma * Vect3_beta);

if (Vect3_gamma.Cross(Vect3_beta) * P4_1.Vect() < 0.) alpha_hadr = 2.*M_PI - alpha_hadr;


/*
///2
a_gamma = sqrt(1./(1-pow((P4_3.Vect().Unit() * V3_anti_z),2)));
b_gamma = -(P4_3.Vect().Unit() * V3_anti_z)*a_gamma;
Vect3_gamma = a_gamma*V3_anti_z +b_gamma*P4_3.Vect().Unit();

a_beta = sqrt(1./(1-pow((P4_3.Vect().Unit() * P4_2.Vect().Unit()),2)));
b_beta = -(P4_3.Vect().Unit() * P4_2.Vect().Unit())*a_beta;
Vect3_beta = a_beta*P4_2.Vect().Unit() + b_beta*P4_3.Vect().Unit();

alpha_hadr = acos(Vect3_gamma * Vect3_beta);

if (Vect3_gamma.Cross(Vect3_beta) * P4_3.Vect() < 0) alpha_hadr = 2.*M_PI - alpha_hadr;
*/
/*
///3
a_gamma = sqrt(1./(1-pow((P4_2.Vect().Unit() * V3_anti_z),2)));
b_gamma = -(P4_2.Vect().Unit() * V3_anti_z)*a_gamma;
Vect3_gamma = a_gamma*V3_anti_z +b_gamma*P4_2.Vect().Unit();

a_beta = sqrt(1./(1-pow((P4_2.Vect().Unit() * P4_1.Vect().Unit()),2)));
b_beta = -(P4_2.Vect().Unit() * P4_1.Vect().Unit())*a_beta;
Vect3_beta = a_beta*P4_1.Vect().Unit() + b_beta*P4_2.Vect().Unit();

alpha_hadr = acos(Vect3_gamma * Vect3_beta);

if (Vect3_gamma.Cross(Vect3_beta) * P4_2.Vect() < 0) alpha_hadr = 2.*M_PI - alpha_hadr;
*/
 return;
 
 };
