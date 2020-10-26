#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include <math.h>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TGObject.h"
#include "TObject.h"
#include "TSystem.h"
#include "TMinuit.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include "global.h"
#include <stdio.h>
#include <dlfcn.h>
#include <sstream>
#include <TLorentzVector.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring> 
#include <cstdlib>

 using namespace std;
 
 TH1F *h_W, *h_Q2, *h_phi_e, *h_zel,*h_nu;
 TH2F *h_Q2vsW, *h_Q2vsW2, *h_dalitz;

 TH1F *h_0_miss, *h_pim_miss, *h_0_miss_2, *h_pim_miss_2;
 TH1F *h_fermi_bonn, *h_0_miss_en, *h_0_miss_en_2,*h_eradgam;
 
 TH1F *h_e_beam_eff;
   
//Single-fold differentil yield for 1.225 < W < 2.15 GeV.
 TH1F *h_odn_inv_m12[37];
 TH1F *h_odn_inv_m23[37];
 TH1F *h_odn_alpha[37];
 TH1F *h_odn_theta[37]; 
 TH1F *h_odn_theta_2[37];
      
//Single-fold differentil yield for 2.1625 < W < 3.0375 GeV.   
 TH1F *h_odn_wwide_inv_m12[18];
 TH1F *h_odn_wwide_inv_m23[18];
 TH1F *h_odn_wwide_alpha[18];
 TH1F *h_odn_wwide_theta[18]; 
 TH1F *h_odn_wwide_theta_2[18];
 
//Single-fold differentil yield for 3.0875 < W < 4.5375 GeV.  
 TH1F *h_odn_wgt3_inv_m12[15];
 TH1F *h_odn_wgt3_inv_m23[15];
 TH1F *h_odn_wgt3_alpha[15];
 TH1F *h_odn_wgt3_theta[15]; 
 TH1F *h_odn_wgt3_theta_2[15];
          
//Integral w-yield for 0.1< Q2 < 1.3 GeV2. 
 TH1F *h_odn_w_dep_t[12];
 TH1F *h_odn_w_dep_l[12];
 TH1F *h_odn_w_dep_l2[12];
 TH1F *h_odn_w_dep_tot[12];
  
//Integral q2-yield for 1.25< Q2 < 2.075 GeV2.    
 TH1F *h_odn_q2_dep_t[33];
 TH1F *h_odn_q2_dep_l[33];
 TH1F *h_odn_q2_dep_l2[33];
 TH1F *h_odn_q2_dep_tot[33];
      
      
//     TH1F *h_int_crsect_t[27];
//     TH1F *h_int_crsect_l[27];
        
      
      
     int hist_def(Float_t E_beam){
     
Float_t phi_e_min = 0;
Float_t phi_e_max = 2*M_PI;
ostringstream qqq;

h_e_beam_eff = new TH1F("h_e_beam_eff","h_e_beam_eff",250,1.5,1.9);

h_W = new TH1F("W","W",100,W_min,W_max);
h_Q2 = new TH1F("Q2","Q2",100,Q2_min,Q2_max);
h_phi_e = new TH1F("phi_e","phi_e",100,phi_e_min,phi_e_max);
h_zel = new TH1F("Z_EL","Z_EL",100,Targ_off-Targ_len/2.-1.,Targ_off+Targ_len/2.+1.);
h_nu = new TH1F("nu","nu",100,-1.*E_beam,E_beam);

h_Q2vsW = new TH2F("Q2vsW","Q2vsWW",100,W_min,W_max,100,Q2_min,Q2_max);
h_Q2vsW2 = new TH2F("Q2vsW2","Q2vsWW2",69,W_min,W_max,4,Q2_min,Q2_max);
h_dalitz = new TH2F("dalitz","dalitz",100,MPIM+MPIP,W_max-MP,100,MPIP+MP,W_max-MPIM);

h_0_miss = new TH1F("h_0_miss","h_0_miss",1000,-0.05,0.05);
h_pim_miss = new TH1F("h_pim_miss","h_pim_miss",200,MPIM*MPIM-0.05,MPIM*MPIM+0.05);
h_0_miss_2 = new TH1F("h_0_miss_2","h_0_miss_2",1000,-0.05,0.05);
h_pim_miss_2 = new TH1F("h_pim_miss_2","h_pim_miss_2",200,MPIM*MPIM-0.05,MPIM*MPIM+0.05);
    
h_fermi_bonn = new TH1F("h_fermi_bonn","h_fermi_bonn",100,-0.1,0.9);
  
h_0_miss_en = new TH1F("h_0_miss_en","h_0_miss_en",200,-0.05,0.05);
h_0_miss_en_2 = new TH1F("h_0_miss_en_2","h_0_miss_en_2",200,-0.05,0.05);
  
h_eradgam = new TH1F("h_eradgam","h_eradgam",500,0.,0.6);

     
 //Single-fold differentil yield for 1.225 < W < 2.15 GeV.    
for (Short_t iy=0;iy<=36;iy++) {
qqq << "h_odn_inv_m12_" << 10000*(1.2375+0.025*iy);
h_odn_inv_m12[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MPIM ,1.25+0.025*iy-MP );
qqq.str("");

qqq << "h_odn_inv_m23_" << 10000*(1.2375+0.025*iy);
h_odn_inv_m23[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MP ,1.25+0.025*iy-MPIM );
qqq.str("");

qqq << "h_odn_alpha_" << 10000*(1.2375+0.025*iy);
h_odn_alpha[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,2*M_PI );
qqq.str("");

qqq << "h_odn_theta_" << 10000*(1.2375+0.025*iy);
h_odn_theta[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

qqq << "h_odn_theta_2_" << 10000*(1.2375+0.025*iy);
h_odn_theta_2[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

};


 //Single-fold differentil yield for 2.1625 < W < 3.0375 GeV.   
 for (Short_t iy=0; iy<18; iy++) {
qqq << "h_odn_wwide_inv_m12_" << 10000*(2.1875+0.05*iy);
h_odn_wwide_inv_m12[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MPIM ,2.2125+0.05*iy-MP );
qqq.str("");

qqq << "h_odn_wwide_inv_m23_" << 10000*(2.1875+0.05*iy);
h_odn_wwide_inv_m23[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MP ,2.2125+0.05*iy-MPIM );
qqq.str("");

qqq << "h_odn_wwide_alpha_" << 10000*(2.1875+0.05*iy);
h_odn_wwide_alpha[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,2*M_PI );
qqq.str("");

qqq << "h_odn_wwide_theta_" << 10000*(2.1875+0.05*iy);
h_odn_wwide_theta[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

qqq << "h_odn_wwide_theta_2_" << 10000*(2.1875+0.05*iy);
h_odn_wwide_theta_2[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");


};

 //Single-fold differentil yield for 3.0875 < W < 4.5375 GeV.
for (Short_t iy=0; iy<15; iy++) {
qqq << "h_odn_wgt3_inv_m12_" << 10000*(3.1375+0.1*iy);
h_odn_wgt3_inv_m12[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MPIM ,3.1875+0.1*iy-MP );
qqq.str("");

qqq << "h_odn_wgt3_inv_m23_" << 10000*(3.1375+0.1*iy);
h_odn_wgt3_inv_m23[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,MPIP+MP ,3.1875+0.1*iy-MPIM );
qqq.str("");

qqq << "h_odn_wgt3_alpha_" << 10000*(3.1375+0.1*iy);
h_odn_wgt3_alpha[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,2*M_PI );
qqq.str("");

qqq << "h_odn_wgt3_theta_" << 10000*(3.1375+0.1*iy);
h_odn_wgt3_theta[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

qqq << "h_odn_wgt3_theta_2_" << 10000*(3.1375+0.1*iy);
h_odn_wgt3_theta_2[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,0 ,M_PI );
qqq.str("");

};



       
//Integral w-yield for 0.1< Q2 < 1.3 GeV2.       
  for (Short_t iy=0; iy<=11; iy++) {

qqq << "h_odn_w_dep_t_" << 100*(0.15+0.1*iy);
h_odn_w_dep_t[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");

qqq << "h_odn_w_dep_l_" << 100*(0.15+0.1*iy);
h_odn_w_dep_l[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");

qqq << "h_odn_w_dep_l2_" << 100*(0.15+0.1*iy);
h_odn_w_dep_l2[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");

qqq << "h_odn_w_dep_tot_" << 100*(0.15+0.1*iy);
h_odn_w_dep_tot[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100, W_min, W_max);
qqq.str("");
     
      };
      
//Integral q2-yield for 1.25< Q2 < 2.075 GeV2.      
for (Short_t iy=0; iy<=32; iy++) { 
     
qqq << "h_odn_q2_dep_t_" << 10000*(1.2625+0.025*iy);
h_odn_q2_dep_t[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");

qqq << "h_odn_q2_dep_l_" << 10000*(1.2625+0.025*iy);
h_odn_q2_dep_l[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");

qqq << "h_odn_q2_dep_l2_" << 10000*(1.2625+0.025*iy);
h_odn_q2_dep_l2[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");


qqq << "h_odn_q2_dep_tot_" << 10000*(1.2625+0.025*iy);
h_odn_q2_dep_tot[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),100,Q2_min,Q2_max);
qqq.str("");
     };


/*
h_int_crsect_t[0] = new TH1F("h_int_crsect_t_0","h_int_crsect_t_0",71, 1.25, 3.025);   
h_int_crsect_l[0] = new TH1F("h_int_crsect_l_0","h_int_crsect_l_0",71, 1.25, 3.025);  

h_int_crsect_t[26] = new TH1F("h_int_crsect_t_26","h_int_crsect_t_26",71, 1.25, 3.025);   
h_int_crsect_l[26] = new TH1F("h_int_crsect_l_26","h_int_crsect_l_26",71, 1.25, 3.025);

for (Short_t iy=1; iy<=25; iy++) {
    
qqq << "h_int_crsect_t_" <<iy;
h_int_crsect_t[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),71, 1.25, 3.025);
qqq.str("");
    
qqq << "h_int_crsect_l_" <<iy;
h_int_crsect_l[iy] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),71, 1.25, 3.025);
qqq.str("");
    };*/  

};

	 
