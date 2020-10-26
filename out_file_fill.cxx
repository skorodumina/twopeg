#include "TFile.h"
#include "TMath.h"
#include <stdio.h>
#include <dlfcn.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <TLorentzVector.h>
#include "global.h"
#include "hist_def.h"
#include <fstream>
#include <iomanip>
 using namespace std; 
 
//This file contains two subroutines: 
//out_file_fill(...) - it fills the desired output file with the generated info + the root-tree with weights.
//hist_fill(..) - it fills the histogram, which are needed for the illustrative and testing puposes.

void out_file_fill(int i,Float_t sigma, Float_t W, Float_t Q2, TLorentzVector &P4_E, TLorentzVector &P4_1, TLorentzVector &P4_2, TLorentzVector &P4_3, Float_t z_EL, Float_t x_EL, Float_t y_EL) {
 
 
#ifdef BOS
 
//%%%%%%%%%I. FILLING BOS OUTPUT WITH DIFFERENT BANKS IF CREATED%%%%%%%
if ((flag_bos == 1)||(flag_bos == 2)){
 
//Create HEAD-bank anyway
clasHEAD_t* HEAD ;
if ( ( HEAD = (clasHEAD_t*)makeBank( &bcs_, "HEAD", 0, sizeof(head_t)/4, 4 ) ) )
    {
      HEAD->head[0].version = 1; 
      HEAD->head[0].nevent = i ;
      HEAD->head[0].nrun = 1;
      HEAD->head[0].evtclass = 7 ;
      HEAD->head[0].type = -2 ;
      HEAD->head[0].time = 0;
      HEAD->head[0].trigbits = 1 ;
    }
    
clasMCEV_t* MCEV ;
if ( ( MCEV = (clasMCEV_t*)makeBank( &bcs_, "MCEV", 0, sizeof(mcev_t)/4, 1 ) ) )
    {
      MCEV->mcev[0].i1 = 100000*((float) rand() / (RAND_MAX)); 
      MCEV->mcev[0].i2 =100000*((float) rand() / (RAND_MAX)) ;
    
    } 
    
//-----------------------------------------------------------------------    
//Create MCTK & MCVX banks if needed  
if (flag_bos==1) {
    
clasMCTK_t* MCTK ;
if ( ( MCTK = (clasMCTK_t*)makeBank( &bcs_, "MCTK", 0, sizeof(mctk_t)/4, 4 ) ) )
    {
    
    //electron
      MCTK->mctk[0].cx = P4_E[0]/((P4_E.Vect()).Mag()); 
      MCTK->mctk[0].cy = P4_E[1]/((P4_E.Vect()).Mag()); 
      MCTK->mctk[0].cz = P4_E[2]/((P4_E.Vect()).Mag()); 
      MCTK->mctk[0].pmom = (P4_E.Vect()).Mag();
      MCTK->mctk[0].mass = 0.001;
      MCTK->mctk[0].charge = -1.;
      MCTK->mctk[0].id = 11;
      MCTK->mctk[0].flag = 0;
      MCTK->mctk[0].beg_vtx = 1;
      MCTK->mctk[0].end_vtx = 0;
      MCTK->mctk[0].parent = 0;
      
      //pi-
     
      MCTK->mctk[1].cx = P4_3[0]/((P4_3.Vect()).Mag()); 
      MCTK->mctk[1].cy = P4_3[1]/((P4_3.Vect()).Mag()); 
      MCTK->mctk[1].cz = P4_3[2]/((P4_3.Vect()).Mag()); 
      MCTK->mctk[1].pmom = (P4_3.Vect()).Mag();
      MCTK->mctk[1].mass = MPIM;
      MCTK->mctk[1].charge = -1.;
      MCTK->mctk[1].id = -211;
      MCTK->mctk[1].flag = 0;
      MCTK->mctk[1].beg_vtx = 1;
      MCTK->mctk[1].end_vtx = 0;
      MCTK->mctk[1].parent = 0;
      
      //pi+
     
      MCTK->mctk[2].cx = P4_2[0]/((P4_2.Vect()).Mag()); 
      MCTK->mctk[2].cy = P4_2[1]/((P4_2.Vect()).Mag()); 
      MCTK->mctk[2].cz = P4_2[2]/((P4_2.Vect()).Mag()); 
      MCTK->mctk[2].pmom = (P4_2.Vect()).Mag();
      MCTK->mctk[2].mass = MPIP;
      MCTK->mctk[2].charge = 1.;
      MCTK->mctk[2].id = 211;
      MCTK->mctk[2].flag = 0;
      MCTK->mctk[2].beg_vtx = 1;
      MCTK->mctk[2].end_vtx = 0;
      MCTK->mctk[2].parent = 0;
      
     //proton
     
      MCTK->mctk[3].cx = P4_1[0]/((P4_1.Vect()).Mag()); 
      MCTK->mctk[3].cy = P4_1[1]/((P4_1.Vect()).Mag()); 
      MCTK->mctk[3].cz = P4_1[2]/((P4_1.Vect()).Mag()); 
      MCTK->mctk[3].pmom = (P4_1.Vect()).Mag();
      MCTK->mctk[3].mass = MP;
      MCTK->mctk[3].charge = 1.;
      MCTK->mctk[3].id = 2212;
      MCTK->mctk[3].flag = 0;
      MCTK->mctk[3].beg_vtx = 1;
      MCTK->mctk[3].end_vtx = 0;
      MCTK->mctk[3].parent = 0;
      
      
      
    } //end MCTK bank creation 

clasMCVX_t* MCVX ;
if ( ( MCVX = (clasMCVX_t*)makeBank( &bcs_, "MCVX", 0, sizeof(mcvx_t)/4, 1 ) ) )
    {
      MCVX->mcvx[0].x = x_EL; 
      MCVX->mcvx[0].y = y_EL;
      MCVX->mcvx[0].z = z_EL;
      MCVX->mcvx[0].tof = 0.;
      MCVX->mcvx[0].flag = 0;
    }    //end MCVX bank creation   
};//end if flag_bos==1

//-------------------------------------------------

//Create PART bank if needed 
if (flag_bos==2) {

clasPART_t* PART ;
if ( ( PART = (clasPART_t*)makeBank( &bcs_, "PART", 0, sizeof(part_t)/4, 4) ) )
    {
    
    //electron
    PART->part[0].pid = 3.; 
    PART->part[0].vert.x = 0.;
    PART->part[0].vert.y = 0.;
    PART->part[0].vert.z = z_EL;
    PART->part[0].p.t = P4_E[3];
    PART->part[0].p.space.x = P4_E[0];
    PART->part[0].p.space.y = P4_E[1];
    PART->part[0].p.space.z =P4_E[2] ;
    PART->part[0].q = -1.;
    PART->part[0].trkid = 0.;
    PART->part[0].qpid = 0.;
    PART->part[0].qtrk = 0.;
    PART->part[0].flags = 0.;
       
       //proton
    PART->part[1].pid = 14.; 
    PART->part[1].vert.x = 0.;
    PART->part[1].vert.y = 0.;
    PART->part[1].vert.z = 0.;
    PART->part[1].p.t = P4_1[3];
    PART->part[1].p.space.x = P4_1[0];
    PART->part[1].p.space.y = P4_1[1];
    PART->part[1].p.space.z =P4_1[2] ;
    PART->part[1].q = 1.;
    PART->part[1].trkid = 0.;
    PART->part[1].qpid = 0.;
    PART->part[1].qtrk = 0.;
    PART->part[1].flags = 0.;   
       
           //pi+
    PART->part[2].pid = 8.; 
    PART->part[2].vert.x = 0.;
    PART->part[2].vert.y = 0.;
    PART->part[2].vert.z = 0.;
    PART->part[2].p.t = P4_2[3];
    PART->part[2].p.space.x = P4_2[0];
    PART->part[2].p.space.y = P4_2[1];
    PART->part[2].p.space.z =P4_2[2] ;
    PART->part[2].q = 1.;
    PART->part[2].trkid = 0.;
    PART->part[2].qpid = 0.;
    PART->part[2].qtrk = 0.;
    PART->part[2].flags = 0.;  
    
    
      //pi-
    PART->part[3].pid = 9.; 
    PART->part[3].vert.x = 0.;
    PART->part[3].vert.y = 0.;
    PART->part[3].vert.z = 0.;
    PART->part[3].p.t = P4_3[3];
    PART->part[3].p.space.x = P4_3[0];
    PART->part[3].p.space.y = P4_3[1];
    PART->part[3].p.space.z =P4_3[2] ;
    PART->part[3].q = -1.;
    PART->part[3].trkid = 0.;
    PART->part[3].qpid = 0.;
    PART->part[3].qtrk = 0.;
    PART->part[3].flags = 0.; 
            
    }  //end PART bank creation     
};//end if flag_bos==2

putBOS( &bcs_, 1, "E" );
dropAllBanks( &bcs_, "E");
cleanBanks( &bcs_ );
};//end if ((flag_bos == 1)||(flag_bos == 2))
#endif
      

//%%%%%%%%%II. FILLING LUND OUTPUT IF CREATED%%%%%%%%%%%%%%%%%%%      
if (flag_lund == 1){   
   
out_lund_stream <<"4  1.  1.  0  0  "<<i<< "  0  "<<std::fixed<<std::setprecision(6)<< W<< "  "<< Q2<< "  "<<sigma << "\n"; 
      
      //electron
       out_lund_stream <<" 1  -1.  1  11  0  0  "<<P4_E[0]<<"  ";
       out_lund_stream <<P4_E[1]<<"  "<<P4_E[2]<<"  "<<P4_E[3]<<"  ";
       out_lund_stream <<"0.0005  0.0000  0.0000  "<<z_EL<<"\n";
       
       //proton
       out_lund_stream <<" 2  1.  1  2212  0  0  "<<P4_1[0]<<"  ";
       out_lund_stream <<P4_1[1]<<"  "<<P4_1[2]<<"  "<<P4_1[3]<<"  ";
       out_lund_stream <<"0.9383  0.0000  0.0000  "<<z_EL<<"\n";
       
        //pi+
       out_lund_stream <<" 3  1.  1  211  0  0  "<<P4_2[0]<<"  ";
       out_lund_stream <<P4_2[1]<<"  "<<P4_2[2]<<"  "<<P4_2[3]<<"  ";
       out_lund_stream <<"0.1396  0.0000  0.0000  "<<z_EL<<"\n";
       
        //pi-
       out_lund_stream <<" 4  -1.  1  -211  0  0  "<<P4_3[0]<<"  ";
       out_lund_stream <<P4_3[1]<<"  "<<P4_3[2]<<"  "<<P4_3[3]<<"  ";
       out_lund_stream <<"0.1396  0.0000  0.0000  "<<z_EL<<"\n";
       
    };       
    
//III. Filling the root tree with the cross section weights and components of the fermi momentum
p_el_test = (P4_E.Vect()).Mag();
t21->Fill();

 };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//This subroutine fills the histogram, which are needed for illustrative and testing puposes. 

void hist_fill(Float_t E_beam,Float_t W,Float_t W_old, Float_t Q2, Float_t Q2_old, Float_t phi_e,Float_t z_EL,Float_t s12,Float_t s23, Float_t th_hadr,Float_t alph_hadr,Float_t ph_hadr,Float_t sigma_t_final,Float_t sigma_l_final,Float_t eps_l,Float_t sigma_total,TLorentzVector  P4_E_prime,TLorentzVector P4_Pfin,TLorentzVector  P4_PIP,TLorentzVector P4_PIM,TLorentzVector P4_Eini_new,TLorentzVector P4_E_prime_new){
    
Float_t xsect_int_test_t, xsect_int_test_l;
//if (W < 3.0125) interpol_int(Q2,W,xsect_int_test_t, xsect_int_test_l);

TLorentzVector P4_0_miss, P4_0_miss_2,P4_pim_miss,P4_pim_miss_2;
TLorentzVector P4_Pini,P4_Pini_fermi,P4_Eini;

P4_Pini.SetXYZT(0.,0.,0.,MP);
P4_Pini_fermi.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
P4_Eini.SetXYZT(0.,0.,E_beam,E_beam);

//P4_Pini - the four-momentum of the initial proton at rest 
//P4_Eini - the four-momentum of the incomimg electron in the Lab frame, no rad eff
//P4_E_prime - the four-momentum of the scattered electron in the Lab frame, no rad eff
//P4_PIP, P4_Pfin, P4_PIM - the four-momenta of the final hadrons in the Lab frame

//P4_Pini_fermi - the four-momentum of the moving initial proton

//P4_Eini_new - the four-momentum of the incomimg electron in the Lab frame, with rad eff
//P4_E_prime_new - the four-momentum of the scattered electron in the Lab frame, with rad eff


P4_0_miss = P4_Pini + P4_Eini - P4_E_prime - P4_PIP - P4_Pfin - P4_PIM;
  
if((flag_radmod == 0)&&(flag_fermi==1))   P4_0_miss_2 =P4_Pini_fermi +P4_Eini     -P4_E_prime     -P4_PIP-P4_Pfin-P4_PIM; 
if((!(flag_radmod == 0))&&(flag_fermi==0))P4_0_miss_2 =P4_Pini       +P4_Eini_new -P4_E_prime_new -P4_PIP-P4_Pfin-P4_PIM; 
if((!(flag_radmod == 0))&&(flag_fermi==1))P4_0_miss_2 =P4_Pini_fermi +P4_Eini_new -P4_E_prime_new -P4_PIP-P4_Pfin-P4_PIM; 


P4_pim_miss = P4_Pini + P4_Eini - P4_E_prime - P4_PIP - P4_Pfin;

if((flag_radmod == 0)&&(flag_fermi==1))    P4_pim_miss_2 = P4_Pini_fermi + P4_Eini     - P4_E_prime     - P4_PIP-P4_Pfin;  if((!(flag_radmod == 0))&&(flag_fermi==0)) P4_pim_miss_2 = P4_Pini       + P4_Eini_new - P4_E_prime_new - P4_PIP-P4_Pfin;  if((!(flag_radmod == 0))&&(flag_fermi==1)) P4_pim_miss_2 = P4_Pini_fermi + P4_Eini_new - P4_E_prime_new - P4_PIP-P4_Pfin;  


h_W->Fill(W_old,sigma_total);
h_Q2->Fill(Q2,sigma_total);
  
h_phi_e->Fill(phi_e,1.); 
h_zel->Fill(z_EL,1.);

h_Q2vsW->Fill(W_old,Q2_old,sigma_total);
h_Q2vsW2->Fill(W_old,Q2_old,1.);

h_nu->Fill((W*W+Q2-MP*MP)/2./MP,1.);
h_dalitz->Fill(sqrt(s12),sqrt(s23),1.);
  
h_0_miss->Fill(P4_0_miss.Mag2(),sigma_total);
h_0_miss_2->Fill(P4_0_miss_2.Mag2(),sigma_total);
 
h_pim_miss->Fill(P4_pim_miss.Mag2(),sigma_total);
h_pim_miss_2->Fill(P4_pim_miss_2.Mag2(),sigma_total);

h_0_miss_en->Fill(P4_0_miss[3],sigma_total);
h_0_miss_en_2->Fill(P4_0_miss_2[3],sigma_total);
h_fermi_bonn->Fill(P4_0_miss.Vect().Mag(),sigma_total);
  
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
//Single-fold differentil yield for 1.225 < W < 2.15 GeV.   
if ((W_old>=1.225)&&(W_old<=2.15)&&(Q2>=0.0001)&&(Q2<=1.4))  {
h_odn_inv_m12[int((W_old-1.225)/0.025)]->Fill(sqrt(s12),sigma_total);
h_odn_inv_m23[int((W_old-1.225)/0.025)]->Fill(sqrt(s23),sigma_total);
h_odn_alpha[int((W_old-1.225)/0.025)]->Fill(alph_hadr,sigma_total);
h_odn_theta[int((W_old-1.225)/0.025)]->Fill(th_hadr,sigma_total);
};

//Single-fold differentil yield for 2.1625 < W < 3.0375 GeV.  
if ((W>=2.1625)&&(W<=3.0375)&&(Q2>=0.0001)&&(Q2<=1.3))  {
h_odn_wwide_inv_m12[int((W-2.1625)/0.05)]->Fill(sqrt(s12),sigma_total);
h_odn_wwide_inv_m23[int((W-2.1625)/0.05)]->Fill(sqrt(s23),sigma_total);
h_odn_wwide_alpha[int((W-2.1625)/0.05)]->Fill(alph_hadr,sigma_total);
h_odn_wwide_theta[int((W-2.1625)/0.05)]->Fill(th_hadr,sigma_total);
};

//Single-fold differentil yield for 3.0875 < W < 4.5375 GeV.  
if ((W>=3.0875)&&(W<=4.5375)&&(Q2>=0.0001))  {
h_odn_wgt3_inv_m12[int((W-3.0875)/0.1)]->Fill(sqrt(s12),sigma_total);
h_odn_wgt3_inv_m23[int((W-3.0875)/0.1)]->Fill(sqrt(s23),sigma_total);
h_odn_wgt3_alpha[int((W-3.0875)/0.1)]->Fill(alph_hadr,sigma_total);
h_odn_wgt3_theta[int((W-3.0875)/0.1)]->Fill(th_hadr,sigma_total);
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
//Integral w-yield for 0.1< Q2 < 1.3 GeV2.   
if ((Q2>=0.1)&&(Q2<=1.3))  {
h_odn_w_dep_t[int((Q2-0.1)/0.1)]->Fill(W,sigma_t_final);
h_odn_w_dep_l[int((Q2-0.1)/0.1)]->Fill(W,eps_l*sigma_l_final);
h_odn_w_dep_l2[int((Q2-0.1)/0.1)]->Fill(W,sigma_l_final);
h_odn_w_dep_tot[int((Q2-0.1)/0.1)]->Fill(W,sigma_total); 
};
   
//Integral q2-yield for 1.25< Q2 < 2.075 GeV2.    
if ((W>1.25)&&(W<=2.075))  { 
h_odn_q2_dep_t[int((W-1.25)/0.025)]->Fill(Q2,sigma_t_final);
h_odn_q2_dep_l[int((W-1.25)/0.025)]->Fill(Q2,eps_l*sigma_l_final);
h_odn_q2_dep_l2[int((W-1.25)/0.025)]->Fill(Q2,sigma_l_final);
h_odn_q2_dep_tot[int((W-1.25)/0.025)]->Fill(Q2,sigma_total);
};
     
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
/*
if ((Q2>=0.00049)&&(Q2<=0.00051))  {
h_int_crsect_t[0]->Fill(W,xsect_int_test_t);
h_int_crsect_l[0]->Fill(W,xsect_int_test_l);
};
          
if ((Q2>=0.0005)&&(Q2<=1.2505))  {
h_int_crsect_t[int((Q2-0.0005)/0.05)+1]->Fill(W,xsect_int_test_t);
h_int_crsect_l[int((Q2-0.0005)/0.05)+1]->Fill(W,xsect_int_test_l);
};
          
if ((Q2>=1.29)&&(Q2<=1.3))  {
h_int_crsect_t[26]->Fill(W,xsect_int_test_t);
h_int_crsect_l[26]->Fill(W,xsect_int_test_l);
};  */ 

};


 
