#include <stdio.h>
#include <stdlib.h>  
#include "TROOT.h"
#include <RQ_OBJECT.h>
#include <iostream>
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include <string> 
#include <dlfcn.h>
#include <sstream>
#include <fstream>
#include "TTree.h"

#ifdef BOS
 //This is needed for BOS output creation
 extern "C" 
{     
#include <signal.h>
#include <errno.h>
#include <ntypes.h>
#include <bostypes.h>
#include <ec.h>
#include <clas_cern.h>
#include <ctype.h>
#include <kinematics.h>
#include <map_manager.h>
#include <trk.h>
#include <clasmdl.h>
#include <utility.h>
#include <pid.h>
#include <makebanks.h>
#include <call.h>
#include <bosddl.h>
#include <tagtnorm.h>
#include <vertex.h>


  void initbos(); 
  int  getBOS(BOSbank *bcs, int lun, char *list);
  void cleanBanks(BOSbank *bcs);
  void dropAllBanks(BOSbank *bcs, char *list);
  void *getBank(BOSbank *,const char *);
  void open_fpack_unit(char *filename,char *dataname,int unitnum);
  void close_fpack_unit(char *dataname);
  BOSbank bcs_ ;
  BOSbank wcs_ ;  
   
}
 //----
#endif  
    using namespace std; 

TTree *t21 ;
 
//variables from input file
Int_t Nevents;
Short_t Targ_Z, Targ_A;
Short_t  Twi_Z, Twf_Z, Twi_A, Twf_A;
Short_t flag_bos, flag_lund, flag_radmod,flag_fermi,flag_flux;


string out_bos_file, out_lund_file;
Float_t sigma_total,p_el_test;
Float_t W_min, W_max, Q2_min,Q2_max,Theta_min,Theta_max,E_eprime_min;
Float_t Targ_rad, Targ_len, Targ_off, Targ_dens, Targ_radlen;
Float_t  Twi_thick, Twf_thick,Twi_dens, Twf_dens,Twi_radlen, Twf_radlen;
Float_t E_beam;
//Particles' masses. They will be taken fom PDG-tables later in the main_prog 
Float_t MP, MPIP, MPIM,Me; 
Float_t px_fermi,py_fermi,pz_fermi;

//Variables from files with cross section
Float_t W_ARR[17];
Float_t Q2_ARR[3];
Float_t S12_ARR[12][17];
Float_t S23_ARR[12][17];
Float_t THETA_ARR[6]; 
Float_t ALPHA_ARR[6];
Float_t SIGMA_ARR[6][3][17][12][12][6][6];


Float_t W_ARR_GOL[30];
Float_t S12_ARR_GOL[16][30];
Float_t S23_ARR_GOL[16][30];
Float_t THETA_ARR_GOL[14]; 
Float_t ALPHA_ARR_GOL[14];
Float_t SIGMA_ARR_GOL[30][16][16][14][14];
 
Float_t W_ARR_FED[12];
Float_t Q2_ARR_FED[7];
Float_t S12_ARR_FED[10][12];
Float_t S23_ARR_FED[10][12];
Float_t THETA_ARR_FED[8]; 
Float_t ALPHA_ARR_FED[8];
Float_t SIGMA_ARR_FED[6][7][12][10][10][8][8];
 
Float_t EPSILON_L_RIPANI[3][17];
Float_t FIT_PARAM_SIGMA_T[4][17];
Float_t FIT_PARAM_SIGMA_L[3][17]; 
Float_t FIT_PARAM_SIGMA_C2F[3][17];  
Float_t FIT_PARAM_SIGMA_CF[3][17];    

Float_t FIT_PARAM_SIGMA_T_FED[2][10];
Float_t FIT_PARAM_SIGMA_L_FED[3][10]; 
Float_t FIT_PARAM_SIGMA_C2F_FED[3][10];  
Float_t FIT_PARAM_SIGMA_CF_FED[3][10]; 



Float_t S12_ARR_FED_THRESH[10][3];
Float_t S23_ARR_FED_THRESH[10][3];
Float_t W_ARR_FED_THRESH[3];
 
Float_t SIGMA_ARR_FED_THRESH[6][7][3][10][10][8][8];
  
Float_t W_ARR_RIP2[21];
Float_t S12_ARR_RIP2[12][21];
Float_t S23_ARR_RIP2[12][21];
Float_t SIGMA_ARR_RIP2[6][21][12][12][6][6];
Float_t EPS_L_RIP2[21];

Float_t W_ARR_RIP3[10];
Float_t S12_ARR_RIP3[16][10];
Float_t S23_ARR_RIP3[16][10]; 
Float_t SIGMA_ARR_RIP3[6][10][16][16][6][6];
 
Float_t W_ARR_gt_3[15];
Float_t S12_ARR_gt_3[16][15];
Float_t S23_ARR_gt_3[16][15];
 
Float_t SIGMA_ARR_gt_3[6][15][16][16][6][6];
Float_t SIGMA_ARR_phot_gt_3[15][16][16][6][6];
 
Float_t SIGMA_ARR_GOL2[10][16][16][6][6]; 

Float_t W_ARR_2pi_INT[71];
Float_t Q2_ARR_2pi_INT[27];
Float_t SIGMA_T_ARR_2pi_INT[27][71];
Float_t SIGMA_L_ARR_2pi_INT[27][71];
    
Float_t SIGMA_T_INT_RIPANI[17];
Float_t SIGMA_T_INT_GOLOVA[17]; 
//Cross-sections interpolated over 4 hadron variables (s12, s23, theta, alpha). They will be interpolated over W and Q2.
Float_t sigma_wright_q2left[6];
Float_t sigma_wright_q2right[6];
Float_t sigma_wleft_q2right[6];
Float_t sigma_wleft_q2left[6];
  
const char *data_dir; 
std::ofstream out_lund_stream;
std::ostringstream data_dir_2pi;
std::ostringstream PATH;

void global(){
//Starting values
px_fermi =0.;
py_fermi =0.;
pz_fermi =0.;
sigma_total = 0.;
p_el_test = 0.; 
  
};
