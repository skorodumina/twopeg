#ifndef GLOBAL_H 
#define GLOBAL_H
#include "TTree.h"  

#ifdef BOS
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


   extern void initbos();
   extern int  getBOS(BOSbank *bcs, int lun, char *list);
   extern void cleanBanks(BOSbank *bcs);
   extern void dropAllBanks(BOSbank *bcs, char *list);
   extern void *getBank(BOSbank *,const char *);
   extern void open_fpack_unit(char *filename,char *dataname,int unitnum);
   extern void close_fpack_unit(char *dataname);
   extern   BOSbank bcs_ ;
   extern  BOSbank wcs_ ; 
}
#endif  
  using namespace std;  

extern TTree*t21 ;
      
extern Int_t Nevents;
extern Short_t Targ_Z, Targ_A;
extern Short_t  Twi_Z, Twf_Z, Twi_A, Twf_A;
extern Short_t flag_bos, flag_lund, flag_radmod,flag_fermi,flag_flux;

extern string out_bos_file, out_lund_file;
extern Float_t sigma_total,p_el_test;
extern Float_t W_min, W_max, Q2_min,Q2_max,Theta_min,Theta_max,E_eprime_min;
extern Float_t Targ_rad, Targ_len, Targ_off, Targ_dens, Targ_radlen;
extern Float_t  Twi_thick, Twf_thick,Twi_dens, Twf_dens,Twi_radlen, Twf_radlen;
  
extern Float_t MP, MPIP, MPIM, Me; 
extern Float_t px_fermi,py_fermi,pz_fermi;
 
 
extern Float_t W_ARR[17];
extern Float_t Q2_ARR[3];
extern Float_t S12_ARR[12][17];
extern Float_t S23_ARR[12][17];
extern Float_t THETA_ARR[6]; 
extern Float_t ALPHA_ARR[6];
extern Float_t SIGMA_ARR[6][3][17][12][12][6][6];
  
extern Float_t W_ARR_GOL[30];
extern Float_t S12_ARR_GOL[16][30];
extern Float_t S23_ARR_GOL[16][30];
extern Float_t THETA_ARR_GOL[14]; 
extern Float_t ALPHA_ARR_GOL[14];
extern Float_t SIGMA_ARR_GOL[30][16][16][14][14];


extern Float_t W_ARR_FED[12];
extern Float_t Q2_ARR_FED[7];
extern Float_t S12_ARR_FED[10][12];
extern Float_t S23_ARR_FED[10][12];
extern Float_t THETA_ARR_FED[8]; 
extern Float_t ALPHA_ARR_FED[8];
extern Float_t SIGMA_ARR_FED[6][7][12][10][10][8][8];

extern Float_t EPSILON_L_RIPANI[3][17];
extern Float_t FIT_PARAM_SIGMA_T[4][17];
extern Float_t FIT_PARAM_SIGMA_L[3][17]; 
extern Float_t FIT_PARAM_SIGMA_C2F[3][17];  
extern Float_t FIT_PARAM_SIGMA_CF[3][17];	
extern Float_t SIGMA_T_INT_RIPANI[17];
extern Float_t SIGMA_T_INT_GOLOVA[17];

 
extern Float_t FIT_PARAM_SIGMA_T_FED[2][10];
extern Float_t FIT_PARAM_SIGMA_L_FED[3][10]; 
extern Float_t FIT_PARAM_SIGMA_C2F_FED[3][10];  
extern Float_t FIT_PARAM_SIGMA_CF_FED[3][10]; 
 
extern Float_t S12_ARR_FED_THRESH[10][3];
extern Float_t S23_ARR_FED_THRESH[10][3];
 
extern Float_t W_ARR_RIP2[21];
extern Float_t S12_ARR_RIP2[12][21];
extern Float_t S23_ARR_RIP2[12][21];
extern Float_t SIGMA_ARR_RIP2[6][21][12][12][6][6];
extern Float_t EPS_L_RIP2[21];
 
extern Float_t W_ARR_RIP3[10];
extern Float_t S12_ARR_RIP3[16][10];
extern Float_t S23_ARR_RIP3[16][10]; 
extern Float_t SIGMA_ARR_RIP3[6][10][16][16][6][6];
    
extern Float_t W_ARR_gt_3[15];
extern Float_t S12_ARR_gt_3[16][15];
extern Float_t S23_ARR_gt_3[16][15];
     
extern Float_t SIGMA_ARR_gt_3[6][15][16][16][6][6];
extern Float_t SIGMA_ARR_phot_gt_3[15][16][16][6][6];
    
extern Float_t SIGMA_ARR_GOL2[10][16][16][6][6]; 
    
extern Float_t W_ARR_2pi_INT[71];
extern Float_t Q2_ARR_2pi_INT[27];
extern Float_t SIGMA_T_ARR_2pi_INT[27][71];
extern Float_t SIGMA_L_ARR_2pi_INT[27][71];
    
         
extern Float_t sigma_wright_q2left[6];
extern Float_t sigma_wright_q2right[6];
extern Float_t sigma_wleft_q2right[6];
extern Float_t sigma_wleft_q2left[6];

extern Float_t SIGMA_ARR_FED_THRESH[6][7][3][10][10][8][8];
extern Float_t W_ARR_FED_THRESH[3];  

extern const char *data_dir;
extern std::ofstream out_lund_stream; 
extern std::ostringstream data_dir_2pi;
extern std::ostringstream PATH;

void global();
#endif
