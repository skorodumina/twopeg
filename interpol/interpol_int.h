#ifndef INTERPOL_INT_H
#define INTERPOL_INT_H


Float_t Func_q2_dep(Float_t Q2);

void interpol_int(Float_t Q2gen, Float_t Wgen,Float_t &sigma_t_int, Float_t &sigma_l_int);

void interpol_int_genev_old(Float_t Q2gen, Float_t Wgen,Float_t &sigma_t_int, Float_t &sigma_l_int);


#endif
