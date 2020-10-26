#ifndef GET_XSECT_14_18_LOWQ2_FIT_H
#define GET_XSECT_14_18_LOWQ2_FIT_H

Float_t getEpsL (Float_t E_beam, Float_t W, Float_t Q2);
Float_t getEpsT (Float_t E_beam, Float_t W, Float_t Q2);
Float_t func_sigma_t (Float_t  x, Short_t wbin) ;
Float_t pol2 (Float_t  x, Short_t wbin, Short_t i);
void get_xsect_14_18_lowq2_fit(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final );
#endif


