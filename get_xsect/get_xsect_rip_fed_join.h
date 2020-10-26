#ifndef GET_XSECT_RIP_FED_JOIN_H
#define GET_XSECT_RIP_FED_JOIN_H


Float_t pol1_fed_0575_065 (Float_t  x, Short_t wbin, Short_t i);
Float_t pol1_fed_065_095 (Float_t  x, Short_t wbin, Short_t i);
Float_t pol1_fed_095_130 (Float_t  x, Short_t wbin, Short_t i);
Float_t func_sigma_t_fed (Float_t  x, Short_t wbin);
Float_t pol2_fed (Float_t  x, Short_t wbin, Short_t i);





void get_xsect_rip_fed_join(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen, Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final );

#endif
