#ifndef GET_XSECT_RIPANI_H
#define GET_XSECT_RIPANI_H



Short_t getWbin (Float_t W);
Short_t getQ2bin (Float_t Q2);
Short_t getsbin (Short_t Wbin, Float_t sgen, Float_t Smax, Float_t Smin);
Short_t getanglebin (Float_t anglegen, Float_t anglemax);
void get_xsect_ripani(Float_t Q2gen, Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Float_t phigen,  Float_t &sigma_t_final, Float_t &sigma_l_final,Float_t  &sigma_c2f_final,Float_t  &sigma_s2f_final,Float_t &sigma_cf_final,Float_t  &sigma_sf_final);

#endif


