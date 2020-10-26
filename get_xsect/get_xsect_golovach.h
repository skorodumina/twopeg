#ifndef GET_XSECT_GOLOVACH_H
#define GET_XSECT_GOLOVACH_H



Short_t getsbin_GOL (Float_t sgen, Float_t Smax, Float_t Smin) ;

void get_xsect_golovach(Float_t Wgen, Float_t s12gen,Float_t s23gen, Float_t thetagen, Float_t alphagen, Short_t &Wleft_bin, Float_t &sigma_wright_gol,Float_t &sigma_wleft_gol);

#endif
