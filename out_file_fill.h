#ifndef OUT_FILE_WRITE_H 
#define OUT_FILE_WRITE_H


void out_file_fill(int i, Float_t sigma, Float_t W, Float_t Q2, TLorentzVector &P4_E,TLorentzVector &P4_1, TLorentzVector &P4_2, TLorentzVector &P4_3,Float_t z_EL, Float_t x_EL,  Float_t y_EL);

void hist_fill(Float_t E_beam, Float_t W,Float_t W_old, Float_t Q2,Float_t Q2_old, Float_t phi_e,Float_t z_EL,Float_t s12,Float_t s23, Float_t th_hadr,Float_t alph_hadr,Float_t ph_hadr,Float_t sigma_t_final,Float_t sigma_l_final,Float_t eps_l,Float_t sigma_total,TLorentzVector  P4_E_prime,TLorentzVector P4_Pfin,TLorentzVector  P4_PIP,TLorentzVector P4_PIM,TLorentzVector P4_Eini_new,TLorentzVector P4_E_prime_new);

#endif

