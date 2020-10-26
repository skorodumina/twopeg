#include "TMath.h"
#include <stdio.h>
#include <dlfcn.h>
#include <TGClient.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "TTree.h"
#include <sys/stat.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
 using namespace std;

int out_file_open() {

#ifdef BOS
//I. Create BOS output if needed
if ((flag_bos == 1)||(flag_bos == 2)){
remove(out_bos_file.c_str());
 
char mess[256];
initbos();
bankList( &bcs_, "E=", "HEADPARTMCTKMCVXMCEV" );
sprintf( mess, "OPEN BOSOUTPUT UNIT=1 FILE=\"%s\" WRITE STATUS=NEW RECL=3600", out_bos_file.c_str() );
     fparm_c( mess );
};//end of BOS output flag check
#endif     
  

//II. Create LUND output if needed   
if (flag_lund == 1){ 
remove(out_lund_file.c_str());
out_lund_stream.open(out_lund_file.c_str(), std::ofstream::out);
};//end of LUND output flag check
  
//III. Create the ROOT TREE with the cross section values by default
t21 = new TTree("h10","Tree h10"); 
t21->SetDirectory(0);
t21->Branch("sigma",&sigma_total); 
t21->Branch("p_el_test",&p_el_test); 
t21->Branch("px_fermi",&px_fermi); 
t21->Branch("py_fermi",&py_fermi);
t21->Branch("pz_fermi",&pz_fermi);

};
