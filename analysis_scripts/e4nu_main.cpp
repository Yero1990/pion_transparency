#include "e4nu_analyzer.cpp"
#include "e4nu_analyzer.h"
#include <iostream>

int e4nu_main(TString inHIPO="", TString outROOT="", TString tgt="", TString det_h=""){

  e4nu_analyzer e4nu(inHIPO, outROOT, tgt, det_h);
  e4nu.run_data_analysis();
  
  return 0;
  
}
