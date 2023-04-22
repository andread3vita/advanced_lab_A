// Copyright 2023 nicolò salimbeni andrea de vita
/*
  this script is made to load properly all the libraries for
  the analysis
  TO DO
*/
#include <RtypesCore.h>

#include <string>

#include "TROOT.h"
#include "TSystem.h"

void run() {
  // load usefull libraries
  std::string lib_path     = "./../../lib/";
  std::string include_path = "./../../include/";

  gROOT->ProcessLine(".L ./../../src/AnUtil.cc");
  gROOT->ProcessLine(".L ./../../src/Event.cc");
  gROOT->ProcessLine(".L ./../../src/InfoAcq.cc");

  // load macros for this analysis
  gROOT->ProcessLine(".L ./src/Find_max.cc");
  gROOT->ProcessLine(".L ./src/ReadTree.C");
}
