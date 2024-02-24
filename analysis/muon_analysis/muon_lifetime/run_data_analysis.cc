// Copyright 2023 nicolò salimbeni andrea de vita
/*
 */
#include <RtypesCore.h>

#include <iostream>
#include <string>

#include "TROOT.h"
#include "TSystem.h"

void run_data_analysis(int xmin_fit = 2000, int xmax_fit = 16000, int total_bins = 75)
{
  // load usefull libraries
  std::string lib_path     = "./../../../lib/";
  std::string include_path = "./../../../include/";
  std::string src_path     = "./../../../src/";

  gROOT->ProcessLine((".L " + src_path + "AnUtil.cc").c_str());
  gROOT->ProcessLine((".L " + src_path + "Event.cc").c_str());
  gROOT->ProcessLine((".L " + src_path + "InfoAcq.cc").c_str());
  gROOT->ProcessLine((".L " + src_path + "StateFile.cc").c_str());
  gROOT->ProcessLine((".L " + src_path + "GrUtil.cc").c_str());

  // load macros for this analysis
  gROOT->ProcessLine(".L ./src/muon_lifetime.cc");

  // create the command to launch the analysis
  std::string command = "FitData(\"./data/total.dat\",";
  command += std::to_string(xmin_fit);
  command += ",";
  command += std::to_string(xmax_fit);
  command += ",";
  command += std::to_string(total_bins);
  command += ")";
  std::cout << command << std::endl;
  gROOT->ProcessLine(command.c_str());
  gROOT->ProcessLine(".q");
}
