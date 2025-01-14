#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# get fit range from user
RANGE_XMIN=$1
RANGE_XMAX=$2
NBINS_NO_B=$3
NBINS_B=$4

cd $SCRIPT_DIR/../../muon_lifetime
root -l -x "run_data_analysis.cc($RANGE_XMIN,$RANGE_XMAX,$NBINS_NO_B)"

cd $SCRIPT_DIR/..

root -l -x "run_data_analysis.cc($NBINS_B)"


