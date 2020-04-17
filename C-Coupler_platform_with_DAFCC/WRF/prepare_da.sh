#!/bin/bash
source source_da_env.sh
echo
echo
echo "************* Prepare for GSI ensmean run ************"
echo
./prepare_gsi_ensmean_online_run.sh
echo
echo "************* Prepare for GSI members run ************"
echo
./prepare_gsi_online_run.sh
echo
echo "************* Finish preparation for GSI run ************"

