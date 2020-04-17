#!/bin/bash

Env=$1
Srclist=$3

source $Env

touch $Srclist
cat > $Srclist << EOF
$CODEROOT/libs/shr/shr_orb_mod.F90
$CODEROOT/libs/shr/shr_sys_mod.F90
$CODEROOT/libs/shr/shr_mpi_mod.F90
$CODEROOT/libs/shr/shr_kind_mod.F90
$CODEROOT/libs/shr/shr_const_mod.F90
$CODEROOT/libs/c_coupler/src/Data_MGT
$CODEROOT/libs/c_coupler/src/Parallel_MGT
$CODEROOT/libs/c_coupler/src/Runtime_MGT
$CODEROOT/libs/c_coupler/src/Driver
$CODEROOT/libs/c_coupler/src/Utils
$CODEROOT/libs/c_coupler/src/XML
$CODEROOT/libs/c_coupler/src/CoR
EOF

exit 0
