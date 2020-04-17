#!/bin/bash

#set -x
ROOT=$EXP_DIR
# =================================

#####################################################
# machine set up (users should change this part)
#####################################################
#

  #ANAL_TIME=2016060106
  ANAL_TIME=$1
  echo "run_time: " $ANAL_TIME
  YY=`echo $ANAL_TIME | cut -c1-4`
  MM=`echo $ANAL_TIME | cut -c5-6`
  DD=`echo $ANAL_TIME | cut -c7-8`
  HH=`echo $ANAL_TIME | cut -c9-10`
  DIR_ROOT=${ROOT}/run/ensemble_01/atm/wrf/data
  WORK_ROOT=${ROOT}/enkf
  GSI_ROOT=/home/sunchao/work/GSI_EnKF/lib_test/comGSIv3.6_EnKFv1.2/dtc
  ENKF_NAMELIST=${GSI_ROOT}/run/enkf_wrf_namelist.sh

# ensemble parameters
#
  NMEM_ENKF=$ENS_NUMS
  NLONS=480
  NLATS=360
  NLEVS=32
  IF_ARW=.true.
  IF_NMM=.false.

  # Build the GSI namelist on-the-fly
. $ENKF_NAMELIST

#  list="conv"
#  list="conv airs_aqua amsua_aqua amsua_metop-a amsua_metop-b amsua_n15 amsua_n18 amsua_n19 gome_metop-a gome_metop-b hirs4_metop-a hirs4_metop-b hirs4_n19 mhs_metop-a mhs_metop-b mhs_n18 mhs_n19 sbuv2_n19"
#

# get mean
memdir=${ROOT}/gsi_ensmean
list=`ls ${memdir}/pe* | cut -f2 -d"." | awk '{print substr($0, 0, length($0)-3)}' | sort | uniq `
rm ${memdir}/diag_*_ges.ensmean.${ANAL_TIME}
for type in $list; do
   count=`ls ${memdir}/pe*${type}_01* | wc -l`
   if [[ $count -gt 0 ]]; then
      cat ${memdir}/pe*${type}_01* > ${memdir}/diag_${type}_ges.ensmean.${ANAL_TIME}
      rm ${memdir}/pe*${type}_01*
   fi
   rm ${WORK_ROOT}/diag_${type}_ges.ensmean
   ln -s ${memdir}/diag_${type}_ges.ensmean.${ANAL_TIME} ${WORK_ROOT}/diag_${type}_ges.ensmean
   #ln -s ${memdir}/diag_${type}_ges.mem001.${ANAL_TIME} ${WORK_ROOT}/diag_${type}_ges.ensmean
   #cp ${memdir}/diag_${type}_ges.ensmean.${ANAL_TIME} ${WORK_ROOT}/diag_${type}_ges.ensmean
   #cp ${memdir}/diag_${type}_ges.mem001.${ANAL_TIME} ${WORK_ROOT}/diag_${type}_ges.ensmean
done

# get each member
imem=1
while [[ $imem -le $NMEM_ENKF ]]; do
   member="mem"`printf %03i $imem`
   memdir_id="ensemble_"`printf %02i $imem`
   memdir=${ROOT}/run/${memdir_id}/atm/wrf/data
   for type in $list; do
      if [ -r "${WORK_ROOT}/diag_${type}_ges.${member}" ]; then
          rm ${WORK_ROOT}/diag_${type}_ges.${member}
      fi
      if [ -r "${memdir}/diag_${type}_ges.${member}.${ANAL_TIME}" ]; then
          rm ${memdir}/diag_${type}_ges.${member}.${ANAL_TIME}
      fi
      count=`ls ${memdir}/pe*${type}_01* | wc -l`
      if [[ $count -gt 0 ]]; then
          cat ${memdir}/pe*${type}_01* > ${memdir}/diag_${type}_ges.${member}.${ANAL_TIME}
          ln -s ${memdir}/diag_${type}_ges.${member}.${ANAL_TIME} ${WORK_ROOT}/diag_${type}_ges.${member}
          rm ${memdir}/pe*${type}_01*
      fi
   done
   BK_FILE=${memdir}/wrfout_d01_${YY}-${MM}-${DD}_${HH}:00:00
   #BK_FILE=${memdir}/wrfinput_d01
   if [ -r "$BK_FILE" ]; then
       ncks -O -d Time,0 $BK_FILE  ${WORK_ROOT}/analysis.${member}
   fi
   #rm pe*.obs_setup
   #rm fort.*
   #rm obs_input.*
   (( imem = $imem + 1 ))
done
nces -O ${WORK_ROOT}/analysis.mem* ${WORK_ROOT}/analysis.ensmean
# Fixed files
# CONVINFO=${FIX_ROOT}/global_convinfo.txt
# SATINFO=${FIX_ROOT}/global_satinfo.txt
# SCANINFO=${FIX_ROOT}/global_scaninfo.txt
# OZINFO=${FIX_ROOT}/global_ozinfo.txt
ANAVINFO=${DIR_ROOT}/anavinfo
CONVINFO=${DIR_ROOT}/convinfo
SATINFO=${DIR_ROOT}/satinfo
SCANINFO=${DIR_ROOT}/scaninfo
OZINFO=${DIR_ROOT}/ozinfo
# LOCINFO=${FIX_ROOT}/global_hybens_locinfo.l64.txt

cp $ANAVINFO        ${WORK_ROOT}/anavinfo
cp $CONVINFO        ${WORK_ROOT}/convinfo
cp $SATINFO         ${WORK_ROOT}/satinfo
cp $SCANINFO        ${WORK_ROOT}/scaninfo
cp $OZINFO          ${WORK_ROOT}/ozinfo
# cp $LOCINFO         ${WORK_ROOT}/hybens_locinfo

cp $DIR_ROOT/satbias_in ${WORK_ROOT}/satbias_in
cp $DIR_ROOT/satbias_pc ${WORK_ROOT}/satbias_pc
