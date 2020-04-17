#!/bin/bash
set -x
ROOT=$EXP_DIR
# =================================

  #ANAL_TIME=2016060106
  ANAL_TIME=$1
  YY=`echo $ANAL_TIME | cut -c1-4`
  MM=`echo $ANAL_TIME | cut -c5-6`
  DD=`echo $ANAL_TIME | cut -c7-8`
  HH=`echo $ANAL_TIME | cut -c9-10`
  DIR_ROOT=/home/sunchao/work/GSI_EnKF
  OBS_ROOT=${DIR_ROOT}/GFS_data/obs/20160601-0607/CONV_GDAS/gdas.${YY}${MM}${DD}
  PREPBUFR=${OBS_ROOT}/gdas1.t${HH}z.prepbufr.nr
  GSI_ROOT=${DIR_ROOT}/lib_test/comGSIv3.6_EnKFv1.2/dtc
  FIX_ROOT=${GSI_ROOT}/../fix
  GSI_NAMELIST=${GSI_ROOT}/run/comgsi_namelist.sh
  ENSMEAN_ROOT=${ROOT}/gsi_ensmean

ensemble_num=$ENS_NUMS
bklevels=32
bklevels_stag=33

##
# if_observer = Yes  : only used as observation operater for enkf
# if_hybrid   = Yes  : Run GSI as 3D/4D EnVar
# if_4DEnVar  = Yes  : Run GSI as 4D EnVar
if_hybrid=No     # Yes, or, No -- case sensitive !
if_4DEnVar=No    # Yes, or, No -- case sensitive (if_hybrid must be Yes)!
if_observer=Yes   # Yes, or, No -- case sensitive !

bk_core=ARW
bkcv_option=NAM
if_clean=clean

BYTE_ORDER=Big_Endian
#BYTE_ORDER=Little_Endian
##
##
for ((i = 1; i <= ensemble_num ; i++))
do
    fix=2
    ens_index=$(printf "%0${fix}d" "$i")
    mem_index=$(printf "%03d" "$i")
    echo "ensemble index: " ${ens_index}
    BK_ROOT=$ROOT/run/ensemble_${ens_index}/atm/wrf/data
    cd ${BK_ROOT}
    echo `pwd`
    # Link to the prepbufr data
    if [ -r "prepbufr" ]; then
       rm prepbufr
    fi
    ln -s ${PREPBUFR} prepbufr
    echo "link source obs file ${PREPBUFR}"

    # ln -s ${OBS_ROOT}/gdas1.t${HH}z.sptrmm.tm00.bufr_d tmirrbufr
    # Link to the radiance data
    srcobsfile[1]=${OBS_ROOT}/gdas1.t${HH}z.satwnd.tm00.bufr_d
    gsiobsfile[1]=satwnd
    srcobsfile[2]=${OBS_ROOT}/gdas1.t${HH}z.1bamua.tm00.bufr_d
    gsiobsfile[2]=amsuabufr
    srcobsfile[3]=${OBS_ROOT}/gdas1.t${HH}z.1bhrs4.tm00.bufr_d
    gsiobsfile[3]=hirs4bufr
    srcobsfile[4]=${OBS_ROOT}/gdas1.t${HH}z.1bmhs.tm00.bufr_d
    gsiobsfile[4]=mhsbufr
    srcobsfile[5]=${OBS_ROOT}/gdas1.t${HH}z.1bamub.tm00.bufr_d
    gsiobsfile[5]=amsubbufr
    srcobsfile[6]=${OBS_ROOT}/gdas1.t${HH}z.ssmisu.tm00.bufr_d
    gsiobsfile[6]=ssmirrbufr
    srcobsfile[7]=${OBS_ROOT}/gdas1.t${HH}z.airsev.tm00.bufr_d
    gsiobsfile[7]=airsbufr
    srcobsfile[8]=${OBS_ROOT}/gdas1.t${HH}z.sevcsr.tm00.bufr_d
    gsiobsfile[8]=seviribufr
    srcobsfile[9]=${OBS_ROOT}/gdas1.t${HH}z.iasidb.tm00.bufr_d
    gsiobsfile[9]=iasibufr
    srcobsfile[10]=${OBS_ROOT}/gdas1.t${HH}z.gpsro.tm00.bufr_d
    gsiobsfile[10]=gpsrobufr
    srcobsfile[11]=${OBS_ROOT}/gdas1.t${HH}z.amsr2.tm00.bufr_d
    gsiobsfile[11]=amsrebufr
    srcobsfile[12]=${OBS_ROOT}/gdas1.t${HH}z.atms.tm00.bufr_d
    gsiobsfile[12]=atmsbufr
    srcobsfile[13]=${OBS_ROOT}/gdas1.t${HH}z.geoimr.tm00.bufr_d
    gsiobsfile[13]=gimgrbufr
    #srcobsfile[14]=${OBS_ROOT}/gdas1.t${HH}z.gome.tm00.bufr_d
    #gsiobsfile[14]=gomebufr
    srcobsfile[15]=${OBS_ROOT}/gdas1.t${HH}z.omi.tm00.bufr_d
    gsiobsfile[15]=omibufr
    srcobsfile[16]=${OBS_ROOT}/gdas1.t${HH}z.osbuv8.tm00.bufr_d
    gsiobsfile[16]=sbuvbufr
    srcobsfile[17]=${OBS_ROOT}/gdas1.t${HH}z.eshrs3.tm00.bufr_d
    gsiobsfile[17]=hirs3bufrears
    srcobsfile[18]=${OBS_ROOT}/gdas1.t${HH}z.esamua.tm00.bufr_d
    gsiobsfile[18]=amsuabufrears
    srcobsfile[19]=${OBS_ROOT}/gdas1.t${HH}z.esmhs.tm00.bufr_d
    gsiobsfile[19]=mhsbufrears
    srcobsfile[20]=${OBS_ROOT}/rap.t${HH}z.nexrad.tm00.bufr_d
    gsiobsfile[20]=l2rwbufr
    srcobsfile[21]=${OBS_ROOT}/rap.t${HH}z.lgycld.tm00.bufr_d
    gsiobsfile[21]=larcglb
    ii=1
    while [[ $ii -le 21 ]]; do
       #echo "link source obs file"
       if [ -r "${srcobsfile[$ii]}" ]; then
          rm ${gsiobsfile[$ii]}
       fi
       ln -s ${srcobsfile[$ii]}  ${gsiobsfile[$ii]}
       echo "link source obs file ${srcobsfile[$ii]}"
       (( ii = $ii + 1 ))
    done
    
    # default is NAM
    # as_op='1.0,1.0,0.5 ,0.7,0.7,0.5,1.0,1.0,'
    vs_op='1.0,'
    hzscl_op='0.373,0.746,1.50,'
    if [ ${bkcv_option} = GLOBAL ] ; then
    #   as_op='0.6,0.6,0.75,0.75,0.75,0.75,1.0,1.0'
       vs_op='0.7,'
       hzscl_op='1.7,0.8,0.5,'
    fi
    if [ ${bk_core} = NMMB ] ; then
       vs_op='0.6,'
    fi

    # default is NMM
       bk_core_arw='.false.'
       bk_core_nmm='.true.'
       bk_core_nmmb='.false.'
       bk_if_netcdf='.true.'
    if [ ${bk_core} = ARW ] ; then
       bk_core_arw='.true.'
       bk_core_nmm='.false.'
       bk_core_nmmb='.false.'
       bk_if_netcdf='.true.'
    fi
    if [ ${bk_core} = NMMB ] ; then
       bk_core_arw='.false.'
       bk_core_nmm='.false.'
       bk_core_nmmb='.true.'
       bk_if_netcdf='.false.'
    fi

    if [ ${if_observer} = Yes ] ; then
      nummiter=0
      if_read_obs_save='.true.'
      if_read_obs_skip='.false.'
    else
      nummiter=2
      if_read_obs_save='.false.'
      if_read_obs_skip='.false.'
    fi
    ifhyb=.false.
    if4d=.false.
    # Build the GSI namelist on-the-fly
    nummiter=0
    if_read_obs_save='.false.'
    if_read_obs_skip='.true.' 
    . $GSI_NAMELIST
    list=`ls ${ENSMEAN_ROOT}/obs_input.*`
    for type in $list; do
        rm obs_input.*
        echo "link obs_input file ${type}"
    done
    cp ${ENSMEAN_ROOT}/obs_input.* .

done
#


