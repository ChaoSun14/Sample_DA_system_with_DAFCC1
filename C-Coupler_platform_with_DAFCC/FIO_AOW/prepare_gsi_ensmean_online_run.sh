#!/bin/bash
#set -x
ROOT=$EXP_DIR
# =================================
DIR_ROOT=/home/sunchao/work/GSI_EnKF
CRTM_ROOT=${DIR_ROOT}/CRTM_2.2.3
GSI_ROOT=${DIR_ROOT}/lib_test/comGSIv3.6_EnKFv1.2/dtc
FIX_ROOT=${GSI_ROOT}/../fix
GSI_NAMELIST=${GSI_ROOT}/run/comgsi_namelist.sh

bklevels=70
bklevels_stag=71

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
echo ${ROOT}
    BK_ROOT=$ROOT/gsi_ensmean
    echo " Create working directory:" ${BK_ROOT}
    if [ ! -d "${BK_ROOT}" ]; then
          mkdir -p ${BK_ROOT}
    fi
    #cd ${BK_ROOT}
    rm ${BK_ROOT}/pe*_01*
    rm ${BK_ROOT}/*obs_input*
    rm ${BK_ROOT}/*obs_setup*
    echo " Copy fixed files and link CRTM coefficient files to working directory"

    # Set fixed files
    #   berror   = forecast model background error statistics
    #   specoef  = CRTM spectral coefficients
    #   trncoef  = CRTM transmittance coefficients
    #   emiscoef = CRTM coefficients for IR sea surface emissivity model
    #   aerocoef = CRTM coefficients for aerosol effects
    #   cldcoef  = CRTM coefficients for cloud effects
    #   satinfo  = text file with information about assimilation of brightness temperatures
    #   satangl  = angle dependent bias correction file (fixed in time)
    #   pcpinfo  = text file with information about assimilation of prepcipitation rates
    #   ozinfo   = text file with information about assimilation of ozone data
    #   errtable = text file with obs error for conventional data (regional only)
    #   convinfo = text file with information about assimilation of conventional data
    #   bufrtable= text file ONLY needed for single obs test (oneobstest=.true.)
    #   bftab_sst= bufr table for sst ONLY needed for sst retrieval (retrieval=.true.)

    if [ ${bkcv_option} = GLOBAL ] ; then
      echo ' Use global background error covariance'
      BERROR=${FIX_ROOT}/${BYTE_ORDER}/nam_glb_berror.f77.gcv
      OBERROR=${FIX_ROOT}/prepobs_errtable.global
      if [ ${bk_core} = NMM ] ; then
         ANAVINFO=${FIX_ROOT}/anavinfo_ndas_netcdf_glbe
      fi
      if [ ${bk_core} = ARW ] ; then
        ANAVINFO=${FIX_ROOT}/anavinfo_arw_netcdf_glbe
      fi
      if [ ${bk_core} = NMMB ] ; then
        ANAVINFO=${FIX_ROOT}/anavinfo_nems_nmmb_glb
      fi
    else
      echo ' Use NAM background error covariance'
      BERROR=${FIX_ROOT}/${BYTE_ORDER}/nam_nmmstat_na.gcv
      OBERROR=${FIX_ROOT}/nam_errtable.r3dv
      if [ ${bk_core} = NMM ] ; then
         ANAVINFO=${FIX_ROOT}/anavinfo_ndas_netcdf
      fi
      if [ ${bk_core} = ARW ] ; then
         ANAVINFO=${FIX_ROOT}/anavinfo_arw_netcdf
      fi
      if [ ${bk_core} = NMMB ] ; then
         ANAVINFO=${FIX_ROOT}/anavinfo_nems_nmmb
      fi
    fi

    SATANGL=${FIX_ROOT}/global_satangbias.txt
    SATINFO=${FIX_ROOT}/global_satinfo.txt
    CONVINFO=${FIX_ROOT}/global_convinfo.txt
    OZINFO=${FIX_ROOT}/global_ozinfo.txt
    PCPINFO=${FIX_ROOT}/global_pcpinfo.txt

    #  copy Fixed fields to working directory
     cp $ANAVINFO ${BK_ROOT}/anavinfo
     cp $BERROR   ${BK_ROOT}/berror_stats
     cp $SATANGL  ${BK_ROOT}/satbias_angle
     cp $SATINFO  ${BK_ROOT}/satinfo
     cp $CONVINFO ${BK_ROOT}/convinfo
     cp $OZINFO   ${BK_ROOT}/ozinfo
     cp $PCPINFO  ${BK_ROOT}/pcpinfo
     cp $OBERROR  ${BK_ROOT}/errtable
    #
    # modify the anavinfo vertical levels based on wrf_inout for WRF ARW and NMM
    if [ ${bk_core} = ARW ] || [ ${bk_core} = NMM ] ; then
        anavlevels=`cat ${BK_ROOT}/anavinfo | grep ' sf ' | tail -1 | awk '{print $2}' `  # levels of sf, vp, u, v, t, etc
        anavlevels_stag=`cat ${BK_ROOT}/anavinfo | grep ' prse ' | tail -1 | awk '{print $2}' `  # levels of prse
        sed -i 's/ '$anavlevels'/ '$bklevels'/g' ${BK_ROOT}/anavinfo
        sed -i 's/ '$anavlevels_stag'/ '$bklevels_stag'/g' ${BK_ROOT}/anavinfo
    fi
    #
    # CRTM Spectral and Transmittance coefficients
    CRTM_ROOT_ORDER=${CRTM_ROOT}/${BYTE_ORDER}
    emiscoef_IRwater=${CRTM_ROOT_ORDER}/Nalli.IRwater.EmisCoeff.bin
    emiscoef_IRice=${CRTM_ROOT_ORDER}/NPOESS.IRice.EmisCoeff.bin
    emiscoef_IRland=${CRTM_ROOT_ORDER}/NPOESS.IRland.EmisCoeff.bin
    emiscoef_IRsnow=${CRTM_ROOT_ORDER}/NPOESS.IRsnow.EmisCoeff.bin
    emiscoef_VISice=${CRTM_ROOT_ORDER}/NPOESS.VISice.EmisCoeff.bin
    emiscoef_VISland=${CRTM_ROOT_ORDER}/NPOESS.VISland.EmisCoeff.bin
    emiscoef_VISsnow=${CRTM_ROOT_ORDER}/NPOESS.VISsnow.EmisCoeff.bin
    emiscoef_VISwater=${CRTM_ROOT_ORDER}/NPOESS.VISwater.EmisCoeff.bin
    emiscoef_MWwater=${CRTM_ROOT_ORDER}/FASTEM6.MWwater.EmisCoeff.bin
    aercoef=${CRTM_ROOT_ORDER}/AerosolCoeff.bin
    cldcoef=${CRTM_ROOT_ORDER}/CloudCoeff.bin

    if [ ! -e ${BK_ROOT}/Nalli.IRwater.EmisCoeff.bin ] ; then  
        ln -s $emiscoef_IRwater ${BK_ROOT}/Nalli.IRwater.EmisCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/NPOESS.IRice.EmisCoeff.bin ] ; then  
        ln -s $emiscoef_IRice ${BK_ROOT}/NPOESS.IRice.EmisCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/NPOESS.IRsnow.EmisCoeff.bin ] ; then  
        ln -s $emiscoef_IRsnow ${BK_ROOT}/NPOESS.IRsnow.EmisCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/NPOESS.IRland.EmisCoeff.bin ] ; then  
        ln -s $emiscoef_IRland ${BK_ROOT}/NPOESS.IRland.EmisCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/NPOESS.VISice.EmisCoeff.bin ] ; then  
        ln -s $emiscoef_VISice ${BK_ROOT}/NPOESS.VISice.EmisCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/NPOESS.VISland.EmisCoeff.bin ] ; then  
        ln -s $emiscoef_VISland ${BK_ROOT}/NPOESS.VISland.EmisCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/NPOESS.VISsnow.EmisCoeff.bin ] ; then  
        ln -s $emiscoef_VISsnow ${BK_ROOT}/NPOESS.VISsnow.EmisCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/NPOESS.VISwater.EmisCoeff.bin ] ; then  
        ln -s $emiscoef_VISwater ${BK_ROOT}/NPOESS.VISwater.EmisCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/FASTEM6.MWwater.EmisCoeff.bin ] ; then  
        ln -s $emiscoef_MWwater ${BK_ROOT}/FASTEM6.MWwater.EmisCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/AerosolCoeff.bin ] ; then  
        ln -s $aercoef  ${BK_ROOT}/AerosolCoeff.bin
    fi
    if [ ! -e ${BK_ROOT}/CloudCoeff.bin ] ; then  
        ln -s $cldcoef  ${BK_ROOT}/CloudCoeff.bin
    fi

    # Copy CRTM coefficient files based on entries in satinfo file
    for file in `awk '{if($1!~"!"){print $1}}' ${BK_ROOT}/satinfo | sort | uniq` ;do
        if [ ! -e ${BK_ROOT}/${file}.SpcCoeff.bin ] ; then  
            ln -s ${CRTM_ROOT_ORDER}/${file}.SpcCoeff.bin ${BK_ROOT}/
        fi
        if [ ! -e ${BK_ROOT}/${file}.TauCoeff.bin ] ; then  
            ln -s ${CRTM_ROOT_ORDER}/${file}.TauCoeff.bin ${BK_ROOT}/
        fi
    done

    # Only need this file for single obs test
     bufrtable=${FIX_ROOT}/prepobs_prep.bufrtable
     cp $bufrtable ${BK_ROOT}/prepobs_prep.bufrtable
    # for satellite bias correction
    # Users may need to modify these to use their own bias corrections
    cp ${FIX_ROOT}/gdas1.t12z.abias ${BK_ROOT}/satbias_in
    cp ${FIX_ROOT}/gdas1.t12z.abias_pc ${BK_ROOT}/satbias_pc
#


