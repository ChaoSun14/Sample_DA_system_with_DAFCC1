#!/bin/csh -f

cd $NAMELIST_DST_DIR

set basedate_num = `echo $RUN_START_DATE | sed -e 's/-//g'`  

if !(-d $DATA_DST_DIR/input_data) mkdir -p $DATA_DST_DIR/input_data
if !(-d $DATA_DST_DIR/input_data/start) mkdir -p $DATA_DST_DIR/input_data/start
if !(-d $DATA_DST_DIR/output_data) mkdir -p $DATA_DST_DIR/output_data
if !(-d $DATA_DST_DIR/output_data/restart) mkdir -p $DATA_DST_DIR/output_data/restart

set start_year = `echo $RUN_START_DATE | awk -F '-' '{print $1}'`
set start_month = `echo $RUN_START_DATE | awk -F '-' '{print $2}'`
set start_day = `echo $RUN_START_DATE | awk -F '-' '{print $3}'`
@ start_hour_num = $RUN_START_SECOND / 3600
@ start_minute_num = ( $RUN_START_SECOND / 60 ) % 60
@ start_second_num = $RUN_START_SECOND % 60
set start_hour = "$start_hour_num"
set start_minute = "$start_minute_num"
set start_second = "$start_second_num"
if ( $start_hour_num < 10 ) set start_hour = "0$start_hour_num"
if ( $start_minute_num < 10 ) set start_minute = "0$start_minute_num"
if ( $start_second_num < 10 ) set start_second = "0$start_second_num"


set stop_year = `echo $RUN_STOP_DATE | awk -F '-' '{print $1}'`
set stop_month = `echo $RUN_STOP_DATE | awk -F '-' '{print $2}'`
set stop_day = `echo $RUN_STOP_DATE | awk -F '-' '{print $3}'`
@ stop_hour_num = $RUN_STOP_SECOND / 3600
@ stop_minute_num = ( $RUN_STOP_SECOND / 60 ) % 60
@ stop_second_num = $RUN_STOP_SECOND % 60
set stop_hour = "$stop_hour_num"
set stop_minute = "$stop_minute_num"
set stop_second = "$stop_second_num"
if ( $stop_hour_num < 10 ) set stop_hour = "0$stop_hour_num"
if ( $stop_minute_num < 10 ) set stop_minute = "0$stop_minute_num"
if ( $stop_second_num < 10 ) set stop_second = "0$stop_second_num"


link_data "$DATA_SRC_DIR/demo/climatic/bc" "$DATA_DST_DIR/input_data/"
link_data "$DATA_SRC_DIR/demo/climatic/csurf-bndry" "$DATA_DST_DIR/input_data/"
link_data "$DATA_SRC_DIR/demo/climatic/initial" "$DATA_DST_DIR/input_data/"
#link_data "$DATA_SRC_DIR/demo/realtime/$start_year$start_month$start_day$start_hour/tide" "$DATA_DST_DIR/input_data/"
#link_data "$DATA_SRC_DIR/demo/realtime/$start_year$start_month$start_day$start_hour/wind" "$DATA_DST_DIR/input_data/"
link_data "$DATA_SRC_DIR/demo/realtime/$start_year$start_month$start_day$start_hour/start/DEMO.res.$start_year$start_month$start_day-$start_hour$start_minute$start_second.nc" "$DATA_DST_DIR/input_data/start"


cat >! runsetocean.dat << EOF
!==============================================================================!
!   INPUT FILE FOR PARAMETERS CONTROLLING EXECUTION OF wave                    !
!   DESCRIPTION OF VARIABLES AND SUGGESTED PARAMETERS CAN BE FOUND AT BOTTOM   !
!                                                                              !
!        FORMAT:			                                       !
!       1.) VARIABLE  = VALUE  (EQUAL SIGN MUST BE USED)                       !
!       2.) FLOATING POINT VARIABLES MUST CONTAIN A PERIOD "." EX: 1.3, 2.,etc !
!       3.) BLANK LINES ARE IGNORED AS ARE LINES BEGINNING WITH ! (F90 COMMENT)!
!       4.) COMMENTS CAN FOLLOW VALUES IF MARKED BY !                          !
!       5.) ORDER OF VARIABLES IS NOT IMPORTANT                                !
!       6.) FOR MULTIPLE VALUE VARIABLES FIRST ENTRY IS NUMBER OF VARIABLES    !
!           TO FOLLOW (OR 0 IF NONE)                                           !
!       7.) DO NOT USE COMMAS TO SEPARATE VARIABLES                            !
!       8.) DO NOT EXCEED EIGHTY CHARACTERS PER LINE                           !
!       9.) FOR LINE CONTINUATION ADD \\\\ TO END OF LINE TO FORCE CONTINUE      !
!           TO NEXT LINE.  MAXIMUM 4 CONTINUATIONS                             !
!       10.) TRUE = T, FALSE = F                                               !
!                                                                              ! 
!  THE PREVIOUS FORMAT OF "VARIABLE: VALUE" IS NO LONGER VALID                 !
!  THE MORE ATTRACTIVE " = " FORMAT WAS SUGGESTED BY Hernan G. Arango          !
!    AND SUBSEQUENTLY ADOPTED                                                  !
!==============================================================================!


!============ Case Title========================================================

ModelVersion =  DEMO  !maximum number of character is 80  
RunType = $RUN_TYPE
Initial_DataFile = DEMO.res.$start_year$start_month$start_day-$start_hour$start_minute$start_second.nc
INPUT = ./input_data/
OUTPUT = ./output_data/
!=========Parameters of model region============================================
IGRID = 2
LON = 99. 150.
LAT = 0. 50. 
SIGMA_LAYER = 30
!=========INITIAL DATA OF TEMPERATRUE AND SALINITY==============================
INITIAL_TS = T     !TRUE(T) MEAN THAT INITIAL DATA OF TEMPERATRUE AND SALINITY
                   !ARE PREPARED, OR NOT 
!=========Parameters controlling parallel size of CPU===========================
NUMX = $num_x_proc
NUMY = $num_y_proc
!========CONTROL TIME===========================================================
ISPLIT = 30
TIME_STEP = 10.  !units is second #DTE's value
REAL_TIME = T    !IF T respresent that the forcing data is real time,or climate
START_TIME = $RUN_START_DATE ${start_hour}:${start_minute}:${start_second}
                 !start run time example: yyyy-mm-dd HH:MM:SS
END_TIME = $RUN_STOP_DATE ${stop_hour}:${stop_minute}:${stop_second}
                 !end run time example: yyyy-mm-dd HH:MM:SS
RUN_TIME = -70   !IF minus units is hour,IF plus units is day. 
                 !BUT all must be integer
!========CONTROL OUTFILE DATA====+=============================================
INTERVAL_OUT = 1          ! units hour
INTERVAL_RESTART = 120000     ! units hour
!========================control open boundary parameter ======================
!OBC_FLG = [S E N W] South East North West
OBC_FLG = T T T T
!==============================================================================
#OUT_DATA_NAME = ctime bv windx windy hs tp tz th
OUT_DATA_FIG =  T     T  T     T     T  T  T  T   
                                 ! OUT_DATA_FIG respresent that each OUT_DATA_NAME 
                                 ! is or not exproted. T is exporting, F is Not
!========CONTROL FORCE DATA=======+============================================
! units is hour
INTERVAL_FORCING = 1
!=======CONTROL TIDE ==========================================================
TIDETYPE = F
!=======CONTROL ISPADV ========================================================
ISPADV = 5
!=======CONTROL WIND STRESS ===================================================
WINDTYPE = F  !IF WINDTYPE IS T MEANS THAT USING CLIMATIC WIND STRESS
              !IF WINDTYPE IS F MEANS THAT USING REAL TIME WIND STRESS
!=======CONTROL HEAT FLUX   ===================================================
HEATTYPE = F  !IF HEATTYPE IS T MEANS THAT USING CLIMATIC HEAT FLUX
              !IF HEATTYPE IS F MEANS THAT USING REAL TIME HEAT FLUX 

!======= whether coupled with wave model using BV =============================
BVTYPE  = F   !IF BVTYPE IS T MEANS THAT USING BV from masnum wave model
              !IF BVTYPE IS F,will set BV=0.0 in the POM model
!=======CONTROL MODEL =========================================================
MODE = 3      !MODE = 2; 2-D CALCULATION (BOTTOM STRESS CALCULATED IN ADVAVE)
              !       3; 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
              !       4; 3-D CALCULATION WITH T AND S HELD FIXED
!========THE VERTICAL SIGMA GRID===============================================
! UNITS M
ZSS = 0. 10. 20. 30. 50. 75. 100. 125. 150. 200. 250.        \\\\ 
      300. 400. 500. 600. 700. 800. 900. 1000. 1100. 1200.   \\\\
      1300. 1400. 1500. 1750. 2000. 2500. 3000. 3500. 4000. 
EOF
