#!/bin/bash
# Runs the official MEGANv3.21 binaries (as shipped in Zenodo 10.5281/zenodo.10526206) on the SE_US test case.
# Usage: run_fortran.sh <outdir> [GAMSM_YN] [GAMBD_YN] [NLAI]
set -e
M=${MGN:-$HOME/scratch/megan321/MEGANv3.21}; OUT=$1; SM=${2:-N}; BD=${3:-N}; NL=${4:-46}
export LD_LIBRARY_PATH=${OMP_LIB:-$HOME/.local/lib}:$LD_LIBRARY_PATH PROMPTFLAG=N EXECUTION_ID=bench
mkdir -p $OUT; INP=$M/Input; JD=2018122
export GRIDDESC=$INP/MAP/GRIDDESC GDNAM3D=SE_US
# DAYMET
( export EPISDATE=$JD STIME=0 RLENG=230000 NDAYS=1 MGNMET001=$INP/MGNMET/MET.MEGAN.SE_US.rad45.$JD.ncf DailyMET=$OUT/DAYMET.ncf
  rm -f $DailyMET; $M/src/DAYMET/daymet > $OUT/log_daymet.txt 2>&1 )
# MEGCAN
( export NLAI=$NL SDATE=$JD STIME=0 RLENG=240000 CANTYPFILE=$M/work/CANTYPFILE LAIFILE=$INP/MAP/LAI3_SE_US.ncf \
         MGNMET=$INP/MGNMET/MET.MEGAN.SE_US.rad45.$JD.ncf CANMET=$OUT/CANMET.ncf
  rm -f $CANMET; $M/src/MEGCAN/megcan > $OUT/log_megcan.txt 2>&1 )
# MEGVEA
( export Layers=5 NLAI=$NL GAMAQ_YN=N GAMCO2_YN=N GAMHW_YN=N GAMHT_YN=N GAMLT_YN=N GAMSM_YN=$SM GAMBD_YN=$BD \
         SDATE=$JD STIME=0 RLENG=240000 N_MaxT=1 N_MinT=1 N_MaxWS=1 LAIFILE=$INP/MAP/LAI3_SE_US.ncf CANMET=$OUT/CANMET.ncf \
         DailyMET=$OUT/DAYMET.ncf SMFILE=${SMFILE:-$M/Output/INT/MGNSEA.SE_US.$JD.ncf} AQFILE=$INP/MAP/W126_SE_US.ncf \
         LDFILE=$INP/MAP/LDF_SE_US.2019b.ncf MGNERS=$OUT/MGNERS.ncf
  rm -f $MGNERS; $M/src/MEGVEA/megvea > $OUT/log_megvea.txt 2>&1 )
echo "fortran done -> $OUT"; tail -3 $OUT/log_megvea.txt
