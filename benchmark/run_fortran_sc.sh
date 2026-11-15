#!/bin/bash
# Official MEGANv3.21 binaries on the Santa Catarina benchmark inputs (lat-lon IOAPI). usage: run_fortran_sc.sh <dir>
set -e
M=${MGN:-$HOME/scratch/megan321/MEGANv3.21}; D=$(cd $1 && pwd); SD=$(python3 -c "import numpy as np;print(int(np.load('$D/py_inputs.npz')['sdate']))")
export LD_LIBRARY_PATH=${OMP_LIB:-$HOME/.local/lib}:$LD_LIBRARY_PATH PROMPTFLAG=N EXECUTION_ID=bench GRIDDESC=$D/GRIDDESC GDNAM3D=SC_BENCH
mkdir -p $D/fort; O=$D/fort
( export EPISDATE=$SD STIME=0 RLENG=240000 NDAYS=1 MGNMET001=$D/MET.ncf DailyMET=$O/DAYMET.ncf; rm -f $DailyMET; $M/src/DAYMET/daymet > $O/log_daymet.txt 2>&1 )
( export NLAI=12 SDATE=$SD STIME=0 RLENG=240000 CANTYPFILE=$D/CANTYP.ncf LAIFILE=$D/LAI.ncf MGNMET=$D/MET.ncf CANMET=$O/CANMET.ncf; rm -f $CANMET; $M/src/MEGCAN/megcan > $O/log_megcan.txt 2>&1 )
( export NLAI=12 SDATE=$SD STIME=0 RLENG=240000 BDSNP_YN=N EPIC=N CANTYPFILE=$D/CANTYP.ncf LAIFILE=$D/LAI.ncf MGNMET=$D/MET.ncf MGNSEA=$O/MGNSEA.ncf; rm -f $MGNSEA; $M/src/MEGSEA/megsea > $O/log_megsea.txt 2>&1 )
( export Layers=5 NLAI=12 GAMAQ_YN=N GAMCO2_YN=N GAMHW_YN=N GAMHT_YN=N GAMLT_YN=N GAMSM_YN=N GAMBD_YN=Y SDATE=$SD STIME=0 RLENG=240000 \
         N_MaxT=1 N_MinT=1 N_MaxWS=1 LAIFILE=$D/LAI.ncf CANMET=$O/CANMET.ncf DailyMET=$O/DAYMET.ncf SMFILE=$O/MGNSEA.ncf \
         LDFILE=$D/LDF.ncf MGNERS=$O/MGNERS.ncf; rm -f $MGNERS; $M/src/MEGVEA/megvea > $O/log_megvea.txt 2>&1 )
for p in daymet megcan megsea megvea; do printf "%s: " $p; grep -c "Normal Completion" $O/log_$p.txt; done
