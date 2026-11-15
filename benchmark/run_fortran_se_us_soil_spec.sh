#!/bin/bash
# Official MEGSEA (Yienger-Levy soil NO, BDSNP_YN=N; with and without soil data) and mgn2mech (RUN_CONVERSION=N,
# SPCTONHR=Y: 201 MEGAN species in tonnes/hr) on the SE_US test case. Needs MGNERS.ncf from run_fortran_se_us.sh in <outdir>.
# usage: run_fortran_se_us_soil_spec.sh <outdir>
set -e
M=${MGN:-$HOME/scratch/megan321/MEGANv3.21}; OUT=$(cd $1 && pwd); HERE=$(cd $(dirname $0) && pwd)
export LD_LIBRARY_PATH=${OMP_LIB:-$HOME/.local/lib}:$LD_LIBRARY_PATH PROMPTFLAG=N EXECUTION_ID=bench
INP=$M/Input; JD=2018122; MET=$INP/MGNMET/MET.MEGAN.SE_US.rad45.$JD.ncf
export GRIDDESC=$INP/MAP/GRIDDESC GDNAM3D=SE_US SDATE=$JD STIME=0 RLENG=240000
python3 -c "import sys; sys.path.insert(0, '$HERE'); from ioapi import drop_vars; drop_vars('$MET', '$OUT/MET_nosoil.ncf', ('SOIM1', 'SOIT1', 'SLTYP'))"
for path in soil nosoil; do
  MM=$MET; [ $path = nosoil ] && MM=$OUT/MET_nosoil.ncf
  ( export NLAI=46 BDSNP_YN=N EPIC=N CANTYPFILE=$M/work/CANTYPFILE LAIFILE=$INP/MAP/LAI3_SE_US.ncf MGNMET=$MM MGNSEA=$OUT/MGNSEA_YL_$path.ncf
    rm -f $MGNSEA; $M/src/MEGSEA/megsea > $OUT/log_megsea_$path.txt 2>&1 )
done
( export TSTEP=10000 RUN_SPECIATE=Y RUN_CONVERSION=N SPCTONHR=Y BDSNP_YN=N MECHANISM=CB6X EFMAPS=$INP/MAP/EFMAP.2019b.SE_US.ncf \
         MGNERS=$OUT/MGNERS.ncf MGNSEA=$OUT/MGNSEA_YL_nosoil.ncf MGNOUT=$OUT/MGN201_tonhr.ncf
  rm -f $MGNOUT; $M/src/MGN2MECH/mgn2mech > $OUT/log_mgn2mech.txt 2>&1 )
for p in megsea_soil megsea_nosoil mgn2mech; do printf "%s: " $p; grep -c "Normal Completion" $OUT/log_$p.txt || true; done
