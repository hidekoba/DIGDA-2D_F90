#!/bin/bash

DIGDADIR=/work/gi55/c24116/digda/DIGDA-2D_F90
BINARY=$DIGDADIR/bin/digda_128x64.par
RUNDIR=$DIGDADIR/RUN
BNDDIR=$DIGDADIR/Bound/data_digda128064
PARAMET=./PARAMET_DIGDA
MPIRUN=mpiexec

SDYRSTR="0"
SDYREND="10000"
SDYRHST="2000"

#SDTSTEP="0.002739726027397260"
#SDTSTEP="0.001369863013698630"

EXPDIR=$DIGDADIR/OUT/3d

INIDIR=$EXPDIR/$SDYRSTR
OUTDIR=$EXPDIR/$SDYREND

mkdir -p $OUTDIR

#cd $RUNDIR
#mv digda.[0-9]*.out digda.[0-9]*.stats exp

cd $OUTDIR
/bin/rm -f $PARAMET ; /bin/rm -rf *

JOB_NAME=`basename $0`
/bin/rm -rf ${JOB_NAME}.y${SDYREND}*
#############################################################################
cat << EOF >> $PARAMET
&nmrun crun='digda test' &end
&nmocbc bound_dir='$BNDDIR' &end
&nmocrt init_dir='$INIDIR',
        rstt_dir='$OUTDIR' &end
&nmsdtm sdtmstp=1, iutstp=3,
        sdyrstr=$SDYR,
        sdyrend=$SDYREND,
        wotstp=t &end
&nmsdhi sdyrhst=$SDYRHST,
        wohist=t,
        hstavg=t &end
&nmsddear deiarc=t &end
&nmsdiso calc13=f, calc14=f &end

&nmsdpm kox=2.00e-9, kax=7.90e-11 &end
&nmsdca rccar=365.d-3, ncar=4.50  &end
&nmsdop rcopa=5.0d-3, roprid=f, roparc=t &end

&nmsdchem alkerr=0.10d0, itrmax=10000, rocn=1.0245d0 &end

&nmsdso sadu=f, sadq=t, sadqf=t, sdcn=f, sdii=t &end
&nmsdli cdcn=f, cdii=t &end

EOF
#############################################################################
cat << EOF > ${JOB_NAME}.y${SDYREND}
#!/bin/tcsh

#PJM -g gi55o
#PJM --fs /work
#PJM -j
#PJM -N digda
#PJM -s

#PJM -L rscgrp=regular-o
#PJM -L elapse=48:00:00
#PJM -L node=1
#PJM -L node-mem=28GiB
#PJM --mpi proc=16
#PJM --mpi "max-proc-per-node=16"
#PJM --omp thread=1

EOF
#############################################################################
cat << EOF >> ${JOB_NAME}.y${SDYREND}
umask 002

module unload intel
module unload impi
module load fj
module load fjmpi
module load hdf5
module load phdf5
module load netcdf
module load netcdf-fortran

setenv PLE_MPI_STD_EMPTYFILE off 
setenv PARALLEL        1 
setenv OMP_NUM_THREADS 1 
setenv FORT90L "-Wl,-T" # big-endian

EOF
#############################################################################
cat << EOF >> ${JOB_NAME}.y${SDYREND}
setenv BINARY    $BINARY
setenv RUNDIR    $RUNDIR
setenv BNDDIR    $BNDDIR
setenv PARAMET   $PARAMET
setenv MPIRUN    $MPIRUN

set SDYRSTR=$SDYRSTR
set SDYREND=$SDYREND

setenv EXPDIR    $EXPDIR

setenv INIDIR    $INIDIR
setenv OUTDIR    $OUTDIR

time $MPIRUN $BINARY >& ERROUT
mkdir ./STDOUT
mv STDOUT.* ./STDOUT

EOF
#######################################################################

VAR=`pjsub --step ${JOB_NAME}.y${SDYREND} -e ERROUT`
echo $VAR
VAR=`echo $VAR | cut -d " " -f 6`
JID=${VAR%_0}
echo $JID
echo "pjsub Job $JID submitted at `date`."

while true; do
  JOB_STATUS=$(pjstat $JID | grep -oP '\b[A-Z]+\b' | tail -1)
 #echo $JOB_STATUS
  if [ "$JOB_STATUS" == "BDEC" ] ; then
    echo "Job $JID has finished with status: $JOB_STATUS"
    break
  fi
  sleep 10
done

echo "post-process `basename $OUTDIR` started at `date`."

bash $RUNDIR/pp_digda_combine_SEDTRR.sh $OUTDIR
bash $RUNDIR/pp_digda_calc_weight_percentage.sh $OUTDIR
