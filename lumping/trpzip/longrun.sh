#PBS -N lumping-jaj-exp-g6
#PBS -q l1
#PBS -l walltime=144:00:00
#PBS -l mem=8gb

cd ${PBS_O_WORKDIR}
#kery=$(date +%Y%m%d-%H%M%S)
#TMP=/scratch/bf269/$PBS_JOBID
#mkdir -p $TMP
#cp -r * $TMP
#cd $TMP
echo PBS assigned me this node:
cat $PBS_NODEFILE
octave -q tauRxnExp.m > output
echo "Job finished. PBS details are:"
qstat -f ${PBS_JOBID}
echo Finished at `date`
#cp -r * ${PBS_O_WORKDIR}/
#cd $PBS_O_WORKDIR
#rm -rf $TMP

