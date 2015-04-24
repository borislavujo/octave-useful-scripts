#PBS -N drps-aladip-test
#PBS -q test
##PBS -l walltime=144:00:00

cd ${PBS_O_WORKDIR}

echo Starting job 
echo PBS assigned me this node:
cat $PBS_NODEFILE
echo "Running DRPS on Ala dipeptide"
echo "short simulation in implicit solvent"
octave shortSimul.m > output
echo "Job finished. PBS details are:"
echo
qstat -f ${PBS_JOBID}
echo
echo Finished at `date`
