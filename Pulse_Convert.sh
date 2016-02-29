#PBS -N Pulse_Convert
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -l walltime=10:00:00
#PBS -q atlas-6
#PBS -k oe
#PBS -m abe

MYDIR=$HOME/data2/testdata

EXE=$MYDIR/converter.exe 
PATH="$MYMPI/bin:$PATH"; export PATH

cd $MYDIR
ulimit -s unlimited

echo "Started on '/bin/hostname'"
echo
echo "PATH is [$PATH]"
echo
echo "Nodes chosen are:"
cat $PBS_NODEFILE
$EXE

