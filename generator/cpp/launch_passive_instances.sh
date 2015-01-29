#!/bin/bash
# To be launched with: oarsub -l /nodes=[number of nodes],walltime=xx:xx:xx launch_passive_instances.sh

binary="MUBQPEval-heuristic-par"
machinefile="machines.txt"
if [ ! -f $OAR_FILE_NODES ]
then
	echo "No reservation made. Launch with: oarsub -l /nodes=[number of nodes],walltime=xx:xx:xx launch_passive_instances.sh"
	exit -1
fi
cat $OAR_FILE_NODES > $machinefile # ok for MUBQPEval-heuristic-par.cpp (will always be more than 2 cores, we reserved 8-cores machines)
procs=$(wc -l $machinefile | cut -d " " -f 1)

########## Checking parameters ##########
if [ ! $# -eq 5 ] 
then
	echo "usage: launch_instances.sh [number of instances] [instance duration] [instance file] [pool size]"
	exit -1
fi
number_of_instances=$1
duration=$2
instance=$3
pool_size=$4

if [ ! -f $instance ]
then
	echo "$instance: file not found"
	exit -1
fi

########## Looking for the binary to execute ##########
if [ ! -f $binary ]
then
	echo "$binary: file not found. Compile the program first."
	exit -1
fi

########## Main ##########
for (( instance_number=1; instance_number<=$number_of_instances; instance_number++ ))
do
	########## Building output file name ##########
	output=$(basename $instance)
	output="${output%.*}_"$duration"_"$procs"_"$pool_size"_"$instance_number".txt"

	########## Launching instance ##########
	mpirun --mca plm_rsh_agent "oarsh" --mca pml ob1 --mca btl tcp,self -machinefile $machinefile ./$binary $instance $pool_size >> $output &

	pid=$!
	sleep $duration
	kill -s SIGUSR1 $pid
done
