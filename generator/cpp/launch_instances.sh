#!/bin/bash

binary="MUBQPEval-heuristic-par"

########## Checking parameters ##########
if [ ! $# -eq 5 ]
then
	echo "usage: launch_instances.sh [number of instances] [instance duration (seconds)] [number of processes] [instance file] [pool size]"
	exit -1
fi
number_of_instances=$1
duration=$2
procs=$3
instance=$4
pool_size=$5

if [ $procs -lt 2 ]
then
	echo "minimum number of processes: 2 (it's a master/slaves schema)"
	exit -1
fi
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
	mpirun -np $procs ./$binary $instance $pool_size > $output &
	pid=$!
	sleep $duration
	kill -s SIGUSR1 $pid
done
