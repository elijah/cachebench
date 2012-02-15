#!/bin/sh

GRAPHS="latency roundtrip bandwidth bibw alltoall broadcast reduce allreduce"

usage()
{
  echo "Usage: $0 <datafile_prefix> [latency,roundtrip,bandwidth,bibw,alltoall,broadcast,reduce,allreduce]"
  exit 1
}

if [ $# -ne 1 -a $# -ne 2 ]; then
  usage;
fi

PREFIX=$1
if [ $# -eq 2 ]; then
  GRAPHS=$2
fi

for g in $GRAPHS
do
  echo "gnuplot < ${PREFIX}_${g}.gp"
  gnuplot < ${PREFIX}_${g}.gp
  if [ $? -ne 0 ]; then
    gperr=1
  else
    echo "Postscript graph is in ${PREFIX}_${g}.ps."
    gperr=0
  fi
done

if [ $gperr -eq 1 ]; then
    echo ""
    echo "GNUplot does not work/exist on this machine. In this directory you"
    echo "will find the file ${PREFIX}-mpbench-datafiles.tar. Please move "
    echo "this file to another machine, un-tar it, and for every <file.gp>:"
    echo ""
    echo "gnuplot < results/<file>.gp"
    echo ""
fi

exit 0
