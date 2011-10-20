#!/bin/bash

cmd="java -cp bin:lib/collections-generic-4.01.jar:lib/java-getopt-1.0.13.jar:lib/jung-algorithms-2.0.1.jar:lib/jung-api-2.0.1.jar:lib/jung-graph-impl-2.0.1.jar:lib/Jaligner.jar TransAssembly_allProbPaths -N 10000 -L 300 -F 300 -C sample_data/RawComps.0/comp0 --stderr -V 20 "

eval $cmd

if [ "$?" -ne "0" ]; then
    echo "Error, command failed: " $cmd
    exit 1
fi

echo
echo
echo "############################################################################"
echo "Butterfly output written to: sample_data/RawComps.0/comp0_allProbPaths.fasta"
echo "############################################################################"
echo
echo
