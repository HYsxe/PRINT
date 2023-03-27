#!/bin/bash
for i in {1..10}
do
   sbatch runTFBS.sh $i
done