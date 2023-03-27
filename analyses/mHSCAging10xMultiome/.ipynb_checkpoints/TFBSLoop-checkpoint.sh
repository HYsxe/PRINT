#!/bin/bash
for i in {1..16}
do
   sbatch runTFBS.sh $i
done
