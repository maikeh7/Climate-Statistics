#!/bin/bash

for i in {1976..2005}
  do
  echo "running year $i"
  nohup Rscript Main_driver_StrawMan.R $i > nohupSM-$i.out &
  echo "done year $i"
  done


