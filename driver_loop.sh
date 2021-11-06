#!/bin/bash

for i in 2002 2003 2004 2005 2006 2007 2008 2009 2010
  do
  echo "running year $i"
  nohup Rscript Main_driver.R $i > nohup-$i.out &
  done

