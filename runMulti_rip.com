#!/bin/bash

#for running years on Rip
nProc=12 #This is the number of processes

#what is name of python script
#runScript=ripRun_multiyear.py
runScript=particleTrack_coastAware.py

#years over which to run
#yearList='2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023'
yearList='2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023'
#yearList='2007'
machine='surge'
depth='10'


for year in $yearList
do
	      echo "####### Year: $year depth: $depth Machine: $machine nProc: $nProc  #######"
	      #nice -n 10 mpirun -np 14 python particleTrack_coastAware.py $year
	      mpirun -np $nProc python particleTrack_coastAware.py $year $machine $depth 2>&1 > runLog.txt
	      #wait
	      echo "done with $year"
done
	      
