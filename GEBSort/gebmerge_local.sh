#!/bin/sh

if [ $# -ne 1 ] 
  then
   echo "Specify only the run number to sort"
  exit 1
fi

RUN=$1

# regenerate all

#rm GTMerge3 mkMap zzip zunzip
#make clean; 
#make GTMerge3 ()&&
#make zzip 
#make zunzip
#make mkMap

#make clean; 
#make DGSSort

# make the map file

#./mkMap > map.dat


# data from DGS, new data format (firmware) 9/25/2012

echo "GEBMerge started at `date`"
  #rm DATA; ln -s runs  DATA
  #ls DATA/merged_run1.gtd* 
  #rm DATA/merged_run1.gtd*
#  ./GEBMerge GEBMerge.chat ../Merged/GEBMerged_run$RUN.gtd `ls ../dfmadata/dfmarun_$RUN.*`  `ls ../dgsdata/dgsrun_$RUN.gt*` `ls ../xadata/xarun_$RUN.gt*` 2> LOG_FILES/GEBMerge_run$RUN.log
 ./GEBMerge GEBMerge.chat ../Merged/GEBMerged_run$RUN.gtd `ls ../data/music40_22Mg_run_$RUN.*` 2> LOG_FILES/GEBMerge_run$RUN.log
# ./GEBMerge GEBMerge.chat ../Merged/GEBMerged_xa_run$RUN.gtd `ls ../xadata/xarun_$RUN.gt*` 2> LOG_FILES/GEBMerge_run$RUN.log
  #ls -lh merged_run1.gtd*
echo "GEBMerge DONE at `date`"

#Create input and output files
#touch ifile1.dat ofile1.dat
#cat >ifile1.dat <<ENDI
#input disk /media/20140610/user/MCP/merged_run$RUN.gtd_000
#ENDI

#cat >ofile1.dat <<ENDO
#output ROOT_FILES/run$RUN.root
#ENDO

#echo "DGSSort started at `date`"
#  #rm ROOT_FILES/run1.root
#  ./DGSSort -chat DGSSort_run.chat > LOG_FILES/DGSSort_run$RUN.log
#  #ls -lh ROOT_FILES/run5.root
#echo "DGSSort DONE at `date`"

exit


