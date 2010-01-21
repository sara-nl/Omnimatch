#!/bin/bash 
# 1.1 - MB - 05/04/08 OMNIMATCH suite based on emutd
#
###################################################
# omnimatch.sh
#
# This script gets executed on the "head" node and is called from the wrapper
# program.  It perfroms the following actions: 0) read in the wanted
# parameters 1) Initialize the connection to SRB 2) Get the software suite and
# input data 3) perfrom the analysis 4) store the output data on SRB
#
set -x

VERSION=1.1
UPDATED="08/04/05 by M. Bouwhuis"

#-----------------------------------------------------------------------
# some useful routines
#
function report_or_error  {
  if [ -z "$2" ] ; then 
    printf "%25s = NOT SET\n" "$1"
    exit 1 
  else 
    printf "%25s = \"%s\"\n" "$1" "$2"
  fi 
}

#-----------------------------------------------------------------------
# read in parameters from omnimatch.conf
#

CONF=omnimatch.conf
if [ -r "$CONF" ] ; then 
  source $CONF 
else
  echo ERROR Can not read $CONF. 
  exit 1 
fi 

if [ ! -z "$GLOBUS_GRAM_MYJOB_CONTACT" ] ; then 
  report_or_error "Number of Processors (grid env)"  "$NP"
  if [ ! -z "$NOPROC" ] ; then 
    echo Variable NOPROC is also set in conffile to $NOPROC ,
    echo but using environment NP = $NP
  fi 
  export NOPROC=$NP
else 
  report_or_error "Number of Processors (conf file)"  "$NOPROC"
fi

report_or_error "Input Data Collection" "$srb_data_dir"
report_or_error "Input Data Files"      "$srb_data_file"
report_or_error "The program command "  "$omnimatch_cmd"
report_or_error "The program arguments" "$omnimatch_args"


if [ ! -z "$PBS_NODEFILE" ] ; then 
  cat $PBS_NODEFILE
fi

#-----------------------------------------------------------------------
# start the SRB setup
#
echo +++++ SRB setup
MDASFILE=${HOME}/.srb/.MdasEnv
MDASDIR=`dirname $MDASFILE` 

if [ ! -z "$GLOBUS_GRAM_MYJOB_CONTACT" ] ; then 
  echo Welcome to a GRID job 
  if [ -f ./mdasEnv ] ; then 
    mkdir -p ${HOME}/.srb
    cp -v ./mdasEnv  $MDASFILE
  else  
    echo ERROR missing ./mdasEnv
    exit 1 
  fi
else 
  echo Welcome to a normal job 
fi



echo +++++ PWD = `pwd`

echo +++++  SRB initilization
which Sinit
Sinit 
if [ $? -ne 0 ] ; then 
  echo ERROR while doing Sinit 
  exit 1 
fi


#------------------------------------------------------------------------
# Getting the OMNIMATCH data out of the SRB
#
echo +++++++ Definition of the wrapper command 
mpi_wrapper="./wrapper"
all_execs="$mpi_wrapper $omnimatch_cmd"


echo +++++++ Getting the data out
if [ "$srb_data_file" = "ALL" ] ; then 
  the_data_file=`Sls ${srb_data_dir}  | grep -v "^  C-/" | grep -v ':$' | tr -d '\n'`
else
  the_data_file="$srb_data_file"
fi


for afile in $the_data_file ; do 
  Sls ${srb_data_dir}/${afile}
  if [ $? -ne 0 ] ; then
    echo ERROR finding srb:${srb_data_dir}/${afile}
    exit 1 
  fi

  Sget ${srb_data_dir}/${afile}
  if [ $? -ne 0 ] ; then
    echo ERROR getting srb:${srb_data_dir}/${afile}
    exit 1
  fi
done


#-----------------------------------------------------------------------
#
#
echo ++++++ Setup Omnimatch software

for afile in $all_execs ; do 
  if [ -r "$afile" ] ; then 
    chmod u+x $afile
  else
    echo ERROR can not find the $afile
    exit 1 
  fi
done


echo ++++ .files
ls -ltra


echo ++++++ Now start the command  `date`
mpirun -np $NOPROC -machinefile $PBS_NODEFILE $omnimatch_cmd $omnimatch_args > omnimatch.outerr 2>&1
echo ++++++ done `date`

ls -ltra 
echo ++++++ Now pack it all up  


dateext=`date +"%y%m%d-%H%M%S"`

mv -v omnimatch.outerr omnimatch.outerr.$dateext

touch jobuid.${dateext}.`basename $EDG_WL_JOBID`

Sput -v omnimatch.outerr.$dateext  ${srb_data_dir}
Sput -v jobuid.${dateext}.`basename $EDG_WL_JOBID`
Sput -v -f ${omnimatch_out}.ang ${srb_data_dir}
Sput -v -f ${omnimatch_out}.ccf ${srb_data_dir}
Schmod -r a emutd groups ${srb_data_dir}
Sls -l ${srb_data_dir}

echo ++++++ done succes

