#!/bin/sh
#$ -q gen4.q,gen5.q,gen6.q
#$ -N Resub
#$ -S /bin/sh
#$ -cwd

. /etc/profile.d/modules.sh

# Disable production of core dump files
ulimit -c 0

echo ""
echo "***********************************************************************"
echo " Starting job:"
echo ""
echo "    Name:              "$JOB_NAME
echo "    ID:                "$JOB_ID
echo "    Hostname:          "$HOSTNAME
echo "    Working directory: "$SGE_O_WORKDIR
echo ""
echo "    Submitted using:   submit -N Resub gen4.q,gen5.q,gen6.q molpro@2010.1.67+mpi"
echo "***********************************************************************"

vulcan load molpro@2010.1.67+mpi

export NSLOTS=4

molpro -n $NSLOTS --nouse-logfile --no-xml-output -o output.dat input.dat
