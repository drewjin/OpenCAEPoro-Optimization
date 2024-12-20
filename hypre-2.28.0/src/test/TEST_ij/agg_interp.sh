#!/bin/sh
# Copyright (c) 1998 Lawrence Livermore National Security, LLC and other
# HYPRE Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

TNAME=`basename $0 .sh`
RTOL=$1
ATOL=$2

#=============================================================================
# compare with baseline case
#=============================================================================

FILES="\
 ${TNAME}.out.1\
 ${TNAME}.out.2\
 ${TNAME}.out.3\
 ${TNAME}.out.4\
 ${TNAME}.out.4.1\
 ${TNAME}.out.4.2\
 ${TNAME}.out.5\
 ${TNAME}.out.6\
 ${TNAME}.out.7\
 ${TNAME}.out.8\
 ${TNAME}.out.8.1\
 ${TNAME}.out.8.2\
 ${TNAME}.out.9\
 ${TNAME}.out.10\
 ${TNAME}.out.10.1\
 ${TNAME}.out.10.2\
 ${TNAME}.out.11\
 ${TNAME}.out.12\
 ${TNAME}.out.12.1\
 ${TNAME}.out.12.2\
 ${TNAME}.out.13\
 ${TNAME}.out.14\
 ${TNAME}.out.15\
 ${TNAME}.out.16\
 ${TNAME}.out.17\
 ${TNAME}.out.18\
 ${TNAME}.out.19\
 ${TNAME}.out.20\
 ${TNAME}.out.21\
 ${TNAME}.out.22\
"

for i in $FILES
do
  echo "# Output file: $i"
  tail -5 $i 
done > ${TNAME}.out

# Make sure that the output file is reasonable
RUNCOUNT=`echo $FILES | wc -w`
OUTCOUNT=`grep "Iterations" ${TNAME}.out | wc -l`
if [ "$OUTCOUNT" != "$RUNCOUNT" ]; then
   echo "Incorrect number of runs in ${TNAME}.out" >&2
fi

#=============================================================================
# remove temporary files
#=============================================================================

rm -f ${TNAME}.testdata*
