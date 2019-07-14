#!/bin/csh -efx
# epi.reader.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1995 Department of Statistics,         *
# *                        Carnegie Mellon University        *
# *                                                          *
# *  This program is distributed in the hope that it will    *
# *  be useful, but WITHOUT ANY WARRANTY; without even the   *
# *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
# *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
# *  nor any of the authors assume any liability for         *
# *  damages, incidental or otherwise, caused by the         *
# *  installation or use of this software.                   *
# *                                                          *
# *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
# *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
# *  FDA FOR ANY CLINICAL USE.                               *
# *                                                          *
# *                                                          *
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#
echo '$Id: ts.reader.csh,v 1.6 2003/07/18 21:27:59 welling Exp $'
#
echo '#'`date`$0

# Will need some contingency filenames
set rdr_output_fname = ${2}_${F_READER_DIMORDER}
set rdr_permute_fname = ${2}_${F_READER_SPACE}

#cat $1.Z | uncompress > trash
reader        -input $1 -dims $F_READER_DIMS \
	-dataout .dat -headerout ${rdr_output_fname}.mri \
	-offset $F_READER_OFFSET -reorder $F_READER_REORDER \
	-skip $F_READER_SKIP -sliceskip $F_READER_SSKIP -type $F_READER_TYPE \
	-veclen $F_READER_VECLEN -dataorder $F_READER_DIMORDER \
	$F_READER_ENDIAN
#rm trash

if ( "$F_READER_DIMORDER" != "vxyzt" ) then
    echo "Permuting data from " $F_READER_DIMORDER " to vxyzt"
    mri_permute -memlimit 32000000 -input ${rdr_output_fname}.mri \
            -headerout ${rdr_permute_fname}.mri -order vxyzt
    rm ${rdr_output_fname}.mri ${rdr_output_fname}.dat
else 
    mv ${rdr_output_fname}.mri ${rdr_permute_fname}.mri
    mv ${rdr_output_fname}.dat ${rdr_permute_fname}.dat
endif

if ("$F_READER_SPACE" == "i") then
    echo "Converting fdata from image space to k-space"
    parallel.run.csh recon -input ${rdr_permute_fname}.mri \
        -headerout $2.mri \
        -dataout .dat -direction forward -recon complex
    rm ${rdr_permute_fname}.mri ${rdr_permute_fname}.dat
else
    mv ${rdr_permute_fname}.mri $2.mri
    mv ${rdr_permute_fname}.dat $2.dat
endif
