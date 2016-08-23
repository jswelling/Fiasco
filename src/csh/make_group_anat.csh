#! /bin/csh -ef
#
# This script makes group anatomy by averaging the anatomy for the
# individuals in the group.
#
echo '#'`date` $0
echo '# $Id: make_group_anat.csh,v 1.4 2007/04/23 19:10:29 welling Exp $'

# Check for help command
if ( $#argv >= 1 ) then 
  if ( junk$argv[1] == junk-help ) then
    if ( $#argv >= 2 ) then
      scripthelp $0 $argv[2]
    else
      scripthelp $0
    endif
    exit
  endif
endif

# Parse command line options
set args = `getopt m:d:h $*`
set help = ""
set mapstuff = ""

while ($#args > 1) 
  switch ( $args[1] )
    case '-m' : 
      if ( ${?F_MAP_PATH} ) then
        setenv F_MAP_PATH $args[2]':'${F_MAP_PATH}
      else
        setenv F_MAP_PATH $args[2]
      endif
      shift args; shift args; breaksw;
    case '-d' :
      set mapstuff = "$mapstuff -d "$args[2]
      shift args; shift args; breaksw;
    case '-h' :
      set help = 1 ; shift args; breaksw;
  endsw
  if (junk${args[1]} == 'junk--') then
    shift args;
    break;
  endif
end
if ($#args < 1 || dummy$help != dummy ) then
  scripthelp $0 usage
  exit -1
endif

#
# Make some scratch space, and go there.
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/make_group_anat_tmp_$$
else
  set tmpdir = ./make_group_anat_tmp_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
echo "Scratch directory will be " $tmpdir
set homedir = `pwd`
echo "Current directory is " $homedir
cd $tmpdir

foreach ds_type ( Axial Inplane StrippedAxial StrippedInplane )
  echo '######## Collecting' $ds_type 'anatomical datasets ########'
  set ds_list = ()

  foreach subj ($args)
    set axial_name = `map_name.py $mapstuff -d cond=1 -d file=$ds_type $subj`
    if ( -e ${axial_name}.mri ) then
      echo "  Found $ds_type anatomical for " $subj
      set ds_list= ( $ds_list $axial_name )
    else
      echo "  $ds_type anatomical for $subj not found- skipping this subject!"
    endif
  end

  if ( ${#ds_list} != 0 ) then
  
    # Build the mean for this dataset
    echo '######## Assembling the mean' $ds_type 'dataset ########'
    set count = 0
    while ( $#ds_list )
      if ( $#ds_list > 5 ) then
        set these_args = ( $ds_list[1-5] )
        set ds_list = ( $ds_list[6-] )
      else
        set these_args = ( $ds_list )
        set ds_list = ( )
      endif
      @ count = $count + $#these_args
      switch ( ${#these_args} )
        case 1 :
          mri_copy_dataset $these_args[1] tmp_sum; 
          breaksw;
        case 2 :
          mri_rpn_math -out tmp_sum '$1,$2,+' $these_args ; 
          breaksw;
        case 3 :
          mri_rpn_math -out tmp_sum '$1,$2,+,$3,+' $these_args ; 
          breaksw;
        case 4 :
          mri_rpn_math -out tmp_sum '$1,$2,+,$3,+,$4,+' $these_args ; 
          breaksw;
        case 5 :
          mri_rpn_math -out tmp_sum '$1,$2,+,$3,+,$4,+,$5,+' $these_args ; 
          breaksw;
      endsw
      if ( -e sum.mri ) then
        mri_rpn_math -out tmp2_sum '$1,$2,+' sum tmp_sum
        mri_copy_dataset tmp2_sum sum
      else
        mri_copy_dataset tmp_sum sum
      endif
      mri_rpn_math -out mean '$1,'$count',/' sum
    end
    mri_destroy_dataset sum
    mri_copy_dataset mean $homedir/mean_${ds_type}

  endif

end

# Clean up
echo '######## Cleaning up #######'
cd $homedir
rm -r $tmpdir

# And that's it.
echo '######## Done! #######'
