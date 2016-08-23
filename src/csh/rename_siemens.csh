#! /bin/csh -ef
#
# This script gets run in a directory containing a bunch of DICOM
# files with names of the form n-n-n.dcm (for 3 different n's), where
# the first n is the sequence number, the second is the image number,
# and the third is the overall slice number.  Links are created with
# names of the form file_n_dddd.dcm where n is the sequence number and
# dddd is the (0-filled) image number.
#
foreach seqnum ( 0 1 2 3 4 5 6 7 8 9 )
  if ( -e ${seqnum}-1-1.dcm ) then
    set nfiles = `ls ${seqnum}-*-*.dcm | wc -l ` 
    echo Found a total of $nfiles files for sequence $seqnum
    if ( $nfiles > 9999 ) then
      echo $nfiles 'is too many for me to deal with!'
      exit -1
    endif
    set time = 1
    while ( $nfiles > 0 )
      set nonomatch
      foreach pad ( '' '0' '00' '000' failure )
        if ( $pad == failure ) then
          echo Unable to find files for sequence $seqnum time $time
          exit -1
        endif
        set files_this_time = ( ${seqnum}-${pad}${time}-*.dcm )
        if ( "$files_this_time" != ${seqnum}-${pad}${time}-\*.dcm ) break
      end
      set num_slices = ${#files_this_time}
      if ( $time < 10 ) then
        set timestring = "000"$time
      else if ( $time < 100 ) then
        set timestring = "00"$time
      else if ( $time < 1000 ) then
        set timestring = "0"$time
      else
        set timestring = $time
      endif  
      if ( $num_slices == 1 ) then
        set newname = "file_"${seqnum}_${timestring}.dcm
	set oldname = ${files_this_time}
	if ( ! -e $oldname ) then
          echo Cannot find file $oldname !
          exit -1
        endif
        ln -s $oldname $newname
      else
        set slice = 1
        while ( $slice <= $num_slices )
          if ( $slice < 10 ) then
            set slicestring = "000"$slice
          else if ( $slice < 100 ) then
            set slicestring = "00"$slice
          else if ( $slice < 1000 ) then
            set slicestring = "0"$slice
          else
            set slicestring = $slice
          endif
	  set oldname = ${seqnum}-${time}-${slice}.dcm
	  set newname = "file_"${seqnum}_${timestring}_${slicestring}.dcm
	  if ( ! -e $oldname ) then
            echo Cannot find file $oldname !
            exit -1
          endif
          ln -s $oldname $newname
          @ slice = $slice + 1
        end
      endif
      @ nfiles = $nfiles - $num_slices
      @ time = $time + 1
    end
  endif
end
