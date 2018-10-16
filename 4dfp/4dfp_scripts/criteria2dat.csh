#!/bin/csh
#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/criteria2dat.csh,v 1.3 2011/01/02 04:34:08 avi Exp $
#$Log: criteria2dat.csh,v $
# Revision 1.3  2011/01/02  04:34:08  avi
# eliminate option -o; output file root taken from input
#
# Revision 1.2  2010/01/10  03:37:22  avi
# option -o
#
# Revision 1.1  2010/01/10  03:14:04  avi
# Initial revision
#

set program = $0; set program = $program:t

#set echo
@ k = 0
@ m = 1
while ($m <= ${#argv})
	set swi = `echo $argv[$m] | awk '{print substr($1,1,2)}'`
	set arg = `echo $argv[$m] | awk '{print substr($0,3)}'`
	switch ($swi)
	case -*:
		echo option $argv[$m] not recognized
		breaksw;
	default:
		switch ($k)
		case 0:
			set F = $argv[$m];
			@ k++; breaksw;
		default:
			breaksw;
		endsw
	endsw
	@ m++
end
if (! $k) goto USAGE

if (! -r $F) then
	echo $program":	"$F not found
	exit -1
endif
set outroot = $F:r

set slines = (`awk '/significance/{printf("%d ", NR)}' $F`)
@ n = ${#slines}
echo $n "("$slines")"

@ l = `wc $F | awk '{print $1}'`
@ l++
set slines = ($slines $l)

@ k = 1
while ($k <= $n)
	set sig = `head -$slines[$k] $F | tail -1 | gawk '{sub(/=/," "); print $2;}'`
	set file = ${outroot}_significance$sig.dat
	echo "Writing: "$file
	head -$slines[$k] $F | tail -1 | gawk '{printf("#%s\n", $0);}'	>! $file
	@ kp1 = $k + 1
	@ nhead = $slines[$kp1] - 1
	@ ntail = $nhead - $slines[$k]
	head -$nhead $F | tail -$ntail					>> $file
	@ k++
end

exit

USAGE:
echo "Usage:	"$program "<criteria text file>"
echo " e.g.,	"$program "PIBpos_surr_v_PIBneg_surr_PIBneg-DAT_ROI_masked_cluster.criteria"
echo "	option"
exit 1
