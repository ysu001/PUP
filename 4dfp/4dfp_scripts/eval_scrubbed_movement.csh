#!/bin/csh
#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/eval_scrubbed_movement.csh,v 1.2 2012/12/14 00:27:30 avi Exp $
#$Log: eval_scrubbed_movement.csh,v $
# Revision 1.2  2012/12/14  00:27:30  avi
# option -D (evaluate ddat results instead or rdat results)
#
# Revision 1.1  2012/12/13  05:07:29  avi
# Initial revision
#
set rcsid = '$Id: eval_scrubbed_movement.csh,v 1.2 2012/12/14 00:27:30 avi Exp $'
echo $rcsid
set program = $0; set program = $program:t

if (${#argv} < 1) goto USAGE

@ debug = 0
set fMRIdir = .
set radius = 50.
@ ddat_flag = 0;
@ m = 1
@ k = 0
while ($m <= ${#argv})
	set swi = `echo $argv[$m] | awk '{print substr($1,1,2)}'`
	set arg = `echo $argv[$m] | awk '{print substr($0,3)}'`
	switch ($swi)
	case -d:
		@ debug++;			breaksw;
	case -D:
		@ ddat_flag++;			breaksw;
	case -H*:
		set set fMRIdir = $arg;		breaksw;
	case -r*:
		set set radius = $arg;		breaksw;
	default:
		switch ($k)
		case 0:
			set format = $argv[$m];
			@ k++; breaksw;
		case 1:
			set lst = $argv[$m];
			@ k++; breaksw;
		default:
			breaksw;
		endsw
	endsw
	@ m++
end
if ($k < 2) goto USAGE

echo "fMRI run path = "$fMRIdir
if (! -e $lst) then
	echo $lst not found
	exit 1
endif
echo "format = "$format

@ nrun = `wc $lst | awk '{printf $1}'`
echo $nrun runs

set string = `format2lst -e $format`
#echo "format = "$string
#condense $string; echo

set out = $lst:t;
if ($ddat_flag) then
	set out = ${out:r}_scrubbed_ddat_rms
	echo | gawk '{printf ("#%9s%10s%10s%10s (rms dmm or ddeg) at %.2f mm\n", "frames", "dtrans", "drot", "dtot", radius)}' radius=$radius >! $out
else
	set out = ${out:r}_scrubbed_rms
	echo | gawk '{printf ("#%9s%10s%10s%10s (rms mm or deg) at %.2f mm\n", "frames", "trans", "rot", "tot", radius)}' radius=$radius >! $out
endif

@ k = 1
@ l = 0
while ($k <= $nrun)
	set file = `head -$k $lst | tail -1`
	if ($file:e == "img")  set file = $file:r
	if ($file:e == "4dfp") set file = $file:r
	set matfile = $fMRIdir/${file}_xr3d.mat
	if (! -e $matfile) then
		echo $matfile not found
		exit -1
	endif
	@ nframe = `cat $matfile | gawk '/s4 frame/{n = $NF}; END {print n}'`
	echo $matfile has $nframe frames
#	echo $nframe $l
	set fmt = `echo $string | gawk '{print substr($1, beg + 1, len)}' beg=$l len=$nframe`
	if ($fmt == "") break
	set fmt = `condense $fmt`
echo	mat2dat -RDL -f$fmt -r$radius $matfile
	mat2dat -RDL -f$fmt -r$radius $matfile | tail -6
	set datfile = ${matfile:t}; set datfile = ${datfile:r}.dat; set rdatfile = ${datfile:r}.rdat; set ddatfile = ${datfile:r}.ddat;
	if ($ddat_flag) then
		set F = $ddatfile;
	else
		set F = $rdatfile;
	endif
	gawk '/counting/{count = $2};/translation/{trans = $NF};/rotation/{rot = $NF};/total/{tot = $7};END{printf("%10d%10.4f%10.4f%10.4f\n", count, trans, rot, tot);}' $F >> $out
	if (! $debug) /bin/rm $datfile $rdatfile $ddatfile
	@ l += $nframe
	@ k++
end
cat $out
gawk -f $RELEASE/eval_scrubbed_movement.awk $out | tail -1 >! ${out}_allruns

exit 0

USAGE:
set format = "5x38+2x65+8(5x105+)"
echo "Usage:	"$program" <format> <xr3d_listfile>"
echo " e.g.:	"$program \"${format}\" "VB16168_xr3d.lst -D/data/nil-bluearc/raichle/avi/NP645/VB16168_TEST"
echo "	option"
echo "	-d	debug mode (preserve local dat and rdat files)"
echo "	-H<dir>	specify directory containing BOLD subdirectories"
echo " 	-r<flt> specify head radius in mm for total motion computation (default=50mm)"
exit 1
