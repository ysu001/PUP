#!/bin/csh
##############################################
# single-subject multi-ROI correlation mapping
##############################################
#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/seed_correl_090115.csh,v 1.15 2012/11/10 04:46:12 avi Exp $
#$Log: seed_correl_090115.csh,v $
# Revision 1.15  2012/11/10  04:46:12  avi
# read $MB flag
#
# Revision 1.14  2011/01/24  04:26:30  avi
# suppress computation of ROI;ROI CCR is $nROI > 256
#
# Revision 1.13  2011/01/19  02:45:12  avi
# optional params file variable ROIimg (alternative to both ROIlist and ROIlistfile)
#
# Revision 1.12  2010/11/02  02:14:12  avi
# always compute ROI-ROI correlation matrix
#
# Revision 1.11  2010/08/06  03:17:01  avi
# remove $ROIs before creating file
#
# Revision 1.10  2010/05/20  04:03:40  avi
# correct operation when $ROIlist not defined
#
# Revision 1.9  2010/05/20  01:29:24  avi
# replace use of qnt_4dfp (and paste) with qntm_4dfp logic
#
# Revision 1.8  2010/04/21  03:02:03  avi
# FCdir optionally specified in params
#
# Revision 1.7  2010/04/18  23:54:16  avi
# unset echo
#
# Revision 1.6  2010/04/17  01:17:43  avi
# ROIlistfile capability
# correct dfnd volume addressing error introduced in last revision
#
# Revision 1.5  2010/04/14  04:51:24  avi
# correct another error when $conc is defined
#
# Revision 1.4  2010/04/14  04:36:17  avi
# correct bug in case $conc is defined
#
# Revision 1.3  2009/10/27  06:54:58  avi
# enable params file format specification to over-ride con2format default
#
# Revision 1.2  2009/05/17  03:55:00  avi
# eliminate wildcards from maskimg_4dfp command
#
# Revision 1.1  2009/01/16  04:17:50  avi
# Initial revision
#
set program = $0
set program = $program:t
set rcsid = '$Id: seed_correl_090115.csh,v 1.15 2012/11/10 04:46:12 avi Exp $'
echo $rcsid

if (${#argv} < 1) then
	echo "Usage:	$program <parameters file> [instructions]"
	echo "e.g.,	$program VB16168.params"
	exit 1
endif 
date
uname -a

set prmfile = $1
if (! -e $prmfile) then
	echo $prmfile not found
	exit -1
endif
source $prmfile
if (${#argv} > 1) then
	set instructions = $2
	if (! -e $instructions) then
		echo $program": "$instructions not found
		exit -1
	endif
	cat $instructions
	source $instructions
endif
if (! ${?MB}) @ MB = 0			# skip slice timing correction and debanding
set MBstr = _faln_dbnd; if ($MB) set MBstr = ""

if (! ${?FCdir}) set FCdir = FCmaps
echo "FCdir 	= $FCdir"
echo "patid 	= $patid"
echo "skip	= $skip"
echo "ROIdir	= $ROIdir"; 		if ($status) exit $status

###########
# ROI image
###########
if (${?ROIimg}) then
	set ROIsrc = $ROIdir/$ROIimg
	if ($ROIsrc:e == "img")  set ROIsrc = $ROIsrc:r
	if ($ROIsrc:e == "4dfp") set ROIsrc = $ROIsrc:r
	@ nROI = `awk '/matrix size \[4\]/{print $NF}' $ROIsrc.4dfp.ifh`
else
	if (${?ROIlistfile}) then
		@ nROI = `wc $ROIlistfile | awk '{print $1}'`
	else
		echo "ROIlist = $ROIlist"; 	if ($status) exit $status
		@ nROI = ${#ROIlist}
	endif
	set ROIs = $cwd/$$ROIs.lst; if (-e $ROIs) /bin/rm $ROIs; touch $ROIs
	@ k = 1
	while ($k <= $nROI)
		if (${?ROIlistfile}) then
			set ROI = `head -$k $ROIlistfile | tail -1`
			echo $ROIdir/$ROI		>> $ROIs
		else
			echo $ROIdir/$ROIlist[$k]	>> $ROIs
		endif
		@ k++
	end
	paste_4dfp -ap1 $ROIs $$ROIs
	if ($status) exit $status
	set ROIsrc = $cwd/$$ROIs
endif

pushd $FCdir		# into $FCdir
if ($status) exit $status

##########################################################
# identify fcMRI preprocessed conc (assume $cwd is $patid)
##########################################################
if (! ${?conc}) then
	set concroot	= ${patid}${MBstr}_xr3d_atl
	set conc	= $concroot.conc
else
	if (! -e $conc) then
		echo $program error: $conc not found
		exit -1
	endif
	set concroot = $conc:r
endif
set residroot	= ${concroot}_g7_bpss_resid
set resid	= $residroot.conc

set lst	= ${patid}_seed_regressors.dat
set rec = ${patid}_seed_regressors.rec; if (-e $rec) /bin/rm $rec;
#####
# rec
#####
if (${?ROIimg}) then
	ln -s $ROIsrc.4dfp.img.rec $rec
else
	touch $rec
	@ k = 1
	while ($k <= $nROI)
		set region = `head -$k $ROIs | tail -1`
		if ($region:e == "img")  set region = $region:r
		if ($region:e == "4dfp") set region = $region:r
		echo Volume $k"	"$region >> $rec
		@ k++
	end
endif

###########
# qntm_4dfp
###########
qntm_4dfp -V $resid $ROIsrc -o$lst
if ($status) exit $status

###########################################################
# compute total correlations and update _tcorr.4dfp.img.rec
###########################################################
if (! ${?format}) set format = `conc2format $conc $skip`
glm_4dfp $format $lst $resid -t
if ($status) exit $status
set l = `wc -l	${residroot}_tcorr.4dfp.img.rec | awk '{print $1}'`
@ l--
head -$l	${residroot}_tcorr.4dfp.img.rec	>! $$.rec
cat $rec					>> $$.rec
tail -1		${residroot}_tcorr.4dfp.img.rec	>> $$.rec
/bin/mv $$.rec	${residroot}_tcorr.4dfp.img.rec

##############################
# mask tcorr by defined voxels
##############################
maskimg_4dfp	${residroot}_tcorr ${concroot}_dfnd -1 ${residroot}_tcorr_dfnd
if ($status) exit $status
ifh2hdr -r-2to2	${residroot}_tcorr_dfnd

######################################
# Fisher z transform total correlation
######################################
rho2z_4dfp	${residroot}_tcorr_dfnd
if ($status) exit $status
ifh2hdr -r-2to2	${residroot}_tcorr_dfnd_zfrm

#############################################################
# compute zero-lag ROI-ROI correlation matrix if $nROI <= 256
#############################################################
if ($nROI <= 256) then
	covariance -uom0 $format $lst 
	if ($status) exit $status
	/bin/rm ${lst:r}_ROI*_CCR.dat
endif

popd			# out of $FCdir
##########
# clean up
##########
if (! ${?ROIimg}) /bin/rm $$ROIs*
exit
