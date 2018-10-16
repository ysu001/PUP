#!/bin/csh
#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/fcMRI_preproc_090115.csh,v 1.11 2012/11/09 23:11:56 avi Exp $
#$Log: fcMRI_preproc_090115.csh,v $
# Revision 1.11  2012/11/09  23:11:56  avi
# understand $MB flag
#
# Revision 1.10  2012/11/06  02:09:22  avi
# *** empty log message ***
#
# Revision 1.9  2012/11/06  01:03:14  avi
# params file option noGSR
#
# Revision 1.8  2011/03/15  22:26:19  avi
# more onestep definition up (prevent script failure if $conc in defined in params)
#
# Revision 1.7  2010/11/02  02:56:54  avi
# additional QA (sd1 before/after fcMRI preproc)
#
# Revision 1.6  2010/10/07  00:28:19  avi
# specify $FCdir in params
#
# Revision 1.5  2009/10/27  06:46:08  avi
# enable params file format specification to over-ride default con2format
#

##############################
# fcMRI-specific preprocessing
##############################
set program = $0
set program = $program:t
set rcsid = '$Id: fcMRI_preproc_090115.csh,v 1.11 2012/11/09 23:11:56 avi Exp $'
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

if (! ${?FCdir}) set FCdir = FCmaps
echo "##### list parameters from $prmfile #####"
echo "patid 	= $patid"
echo "srcdir	= $srcdir"	# containes bold directories
echo "workdir	= $workdir"	# contains $FCdir		
echo "fcbolds	= $fcbolds"
echo "TR_vol	= $TR_vol"
echo "skip	= $skip"
echo "FCdir 	= $FCdir"
echo "##### end parameters list #####"

@ onestep = 0 		# if set then $program will exit at the end of each step
set blur	= .735452	# = .4413/6, i.e., 6 mm blur that is conservative relative to fidl's Monte Carlo
set WB_region	= $REFDIR/glm_atlas_mask_333
set VENT_WM_reg = ($REFDIR/711-2B_lat_vent_333 $REFDIR/small_WM)
if (! -e $WB_region.4dfp.img	|| ! -e $WB_region.4dfp.ifh		\
||  ! -e $VENT_WM_reg[1].4dfp.img	|| ! -e $VENT_WM_reg[1].4dfp.ifh	\
||  ! -e $VENT_WM_reg[2].4dfp.img	|| ! -e $VENT_WM_reg[2].4dfp.ifh) then
	echo $program error: standard nuisance regressor ROI 4dfp images not found
	exit -1
endif

#set echo
if (! ${?MB}) @ MB = 0			# skip slice timing correction and debanding
set MBstr = _faln_dbnd; if ($MB) set MBstr = ""

if (! -e $FCdir) mkdir $FCdir
if (! ${?conc}) then
	set concroot	= ${patid}${MBstr}_xr3d_atl
	set conc	= $concroot.conc
else
	set concroot = $conc:r
	if (-e $conc) then
		/bin/cp $conc* $FCdir
		if ($status) exit $status
		goto COMPUTE_DEFINED
	endif
endif
####################################
# make conc file and mv it to $FCdir
####################################
CONC:
touch $$.lst
foreach run ($fcbolds)
	set file = $srcdir/bold$run/$patid"_b"$run${MBstr}"_xr3d_atl"			
	echo $file >> $$.lst
end
conc_4dfp $concroot -l$$.lst -w
if ($status) exit $status
/bin/rm $$.lst
/bin/mv $conc* $FCdir

################
# debugging code
################
if ($onestep) then
	set echo
	pushd $FCdir; goto GLM
endif

##########################
# run compute_defined_4dfp
##########################
COMPUTE_DEFINED:
pushd $FCdir		# into $FCdir
compute_defined_4dfp $conc
if ($status) exit $status
maskimg_4dfp ${concroot}_dfnd $WB_region ${concroot}_dfndm
if ($status) exit $status

#####################
# compute initial sd1
#####################
if (! ${?format}) set format = `conc2format $conc $skip`
var_4dfp -s -f$format	${concroot}.conc
ifh2hdr -r20		${concroot}_sd1
set sd1_WB0 = `qnt_4dfp ${concroot}_sd1 ${concroot}_dfndm | awk '$1~/Mean/{print $NF}'`
if ($onestep) exit

##############
# spatial blur
##############
BLUR:
gauss_4dfp $conc $blur
if ($status) exit $status
if ($onestep) exit

##########################
# temporal bandpass filter
##########################
BANDPASS:
bandpass_4dfp ${concroot}_g7.conc $TR_vol -bh.1 -oh2 -E 
if ($status) exit $status
if ($onestep) exit

############################################################
# make movement regressors for each BOLD run
# convert rdat (within-run) and ddat (differentiated) output
# of mat2dat for use as regressors in glm_4dfp
############################################################
MOVEMENT:
set regr_output = $workdir/$FCdir/$patid"_mov_regressor".dat
if (-e $regr_output) /bin/rm $regr_output; touch $regr_output

gawk '/file:/{l=index($0,"file:"); print substr($0,l+5);}' $conc > $$.2
gawk '{gsub(/_atl.4dfp.img/,""); print}' $$.2 >! $$.3
set mats = (`cat $$.3`)
set regr_output = $patid"_mov_regressor".dat
if (-e $regr_output) /bin/rm $regr_output; touch $regr_output
if (! -e $RELEASE/trendout.awk) then
	echo $program error: $RELEASE/trendout.awk not found
	exit -1
endif
foreach F ($mats)
	mat2dat $F.mat -R -D
	if ($status) exit $status
	set out = $F"_movement_regressor.dat"
	awk '$1!~/#/{for (i = 2; i <= 7; i++) printf ("%10s", $i); printf ( "\n" );}' $F.rdat >! $$.0
	awk '$1!~/#/{for (i = 2; i <= 7; i++) printf ("%10s", $i); printf ( "\n" );}' $F.ddat >! $$.1
	paste $$.0 $$.1 >! $out		
	cat $out | gawk -f $RELEASE/trendout.awk >> $regr_output
end
/bin/rm $$.*
@ nframe = `wc $patid"_mov_regressor".dat | awk '{print $1}'`

#############################################################
# make the whole brain regressor including the 1st derivative
#############################################################
WB:
if (! ${?format}) set format = `conc2format $conc $skip`
qnt_4dfp -s -d -f$format ${concroot}_g7_bpss.conc $WB_region \
	| awk '$1!~/#/{printf("%10.4f%10.4f\n", $2, $3)}' >! ${patid}_WB_regressor_dt.dat
@ n = `wc ${patid}_WB_regressor_dt.dat | awk '{print $1}'`
if ($n != $nframe) then
	echo $patid"_mov_regressor".dat ${patid}_WB_regressor_dt.dat length mismatch
	exit -1
endif

############################################################################
# make ventricle and bilateral white matter regressors and their derivatives
############################################################################
VENT_WM:
set output = ${patid}_vent_wm_dt.dat; if (-e $output) /bin/rm $output; touch $output
qnt_4dfp -s -d -f$format ${concroot}_g7_bpss.conc $VENT_WM_reg[1] \
	| awk '$1!~/#/{printf("%10.4f%10.4f\n", $2, $3)}' >! $output"1"
qnt_4dfp -s -d -f$format ${concroot}_g7_bpss.conc $VENT_WM_reg[2] \
	| awk '$1!~/#/{printf("%10.4f%10.4f\n", $2, $3)}' >! $output"2"
paste $output"1" $output"2" >! $output
/bin/rm $output*[12]
@ n = `wc ${patid}_vent_wm_dt.dat | awk '{print $1}'`
if ($n != $nframe) then
	echo $patid"_mov_regressor".dat ${patid}_vent_wm_dt.dat.dat length mismatch
	exit -1
endif

#############################################
# optional externally supplied task regressor
#############################################
TASK:
if (! ${?task_regressor}) set task_regressor = ""
if ($task_regressor != "") then
	if (! -r $task_regressor) then
		echo $task_regressor not accessible
		exit -1
	endif
	@ n = `wc $task_regressor | awk '{print $1}'`
	if ($n != $nframe) then
		echo $patid"_mov_regressor".dat $task_regressor length mismatch
		exit -1
	endif
endif

####################################
# paste nuisance regressors together
####################################
PASTE:
set WB = ${patid}_WB_regressor_dt.dat
if (${?noGSR}) then
	if ($noGSR) set WB = ""
endif
paste ${patid}_mov_regressor.dat $WB ${patid}_vent_wm_dt.dat $task_regressor >! ${patid}_nuisance_regressors.dat	

##########################################################################
# run glm_4dfp to remove nuisance regressors out of volumetric time series
##########################################################################
GLM:
if (! ${?format}) set format = `conc2format $conc $skip`
glm_4dfp $format ${patid}_nuisance_regressors.dat	${concroot}_g7_bpss.conc -rresid -o
if ($status) exit $status
ifh2hdr -r-20to20					${concroot}_g7_bpss_coeff
var_4dfp -s -f$format					${concroot}_g7_bpss_resid.conc
if ($status) exit $status
ifh2hdr -r20						${concroot}_g7_bpss_resid_sd1
set sd1_WB1 = `qnt_4dfp ${concroot}_g7_bpss_resid_sd1 ${concroot}_dfndm | awk '$1~/Mean/{print $NF}'`

echo $sd1_WB0 | awk '{printf("whole brain mean sd1 before fcMRI preprocessing = %8.4f\n",$1)}'
echo $sd1_WB1 | awk '{printf("whole brain mean sd1  after fcMRI preprocessing = %8.4f\n",$1)}'

popd		# out of $FCdir
exit
