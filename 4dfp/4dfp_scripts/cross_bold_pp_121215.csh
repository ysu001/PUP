#!/bin/csh
#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/cross_bold_pp_121215.csh,v 1.1 2012/12/18 05:09:09 avi Exp $
#$Log: cross_bold_pp_121215.csh,v $
# Revision 1.1  2012/12/18  05:09:09  avi
# Initial revision
#

set idstr = '$Id: cross_bold_pp_121215.csh,v 1.1 2012/12/18 05:09:09 avi Exp $'
echo $idstr
set program = $0; set program = $program:t

if (${#argv} < 1) then
	echo "usage:	"$program" <params_file> [instructions]"
	echo " e.g.:	"$program" APsub1.params ../butterfly.params"
	exit 1
endif
set prmfile = $1
echo "prmfile="$prmfile

if (! -e $prmfile) then
	echo $program": "$prmfile not found
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

if ($target:h != $target) then
	set tarstr = -T$target
else
	set tarstr = $target
endif

@ runs = ${#irun}
if ($runs != ${#fstd}) then
	echo "irun fstd mismatch - edit "$prmfile
	exit -1
endif

if (! ${?scrdir}) set scrdir = ""
@ usescr = `echo $scrdir | awk '{print length ($1)}'`
if ($usescr) then 
	if (! -e $scrdir) mkdir $scrdir
	if ($status) exit $status
endif
set sourcedir = $cwd

date
###################
# process BOLD data
###################
set interleave = ""
if (${?Siemens_interleave}) then
	if ($Siemens_interleave) set interleave = "-N"
endif
if (! ${?to_MNI152}) @ to_MNI152 = 0
if (! ${?day1_patid}) set day1_patid = "";

@ err = 0
@ k = 1
while ($k <= $runs)
	if ($usescr) then		# test to see if user requested use of scratch disk
		if (-e bold$irun[$k]) /bin/rm bold$irun[$k]	# remove existing link
		if (! -d $scrdir/bold$irun[$k]) mkdir $scrdir/bold$irun[$k]
		ln -s $scrdir/bold$irun[$k] bold$irun[$k]
	else
		if (! -d bold$irun[$k]) mkdir bold$irun[$k]
	endif
	pushd bold$irun[$k]
	set y = $patid"_b"$irun[$k]"_faln_dbnd"
	if (-e $y.4dfp.img && -e $y.4dfp.ifh) goto POP
	if ($sorted) then
		echo		dcm_to_4dfp -q -b study$fstd[$k] $inpath/study$fstd[$k]
		if ($go)	dcm_to_4dfp -q -b study$fstd[$k] $inpath/study$fstd[$k]
	else
		echo		dcm_to_4dfp -q -b study$fstd[$k] $inpath/$dcmroot.$fstd[$k]."*"
		if ($go)	dcm_to_4dfp -q -b study$fstd[$k] $inpath/$dcmroot.$fstd[$k].*
	endif
	echo		unpack_4dfp -V study$fstd[$k] $patid"_b"$irun[$k] -nx$nx -ny$ny
	if ($go)	unpack_4dfp -V study$fstd[$k] $patid"_b"$irun[$k] -nx$nx -ny$ny
	if ($status) then
		@ err++
		/bin/rm $patid"_b"$irun[$k]*
		goto POP
	endif
	echo		/bin/rm  study$fstd[$k]."*"
	if ($go)	/bin/rm  study$fstd[$k].*

	echo		frame_align_4dfp $patid"_b"$irun[$k] $skip -TR_vol $TR_vol -TR_slc $TR_slc -d $epidir $interleave
	if ($go)	frame_align_4dfp $patid"_b"$irun[$k] $skip -TR_vol $TR_vol -TR_slc $TR_slc -d $epidir $interleave

	echo		deband_4dfp -n$skip $patid"_b"$irun[$k]"_faln"
	if ($go)	deband_4dfp -n$skip $patid"_b"$irun[$k]"_faln"
	if ($status)	exit $status

	if ($economy > 2) then
		echo		/bin/rm $patid"_b"$irun[$k].4dfp."*"
		if ($go)	/bin/rm $patid"_b"$irun[$k].4dfp.*
	endif
	if ($economy > 3) then
		echo		/bin/rm $patid"_b"$irun[$k]"_faln".4dfp."*"
		if ($go)	/bin/rm $patid"_b"$irun[$k]"_faln".4dfp.*
	endif
POP:
	popd	# out of bold$irun[$k]
	@ k++
end
if ($err) then
	echo $program": one or more BOLD runs failed preliminary processing"
	exit -1
endif
if ($epi2atl == 2) goto ATL

if (-e  $patid"_xr3d".lst) /bin/rm $patid"_xr3d".lst; touch $patid"_xr3d".lst
if (-e  $patid"_anat".lst) /bin/rm $patid"_anat".lst; touch $patid"_anat".lst
@ k = 1
while ($k <= $runs)
	echo bold$irun[$k]/$patid"_b"$irun[$k]"_faln_dbnd" >>			$patid"_xr3d".lst
	echo bold$irun[$k]/$patid"_b"$irun[$k]"_faln_dbnd_xr3d_norm" 1 >>	$patid"_anat".lst
	@ k++
end

echo cat	$patid"_xr3d".lst
cat		$patid"_xr3d".lst
echo		cross_realign3d_4dfp -n$skip -qv$normode -l$patid"_xr3d".lst
if ($go)	cross_realign3d_4dfp -n$skip -qv$normode -l$patid"_xr3d".lst
if ($status)	exit $status

date
#################################
# compute mode 1000 normalization
#################################
@ k = 1
while ($k <= $runs)
	pushd bold$irun[$k]
	echo 		normalize_4dfp $patid"_b"$irun[$k]"_faln_dbnd_r3d_avg" -h
	if ($go)	normalize_4dfp $patid"_b"$irun[$k]"_faln_dbnd_r3d_avg" -h
	if ($economy > 4 && $epi2atl == 0) then
		echo		/bin/rm $patid"_b"$irun[$k]"_faln_dbnd".4dfp."*"
		if ($go)	/bin/rm $patid"_b"$irun[$k]"_faln_dbnd".4dfp.*
	endif
	popd	# out of bold$irun[$k]
	@ k++
end

date
###############################
# apply mode 1000 normalization
###############################
@ k = 1
while ($k <= $runs)
	pushd bold$irun[$k]
	set file = $patid"_b"$irun[$k]"_faln_dbnd_r3d_avg_norm".4dfp.img.rec
	set f = 1.0; if (-e $file) set f = `head $file | awk '/original/{print 1000/$NF}'`
	echo		scale_4dfp $patid"_b"$irun[$k]"_faln_dbnd_xr3d" $f -anorm
	if ($go)	scale_4dfp $patid"_b"$irun[$k]"_faln_dbnd_xr3d" $f -anorm
	echo		/bin/rm $patid"_b"$irun[$k]"_faln_dbnd_xr3d".4dfp."*"
	if ($go)	/bin/rm $patid"_b"$irun[$k]"_faln_dbnd_xr3d".4dfp.*
	popd	# out of bold$irun[$k]
	@ k++
end

date
###################
# movement analysis
###################
if (! -d movement) mkdir movement
@ k = 1
while ($k <= $runs)
	echo		mat2dat bold$irun[$k]/"*_xr3d".mat -RD -n$skip
	if ($go)	mat2dat bold$irun[$k]/*"_xr3d".mat -RD -n$skip
	echo		/bin/mv bold$irun[$k]/"*_xr3d.*dat"	movement
	if ($go)	/bin/mv bold$irun[$k]/*"_xr3d".*dat	movement
	@ k++
end

date
######################################
# make EPI first frame (anatomy) image
######################################
echo cat	$patid"_anat".lst
cat		$patid"_anat".lst
echo		paste_4dfp -p1 $patid"_anat".lst	$patid"_anat_ave"
if ($go)	paste_4dfp -p1 $patid"_anat".lst	$patid"_anat_ave"
echo		ifh2hdr	-r2000				$patid"_anat_ave"
if ($go)	ifh2hdr	-r2000				$patid"_anat_ave"

######################
# atlas transformation
######################
if (! -d atlas) mkdir atlas
echo		mv $patid"_anat*" atlas
if ($go)	mv $patid"_anat"* atlas

pushd atlas
######################
# make MP-RAGE average
######################
if ($day1_patid != "") then
	echo		cross_day_imgreg_4dfp $patid $day1_path $day1_patid $tarstr
	if ($go)	cross_day_imgreg_4dfp $patid $day1_path $day1_patid $tarstr
	if ($status) exit $status
	goto EPI_to_ATL
endif

@ nmpr = ${#mprs}
if ($nmpr < 1) exit 0
set mprave = $patid"_mpr_n"$nmpr
set mprlst = ()
@ k = 1
while ($k <= $nmpr)
	if ($sorted) then
		echo		dcm_to_4dfp -b $patid"_mpr"$k $inpath/study$mprs[$k]
		if ($go)	dcm_to_4dfp -b $patid"_mpr"$k $inpath/study$mprs[$k]
	else
		echo		dcm_to_4dfp -b $patid"_mpr"$k $inpath/$dcmroot.$mprs[$k]."*"
		if ($go)	dcm_to_4dfp -b $patid"_mpr"$k $inpath/$dcmroot.$mprs[$k].*
	endif
	if ($status) exit $status
	set mprlst = ($mprlst $patid"_mpr"$k)
	@ k++
end
echo		avgmpr_4dfp $mprlst $mprave $tarstr useold 
if ($go)	avgmpr_4dfp $mprlst $mprave $tarstr useold
if ($status) exit $status

date
#########################
# compute atlas transform
#########################
if (! ${?tse})	set tse = ()
if (! ${?t1w})	set t1w = ()
if (! ${?pdt2})	set pdt2 = ()

echo		avgmpr_4dfp $mprlst $mprave $tarstr useold
if ($go)	avgmpr_4dfp $mprlst $mprave $tarstr useold
if ($status) exit $status
set episcript = epi2t2w2mpr2atl2_4dfp;

@ ntse = ${#tse}
if ($ntse) then
	set tselst = ()
	@ k = 1
	while ($k <= $ntse)
		set filenam = $patid"_t2w"
		if ($ntse > 1) set filenam = $filenam$k
		if ($sorted) then
			echo		dcm_to_4dfp -b $filenam $inpath/study$tse[$k]
			if ($go)	dcm_to_4dfp -b $filenam $inpath/study$tse[$k]
		else
			echo		dcm_to_4dfp -b $filenam $inpath/$dcmroot.$tse[$k]."*"
			if ($go) 	dcm_to_4dfp -b $filenam $inpath/$dcmroot.$tse[$k].*
		if ($status) exit $status
		endif
		set tselst = ($tselst $filenam)
		@ k++
	end
	if ($ntse  > 1) then
		echo		collate_slice_4dfp $tselst $patid"_t2w"
		if ($go)	collate_slice_4dfp $tselst $patid"_t2w"
	endif
	echo		$episcript $patid"_anat_ave" $patid"_t2w" $patid"_mpr1" useold $tarstr
	if ($go)	$episcript $patid"_anat_ave" $patid"_t2w" $patid"_mpr1" useold $tarstr
else if (${#t1w}) then
	if ($sorted) then
		echo		dcm_to_4dfp -b $patid"_t1w" $inpath/study$t1w
		if ($go)	dcm_to_4dfp -b $patid"_t1w" $inpath/study$t1w
	else
		echo		dcm_to_4dfp -b $patid"_t1w" $inpath/$dcmroot.$t1w."*"
		if ($go) 	dcm_to_4dfp -b $patid"_t1w" $inpath/$dcmroot.$t1w.*
	endif
	echo		t2w2mpr_4dfp $patid"_mpr1" $patid"_t1w" $tarstr
	if ($go)	t2w2mpr_4dfp $patid"_mpr1" $patid"_t1w" $tarstr

	echo		epi2t1w_4dfp $patid"_anat_ave" $patid"_t1w" $tarstr
	if ($go)	epi2t1w_4dfp $patid"_anat_ave" $patid"_t1w" $tarstr
	if ($status) exit $status

	echo		t4_mul $patid"_anat_ave_to_"$patid"_t1w_t4" $patid"_t1w_to_"$target:t"_t4"
	if ($go)	t4_mul $patid"_anat_ave_to_"$patid"_t1w_t4" $patid"_t1w_to_"$target:t"_t4"
else if (${#pdt2}) then
	if ($sorted) then
		echo		dcm_to_4dfp -b $patid"_pdt2" $inpath/study$pdt2
		if ($go)	dcm_to_4dfp -b $patid"_pdt2" $inpath/study$pdt2
	else
		echo		dcm_to_4dfp -b $patid"_pdt2" $inpath/$dcmroot.$pdt2."*"
		if ($go) 	dcm_to_4dfp -b $patid"_pdt2" $inpath/$dcmroot.$pdt2.*
	endif
	if ($status) exit $status
	echo		extract_frame_4dfp $patid"_pdt2" 2 -o$patid"_t2w"
	if ($go)	extract_frame_4dfp $patid"_pdt2" 2 -o$patid"_t2w"
	if ($status) exit $status

	echo		$episcript $patid"_anat_ave" $patid"_t2w" $patid"_mpr1" useold $tarstr
	if ($go)	$episcript $patid"_anat_ave" $patid"_t2w" $patid"_mpr1" useold $tarstr
else
	echo		epi2mpr2atl2_4dfp $patid"_anat_ave" $patid"_mpr1" useold $tarstr
	if ($go)	epi2mpr2atl2_4dfp $patid"_anat_ave" $patid"_mpr1" useold $tarstr
endif
if ($status) exit $status

EPI_to_ATL:
if ($to_MNI152) then
	foreach x (*to_${target:t}_t4)
		set out = `echo $x | gawk '{gsub(T,"MNI152lin"); print;}' T=${target:t}`
		t4_mul $x $RELEASE/711-2B_to_MNI152lin_T1_t4 $out
	end
	set A = $RELEASE/MNI152_T1_1mm.4dfp.ifh
	foreach img (anat_ave mpr1 t1w t2w)
		if (-e ${patid}_${img}_to_MNI152lin_t4) then
			t4img_4dfp ${patid}_${img}_to_MNI152lin_t4 ${patid}_${img} ${patid}_${img}_on_MNI152lin_1mm -O$A
		endif
	end
else
##################################################################
# make atlas transformed EPI anat_ave in images 711-2B atlas space
##################################################################
	set t4file = ${patid}_anat_ave_to_${target:t}_t4
	foreach O (111 222 333)
		echo		t4img_4dfp $t4file  $patid"_anat_ave"	$patid"_anat_ave_t88_"$O -O$O
		if ($go)	t4img_4dfp $t4file  $patid"_anat_ave"	$patid"_anat_ave_t88_"$O -O$O
		echo		ifh2hdr	 -r2000				$patid"_anat_ave_t88_"$O
		if ($go)	ifh2hdr	 -r2000				$patid"_anat_ave_t88_"$O
	end
endif
if ($status) exit $status

/bin/rm *t4% >& /dev/null
popd		# out of atlas

if (! $epi2atl) exit 0
ATL:
date
###################################################################
# make cross-realigned atlas-transformed resampled BOLD 4dfp stacks
###################################################################
@ k = 1
while ($k <= $runs)
	pushd bold$irun[$k]
	set file = ${patid}_b$irun[$k]_faln_dbnd_r3d_avg_norm.4dfp.img.rec
	set f = 1.0; if (-e $file) set f = `head $file | awk '/original/{print 1000/$NF}'`
	if ($to_MNI152) then
		set A = $RELEASE/MNI152_T1_3mm.4dfp.ifh
		echo		t4_xr3d_4dfp $sourcedir/atlas/${patid}_anat_ave_to_MNI152lin_t4	  ${patid}_b$irun[$k]_faln_dbnd -axr3d_MNI152_3mm -v$normode -c$f -O$A
		if ($go)	t4_xr3d_4dfp $sourcedir/atlas/${patid}_anat_ave_to_MNI152lin_t4   ${patid}_b$irun[$k]_faln_dbnd -axr3d_MNI152_3mm -v$normode -c$f -O$A
	else
		echo		t4_xr3d_4dfp $sourcedir/atlas/${patid}_anat_ave_to_${target:t}_t4 ${patid}_b$irun[$k]_faln_dbnd -axr3d_atl        -v$normode -c$f -O333
		if ($go)	t4_xr3d_4dfp $sourcedir/atlas/${patid}_anat_ave_to_${target:t}_t4 ${patid}_b$irun[$k]_faln_dbnd -axr3d_atl        -v$normode -c$f -O333
	endif
	if ($status) exit $status
	if ($economy > 6) then
		echo		/bin/rm ${patid}_b$irun[$k]_faln_dbnd_xr3d_norm.4dfp."*"
		if ($go)	/bin/rm ${patid}_b$irun[$k]_faln_dbnd_xr3d_norm.4dfp.*
	endif
	popd	# out of bold$irun[$k]
	@ k++
end

echo $program done status=$status
exit
