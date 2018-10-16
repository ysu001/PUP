#!/bin/csh

set idstr = '$Id: generic_RXF1_permute.csh,v 1.2 2012/06/22 04:57:06 avi Exp $'
echo $idstr
set program = $0; set program = $program:t

if (${#argv} < 1) then
        echo "Usage:    "$program" <params_file>"
        echo " e.g.:    "$program "NP317_PM_vs_AM_tbeta_pcc_N25.params"
	echo "A model params file listing follows:"
	echo "------------------------------------"
	echo "set mask	= 711-2B_mask_g5_t500_333m		# 4dfp mask of voxels to be included in the analysis"
	echo "set grp1	= NP317_PM_vs_AM_tbeta_pcc		# contrast label (will be used to name created files)"
	echo "set grp1_zfrm	= NP317_PM_vs_AM_tbeta_pcc_N25.lst	# list of 4dfp data images"
	echo "set grp1_dfnd	= NP317_dfnd_N25.lst			# list of 4dfp defined voxel images"
	echo "@ niter		= 10000					# permutation count (set to a low number to test script)"
	echo "------------------------------------"
	echo "N.B.: "$program" currently accepts only single volume result images"
	echo "N.B.: all 4dfp image specifications may include a path with extensions optional"
        exit 1
endif

set prmfile = $1
echo "prmfile="$prmfile

if (! -e $prmfile) then
        echo $program": "$prmfile not found
        exit -1
endif
source $prmfile

if ($mask:e ==  "img") set mask = $mask:r
if ($mask:e == "4dfp") set mask = $mask:r
set grp1_zfrm	= (`cat $grp1_zfrm`)
set grp1_dfnd	= (`cat $grp1_dfnd`)

set log		= ${grp1}_${mask:t}_cluster.log
if (-e $log) /bin/rm $log
touch $log
############
# make lists
############
set lst1 = ${grp1}.lst
if (-e $lst1) /bin/rm $lst1
touch $lst1
set nst1 = ${grp1}_N.lst
if (-e $nst1) /bin/rm $nst1
touch $nst1
@ n = $#grp1_zfrm
@ k = 1
while ($k <= $n)
	set F = $grp1_zfrm[$k]; if ($F:e == "img") set F = $F:r; if ($F:e == "4dfp") set F = $F:r
	echo $F	>> $lst1
	set F = $grp1_dfnd[$k]; if ($F:e == "img") set F = $F:r; if ($F:e == "4dfp") set F = $F:r
	echo $F	>> $nst1
	@ k++
end
echo "n="$n; 
echo $lst1; cat $lst1
echo $nst1; cat $nst1

imgopr_4dfp -R	-a${grp1}_N -l$nst1		> /dev/null
if ($status) exit $status

RFX_4dfp $lst1 ${grp1}_N -R			> /dev/null
if ($status) exit $status
set spmt = RFX_${lst1:r}_tmap
t2z_4dfp $spmt -n${grp1}_N
if ($status) exit $status
set spmz = ${spmt}_z
maskimg_4dfp -R $spmz $mask ${spmz}_${mask:t}
if ($status) exit $status

@ it = 20
while ($it <= 55)
	set thresh = `echo $it/10 | bc -l | awk '{printf("%.2f",$1)}'`
	echo "thresh= "$thresh							>> $log
	cluster_4dfp -At$thresh ${spmz}_${mask:t} | awk '$2 >= 5 {print}'	>> $log
	if ($status) exit $status
	@ it++
end

img_hist_4dfp -m$mask -r-5to5 -b100 -hu ${spmz}_${mask:t}

SURR:
if (! ${?niter}) exit

#####################################################
# make surrograte group list (same in all iterations)
#####################################################
set lstc = ${grp1}_surr.lst
if (-e $lstc) /bin/rm $lstc
touch $lstc
@ k = 1
while ($k <= $n)
	set F = $grp1_zfrm[$k]:t		# surrogate data will be on /tmp
	if ($F:e == "img") set F = $F:r; if ($F:e == "4dfp") set F = $F:r
	echo /tmp/${F}_surr	>> $lstc	# surrogate data filenames remain constant
	@ k++
end
echo cat $lstc
cat $lstc

set log	 = ${grp1}_surr_${mask:t}_cluster.log
if (-e $log) /bin/rm $log
touch $log
@ iter = 1
while ($iter <= $niter)
	@ m = $iter - 1
	@ m = $m * $n				# $m cumulatively increments from 0 over iterations and data
	echo iter $iter			>> $log

################################################
# copy data to /tmp and pseudorandomly flip sign
################################################
	@ k = 1
	while ($k <= $n)
		set F = $grp1_zfrm[$k]
		if ($F:e == "img") set F = $F:r; if ($F:e == "4dfp") set F = $F:r
		foreach e (img ifh)
			/bin/cp $F.4dfp.$e /tmp/${F:t}_surr.4dfp.$e
		end
		@ seed = $m + $k
		set perm = (`permuteN 2 -s$seed`)
		if ($perm[1] % 2) then
			scale_4dfp /tmp/${F:t}_surr -1	> /dev/null
			if ($status) exit $status
		endif
		@ k++
	end

###########################
# run RFX on surrogate data
###########################
	RFX_4dfp $lstc ${grp1}_N -R			> /dev/null
	if ($status) exit $status
	set spmt = RFX_${lstc:r}_tmap
	t2z_4dfp $spmt -n${grp1}_N
	if ($status) exit $status
	set spmz = ${spmt}_z
	maskimg_4dfp -R $spmz $mask ${spmz}_${mask:t}
	if ($status) exit $status

###############
# clean up /tmp
###############
	@ k = 1
	while ($k <= $n)
		set F = $grp1_zfrm[$k]
		if ($F:e == "img") set F = $F:r; if ($F:e == "4dfp") set F = $F:r
		/bin/rm /tmp/${F:t}_surr.4dfp.*
		@ k++
	end
CLUS:
	@ it = 20
	while ($it <= 55)
		set thresh = `echo $it/10 | bc -l | awk '{printf("%.2f",$1)}'`
		echo "thresh= "$thresh							>> $log
		cluster_4dfp -At$thresh ${spmz}_${mask:t} | awk '$2 >= 5 {print}'	>> $log
		@ it++
	end

HIST:
	set H = ${grp1}_surr_${mask:t}.hist
	img_hist_4dfp -m$mask -r-5to5 -b100 -h ${spmz}
	if ($iter == 1) then
		gawk '$1!~/#/{print}' ${spmz}.hist >! $H
	else
		gawk '$1!~/#/{print}' ${spmz}.hist >! $$.hist
		paste $$.hist $H >! $$1.hist
		gawk '{printf("%15.6f%20d\n", $1, $2 + $4)}' $$1.hist >! $H
	endif

	@ iter++
end

/bin/rm $$*.hist ${spmz}.hist
@ nvox = `ROI_voxcount $mask`
set v = `echo "$nvox * $niter / 10.0" | bc -l`
gawk '{printf("%15.6f%15.6f\n", $1, $2/v)}' v=$v $H >! $H:r"_norm".hist

gawk -f $RELEASE/read_bootstrap_log.awk $log >! ${log:r}.criteria
criteria2dat.csh 				${log:r}.criteria

exit
