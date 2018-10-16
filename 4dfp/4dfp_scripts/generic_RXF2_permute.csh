#!/bin/csh

set idstr = '$Id: generic_RXF2_permute.csh,v 1.2 2011/08/20 00:15:07 avi Exp $'
echo $idstr
set program = $0; set program = $program:t

if (${#argv} < 1) then
        echo "Usage:    "$program" <params_file>"
        echo " e.g.:    "$program "CDR05_vs_CDR1_pcc.params"
	echo "A model params file listing follows:"
	echo "------------------------------------"
	echo "set mask	= CAPIIO_gfc_GM_333		# 4dfp mask of voxels to be included in the analysis"
	echo "set grp1	= CDR05_pcc			# group 1 label (will be used to name created files)"
	echo "set grp2	= CDR1_pcc			# group 2 label (will be used to name created files)"
	echo "set grp1_zfrm	= pcc_CDR05_N91.lst		# list of group 1 4dfp result images"
	echo "set grp2_zfrm	= pcc_CDR1_N33.lst		# list of group 2 4dfp result images"
	echo "set grp1_dfnd	= CDR05_dfnd_N91.lst		# list of group 1 4dfp defined voxel images"
	echo "set grp2_dfnd	= CDR1_dfnd_N33.lst		# list of group 2 4dfp defined voxel images"
	echo "@ niter		= 10000				# permutation count (set to a low number to test script)"
	echo "------------------------------------"
	echo "N.B.: "$program" currently accepts only single volume result images"
	echo "N.B.: all 4dfp image specifications may include a path; extensions are optional"
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
set grp2_zfrm	= (`cat $grp2_zfrm`)
set grp1_dfnd	= (`cat $grp1_dfnd`)
set grp2_dfnd	= (`cat $grp2_dfnd`)

set log		= ${grp1}_vs_${grp2}_${mask:t}_cluster.log
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
@ n1 = $#grp1_zfrm
@ k = 1
while ($k <= $n1)
	echo $grp1_zfrm[$k]	>> $lst1
	echo $grp1_dfnd[$k]	>> $nst1
	@ k++
end
echo "n1="$n1; 
echo $lst1; cat $lst1
echo $nst1; cat $nst1

set lst2 = ${grp2}.lst
if (-e $lst2) /bin/rm $lst2
touch $lst2
set nst2 = ${grp2}_N.lst
if (-e $nst2) /bin/rm $nst2
touch $nst2
@ n2 = $#grp2_zfrm
@ k = 1
while ($k <= $n2)
	echo $grp2_zfrm[$k]	>> $lst2
	echo $grp2_dfnd[$k]	>> $nst2
	@ k++
end
echo "n2="$n2; 
echo $lst2; cat $lst2
echo $nst2; cat $nst2
@ ntot = $n1 + $n2
echo "ntot="$ntot; 

imgopr_4dfp -R	-a${grp1}_N -l$nst1		> /dev/null
if ($status) exit $status
imgopr_4dfp -R	-a${grp2}_N -l$nst2		> /dev/null
if ($status) exit $status

RFX_4dfp $lst1 ${grp1}_N $lst2 ${grp2}_N -R	> /dev/null
if ($status) exit $status
set spmt = RFX_${lst1:r}_vs_${lst2:r}_tmap
t2z_4dfp $spmt -n${grp1}_N -n${grp2}_N
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

set lstc = ${grp1}_${grp2}.lst
cat $lst1 $lst2 >! $lstc
set nstc = ${grp1}_${grp2}_N.lst
cat $nst1 $nst2 >! $nstc

set log		= ${grp1}_surr_vs_${grp2}_surr_${mask:t}_cluster.log
if (-e $log) /bin/rm $log
touch $log
@ iter = 1
while ($iter <= $niter)
############################
# make surrograte goup lists
############################
	set list = (`permuteN $ntot -s$iter`)
	echo iter $iter			>> $log
	echo "permutation "$list	>> $log

	set lst1 = ${grp1}_surr.lst
	if (-e $lst1) /bin/rm $lst1
	touch $lst1
	set nst1 = ${grp1}_surr_N.lst
	if (-e $nst1) /bin/rm $nst1
	touch $nst1
	@ k = 1
	while ($k <= $n1)
		@ j = $list[$k]
		head -$j $lstc | tail -1	>> $lst1
		head -$j $nstc | tail -1	>> $nst1
		@ k++
	end

	set lst2 = ${grp2}_surr.lst
	if (-e $lst2) /bin/rm $lst2
	touch $lst2
	set nst2 = ${grp2}_surr_N.lst
	if (-e $nst2) /bin/rm $nst2
	touch $nst2
	while ($k <= $ntot)
		@ j = $list[$k]
		head -$j $lstc | tail -1	>> $lst2
		head -$j $nstc | tail -1	>> $nst2
		@ k++
	end

	imgopr_4dfp -R	-a${grp1}_surr_N -l$nst1		> /dev/null
	if ($status) exit $status
	imgopr_4dfp -R	-a${grp2}_surr_N -l$nst2		> /dev/null
	if ($status) exit $status

	RFX_4dfp $lst1 ${grp1}_surr_N $lst2 ${grp2}_surr_N -R	> /dev/null
	if ($status) exit $status
	set spmt = RFX_${lst1:r}_vs_${lst2:r}_tmap
	t2z_4dfp $spmt -n${grp1}_surr_N -n${grp2}_surr_N
	if ($status) exit $status
	set spmz = ${spmt}_z
	maskimg_4dfp -R $spmz $mask ${spmz}_${mask:t}
	if ($status) exit $status

CLUS:
	@ it = 20
	while ($it <= 55)
		set thresh = `echo $it/10 | bc -l | awk '{printf("%.2f",$1)}'`
		echo "thresh= "$thresh							>> $log
		cluster_4dfp -At$thresh ${spmz}_${mask:t} | awk '$2 >= 5 {print}'	>> $log
		@ it++
	end

HIST:
	set H = ${grp1}_surr_vs_${grp2}_surr_${mask:t}.hist
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
