#!/bin/csh
set idstr = '$Id: bval_bvec_to_prm.csh,v 1.1 2013/01/16 22:31:06 avi Exp $'
echo $idstr
set program = $0; set program = $program:t

if (${#argv} < 1) then
        echo "Usage:    "$program" <bval|bvec root> [outroot]"
        echo " e.g.:    "$program "x20120905_15301513466705891138085520115936WIPDTI32dirTE108TE12300TimBarretts701a1007 BCH32dir1400"
        exit 1
endif
set inroot = $1
if (! -e $inroot.bval) then
	echo $inroot.bval not found
	exit -1
endif
if (! -e $inroot.bvec) then
	echo $inroot.bvec not found
	exit -1
endif
if (${#argv} > 1) then
	set outroot = $2;
else
	set outroot = $inroot:t;
endif
set bvals = (`cat $inroot.bval`)
@ n = $#bvals

echo "#"$program ${argv[1-]} $USER `date`	>! $outroot.prm
echo "#"$idstr					>> $outroot.prm
echo $n						>> $outroot.prm
set bvecx = (`head -1 $inroot.bvec | tail -1`)
set bvecy = (`head -2 $inroot.bvec | tail -1`)
set bvecz = (`head -3 $inroot.bvec | tail -1`)
@ k = 1
while ($k <= $n)
	echo $bvals[$k] $bvecx[$k] $bvecy[$k] $bvecz[$k] | \
		awk '{if ($1 > 0) printf("%-10.1f%10.6f%10.6f%10.6f\n", $1, $2, $3, -$4); else printf("0.0\n");}' \
		>> $outroot.prm
	@ k++
end
cat $outroot.prm

exit 0
