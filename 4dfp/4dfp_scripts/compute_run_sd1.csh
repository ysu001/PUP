#!/bin/csh
#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/compute_run_sd1.csh,v 1.3 2004/08/23 23:48:07 avi Exp $
#$Log: compute_run_sd1.csh,v $
# Revision 1.3  2004/08/23  23:48:07  avi
# $patid now a required argument
#
# Revision 1.2  2004/07/31  01:12:25  avi
# run var_4dfp -sn4 on both *xr3d_norm and *xr3d_atl
# automatically make movies
#
# Revision 1.1  2003/07/08  02:21:37  avi
# Initial revision
#

set program = $0
set program = $program:t

if (${#argv} < 1) then
	echo "Usage:	"$program" <patid>"
	echo "e.g.,	"$program" VB15792"
	echo "Note:	"$program" runs var_4dfp -sn4 on both *xr3d_norm and *xr3d_atl and makes movies"
	exit 1
endif
set patid = $1
foreach d (bold*)
	if (! -d $d) continue
	pushd $d
	set files = (`ls *_xr3d_norm.4dfp.img *_xr3d_atl.4dfp.img`)
	if (${#files} < 1) goto POP
	set files = ($files:gr); set files = ($files:gr);
	foreach f ($files)
		var_4dfp -sn4 $f
		ifh2hdr $f"_sd1".4dfp.img -r50
	end
POP:
	popd
end
if (! -e var_maps) mkdir var_maps
/bin/mv bold*/*sd1* var_maps

pushd var_maps
foreach type (xr3d_norm_sd1 xr3d_atl_sd1)
	ls *$type.4dfp.img | awk '{printf ("%s\t1\n", $1)}'  >! $patid"_"$type"_movie".lst
	cat							$patid"_"$type"_movie".lst
	@ n = `wc						$patid"_"$type"_movie".lst | awk '{print $1}'`
	if ($n < 2) then
		/bin/rm						$patid"_"$type"_movie".lst
		continue
	endif
	paste_4dfp -ap1						$patid"_"$type"_movie".lst $patid"_"$type"_movie"
	ifh2hdr -r50						$patid"_"$type"_movie"
end

popd
exit
