#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/RFX2.csh,v 1.3 2010/04/16 03:11:30 avi Exp $
#RFX
#Initial implemention M. Fox, J. L. Vincent, 2005
#$Log: RFX2.csh,v $
# Revision 1.3  2010/04/16  03:11:30  avi
# correct bug preventing operation in single group mode
#
# Revision 1.2  2009/08/10  22:07:08  avi
# awk -> gawk
#
# Revision 1.1  2009/08/07  02:34:47  avi
# Initial revision
#
set program = $0
set program = $program:t
set rcsid = '$Id: RFX2.csh,v 1.3 2010/04/16 03:11:30 avi Exp $'
echo $rcsid

@ twosample = 0
set recstr = ""
@ debug = 0
@ k = 0
@ i = 1
while ($i <= ${#argv})
	set swi = `echo $argv[$i] | gawk '{print substr($1,1,2)}'`
	set arg = `echo $argv[$i] | gawk '{print substr($0,3)}'`
	switch ($swi)
		case -d:
			@ debug++;		breaksw;
		case -R:
			set recstr = "-R";	breaksw;
		default:
		switch ($k)
			case 0:
				set listfile1	= $argv[$i];	@ k++; breaksw;
			case 1:
				set Nimage1	= $argv[$i];	@ k++; breaksw;
			case 2:
				set listfile2	= $argv[$i];	@ k++; breaksw;
			case 3:
				set Nimage2	= $argv[$i];	@ k++; breaksw;
			default:
				breaksw;
		endsw
	endsw
	@ i++
end
if ($k < 2 || $k > 4) goto USAGE
if ($k == 4) @ twosample++

	set lf1 = $listfile1:t; set lf1 = $lf1:r
	echo $lf1 $Nimage1
if ($twosample) then
	set lf2 = $listfile2:t; set lf2 = $lf2:r
	echo $lf2 $Nimage2
endif

#######
# mean1
#######
imgopr_4dfp -u$recstr -eRFX_mean_temp1 -l$listfile1
if ($status) exit $status

###########
# variance1
###########
imgopr_4dfp -u$recstr -vRFX_var_temp1  -l$listfile1
if ($status) exit $status

#############
# var_over_N1
#############
imgopr_4dfp -rRFX_var_over_N_temp1 RFX_var_temp1 $Nimage1 -u

if ($twosample) then
	imgopr_4dfp -u$recstr -eRFX_mean_temp2 -l$listfile2
	if ($status) exit $status
	imgopr_4dfp -u$recstr -vRFX_var_temp2  -l$listfile2
	imgopr_4dfp -u$recstr -rRFX_var_over_N_temp2 RFX_var_temp2 $Nimage2
	if ($status) exit $status

#############################
# compute difference of means
#############################
	imgopr_4dfp -sRFX_num_temp RFX_mean_temp1 RFX_mean_temp2 -u

###################################
# sum standard errors in quadrature
###################################
	imgopr_4dfp -aRFX_add_var_over_N_temp RFX_var_over_N_temp1 RFX_var_over_N_temp2 -u
	sqrt_4dfp RFX_add_var_over_N_temp RFX_denom_temp

################################
# evalute Welch's approximate t'
################################
	imgopr_4dfp -rRFX_$lf1"_vs_"$lf2"_tmap" RFX_num_temp RFX_denom_temp -u
    	t2z_4dfp RFX_$lf1"_vs_"$lf2"_tmap" -n$Nimage1 -n$Nimage2
else
######################
# evalute one sample t
######################
	sqrt_4dfp RFX_var_over_N_temp1 RFX_denom_temp
	imgopr_4dfp -rRFX_$lf1"_tmap" RFX_mean_temp1 RFX_denom_temp -u
	t2z_4dfp RFX_$lf1"_tmap" -n$Nimage1
endif	

if (! $debug) /bin/rm RFX*temp*
exit 0

USAGE:
echo "Usage:	$program <list_group1> <Nimage_group1> [<list_group2> <Nimage_group2>]"
echo "	options"
echo "	-d	debug mode"
echo "	-R	suppress creation of large rec files (bootstrap mode)"
echo "N.B.:	<list_group[12]> name 4dfp images on which to run the t-test."
echo "N.B.:	<Nimage_group[12]> are 4dfp 'n' images (number of subjects for which each voxel is defined)."
echo "N.B.:	If one group is entered a t-test will be run on this group against the null hypothesis of 0."
echo "N.B.:	If two groups are entered a t-test will be run comparing the two groups and"
echo "	the computed statistic is Welch's approximate t' (Eqn. 8.11, p. 129 in Zar.)"  
exit 1
