#!/bin/csh

echo "NILSRC="$NILSRC
echo "RELEASE="$RELEASE
cd $NILSRC
if (! -e $RELEASE) mkdir $RELEASE

set FC = "gcc -O -ffixed-line-length-132 -fcray-pointer -c"

set echo

pushd librms
make -f librms.mak
popd

pushd imgblur_4dfp
$FC fimgblur.f
popd

pushd img_hist_4dfp
cc -O -c fimg_mode.c
popd

pushd t4imgs_4dfp
cc -O -c cvrtflip.c
popd

pushd flip_4dfp
cc -O -c cflip.c -I${NILSRC}/TRX
popd

pushd parzen_4dfp
$FC imgpac.f
$FC mskpac.f
$FC parzen.f
popd

pushd t4_actmapf_4dfp
$FC ft4ixyz.f
popd

pushd imglin
cc -O -c t4scale.c
cc -O -c dnormal.c
$FC t4_sub.f
make -f t4_inv.mak			release
if ($status) exit $status
make -f t4_mul.mak			release
if ($status) exit $status
make -f t4_ident.mak			release
if ($status) exit $status
ln -s $RELEASE/t4_ident $RELEASE/t4_null
make -f t4_factor.mak			release
if ($status) exit $status
make -f t4_opr.mak			release
if ($status) exit $status
make -f t4_pts.mak			release
if ($status) exit $status
make -f stretch_out.mak			release
if ($status) exit $status
popd

pushd TRX
make -f endian_4dfp.mak			release
if ($status) exit $status
make -f ifh2hdr.mak			release
if ($status) exit $status
make -f asciito4dfp.mak			release
if ($status) exit $status
make -f spline23dvgh_test.mak
$FC fomega.f
$FC etadeta.f
$FC fimgeta.f
$FC fimgetae.f
$FC twoecal.f
popd

pushd JSSutil
make -f JSSutil.mak
popd

pushd actmapf_4dfp
make -f actmapf_4dfp.mak		release
if ($status) exit $status
make -f compute_defined_4dfp.mak	release
if ($status) exit $status
make -f conc_test.mak			release
if ($status) exit $status
make -f condense.mak			release
if ($status) exit $status
make -f covariance.mak			release
if ($status) exit $status
make -f covariance_analysis.mak		release
if ($status) exit $status
make -f format2lst.mak			release
if ($status) exit $status
make -f glm_4dfp.mak			release
if ($status) exit $status
make -f GC_4dfp.mak			release
if ($status) exit $status
make -f GC_dat.mak			release
if ($status) exit $status
make -f mat_algebra.mak			release
if ($status) exit $status
make -f LOWESS.mak			release
if ($status) exit $status
popd

pushd 4dfptoanalyze
make -f 4dfptoanalyze.mak		release
if ($status) exit $status
popd

pushd S2T_4dfp
make -f S2T_4dfp.mak			release
if ($status) exit $status
make -f C2T_4dfp.mak			release
if ($status) exit $status
make -f T2C_4dfp.mak			release
if ($status) exit $status
make -f T2S_4dfp.mak			release
if ($status) exit $status
popd

pushd algebra_4dfp
make -f imgopr_4dfp.mak			release
if ($status) exit $status
make -f ROI_resolve_4dfp.mak		release
if ($status) exit $status
make -f spatial_corr_4dfp.mak		release
if ($status) exit $status
make -f spatial_cov_multivol_4dfp.mak	release
if ($status) exit $status
make -f RFX_4dfp.mak			release
if ($status) exit $status
popd

pushd analyzeto4dfp
make -f analyzeto4dfp.mak		release
if ($status) exit $status
make -f hdr2txt.mak			release
if ($status) exit $status
popd

pushd cluster_4dfp
make -f cluster_4dfp.mak		release
if ($status) exit $status
popd

pushd collate_slice_4dfp
make -f collate_slice_4dfp.mak		release
if ($status) exit $status
popd

pushd crop_4dfp
make -f crop_4dfp.mak			release
if ($status) exit $status
make -f zero_slice_4dfp.mak		release
if ($status) exit $status
popd

pushd cross_realign3d_4dfp
make -f cross_realign3d_4dfp.mak	release
if ($status) exit $status
make -f mat2dat.mak			release
if ($status) exit $status
popd

pushd cs2ap_4dfp
make -f cs2ap_4dfp.mak			release
if ($status) exit $status
popd

pushd diff4dfp
make -f diff_4dfp.mak			release
if ($status) exit $status
make -f diffRGB_4dfp.mak		release
if ($status) exit $status
make -f whisker_4dfp.mak		release
if ($status) exit $status
popd

pushd dwi_xalign_4dfp
make -f dwi_xalign3d_4dfp.mak		release
if ($status) exit $status
make -f dwi_cross_xalign3d_4dfp.mak	release
if ($status) exit $status
popd

pushd dcm_dump_file
make -f dcm_dump_file.mak		release
if ($status) exit $status
popd

pushd deband_4dfp
make -f deband_4dfp.mak			release
if ($status) exit $status
popd

pushd flip_4dfp
make -f flip_4dfp.mak			release
if ($status) exit $status
popd

pushd frame_align
make -f frame_align_4dfp.mak		release
if ($status) exit $status
popd

pushd gauss_4dfp
make -f gauss_4dfp.mak			release
if ($status) exit $status
make -f hsphere_4dfp.mak		release
if ($status) exit $status
popd

pushd frame_align
make -f frame_align_4dfp.mak		release
if ($status) exit $status
popd

pushd imato4dfp2
make -f unpack_4dfp.mak			release
if ($status) exit $status
popd

pushd img2msk_4dfp
make -f img2msk_4dfp.mak		release
if ($status) exit $status
popd

pushd imgblur_4dfp
make -f imgblur_4dfp.mak		release
if ($status) exit $status
popd

pushd imgreg_4dfp
make -f imgreg_4dfp.mak			release
if ($status) exit $status
popd

pushd img_hist_4dfp
make -f img_hist_4dfp.mak		release
if ($status) exit $status
popd

pushd imgmax_4dfp
make -f imgmax_4dfp.mak			release
if ($status) exit $status
popd

pushd interp_4dfp
make -f interp_4dfp.mak			release
if ($status) exit $status
make -f bandpass_4dfp.mak		release
if ($status) exit $status
popd

pushd maskimg_4dfp
make -f maskimg_4dfp.mak		release
if ($status) exit $status
popd

pushd mpetto4dfp
make -f mpetto4dfp.mak			release
if ($status) exit $status
popd

pushd mprtir
make -f 2Dhist.mak			release
if ($status) exit $status
make -f partitiond_gfc_4dfp.mak		release
if ($status) exit $status
popd

pushd nifti_4dfp
make -f nifti_4dfp.mak			release
if ($status) exit $status
popd

pushd nonlin
make -f morphc_4dfp.mak			release
if ($status) exit $status
popd

pushd normalize_4dfp
make -f normalize_4dfp.mak		release
if ($status) exit $status
popd

pushd permute
make -f permuteN.mak			release
if ($status) exit $status
popd

pushd paste_4dfp
make -f paste_4dfp.mak			release
if ($status) exit $status
popd

pushd qnt_4dfp
make -f qnt_4dfp.mak			release
if ($status) exit $status
make -f qntm_4dfp.mak			release
if ($status) exit $status
make -f qntv_4dfp.mak			release
if ($status) exit $status
make -f qntw_4dfp.mak			release
if ($status) exit $status
popd

pushd reindex_4dfp
make -f reindex_4dfp.mak		release
if ($status) exit $status
popd

pushd rmspike_4dfp
make -f rmspike_4dfp.mak		release
if ($status) exit $status
popd

pushd scale_4dfp
make -f scale_4dfp.mak			release
if ($status) exit $status
popd

pushd sqrt_4dfp
make -f sqrt_4dfp.mak			release
if ($status) exit $status
make -f rho2z_4dfp.mak			release
if ($status) exit $status
make -f t2z_4dfp.mak			release
if ($status) exit $status
make -f z2logp_4dfp.mak			release
if ($status) exit $status
popd

pushd t4imgs_4dfp
make -f t4imgs_4dfp.mak			release
if ($status) exit $status
popd

pushd peak_4dfp				# must follow make -f t4imgs_4dfp.mak
make -f read_4dfp.mak			release
if ($status) exit $status
make -f peak_4dfp.mak			release
if ($status) exit $status
make -f index2atl.mak			release
if ($status) exit $status
make -f burn_sphere_4dfp.mak		release
if ($status) exit $status
popd

pushd imglin				# must follow make -f t4imgs_4dfp.mak
make -f t4_resolve.mak			release
if ($status) exit $status
make -f imgsurf_4dfp.mak		release
if ($status) exit $status
popd

pushd t4_actmapf_4dfp
make -f t4_xr3d_4dfp.mak		release
if ($status) exit $status
popd

pushd var_4dfp
make -f var_4dfp.mak			release
if ($status) exit $status
make -f dvar_4dfp.mak			release
if ($status) exit $status
popd

pushd wrpsmg_4dfp			# must follow make -f t4imgs_4dfp.mak
make -f wrpsmg_4dfp.mak			release
if ($status) exit $status
popd

pushd zero_lt_4dfp
make -f zero_lt_4dfp.mak		release
if ($status) exit $status
make -f zero_gt_4dfp.mak		release
if ($status) exit $status
make -f zero_ltgt_4dfp.mak		release
if ($status) exit $status
make -f zero_gtlt_4dfp.mak		release
if ($status) exit $status
popd

pushd aff_conv
make -f aff_conv.mak			release
if ($status) exit $status
popd

pushd dcm_pet
make -f dcm_pet.mak			release
if ($status) exit $status
popd

pushd sum_pet_4dfp
make -f sum_pet_4dfp.mak			release
if ($status) exit $status
popd

#############
# dcm_to_4dfp
#############
wget --help > /dev/null	# test for existence of wget executable
if (! $status) then	# if wget exists download latest version of source code
	wget ftp://ftp.nrg.wustl.edu/pub/dcm_to_4dfp/dcm_to_4dfp.tar.gz
	if ($status) exit $status
	/bin/rm dcm_to_4dfp.tar
	gzip -d dcm_to_4dfp.tar.gz
	if ($status) exit $status
	/bin/rm -rf dcm_to_4dfp
	tar xvf dcm_to_4dfp.tar
endif
pushd dcm_to_4dfp
if ($status) exit $status
./configure --with-TRX=$NILSRC/TRX --exec-prefix=$RELEASE
make
if ($status) exit $status
/bin/mv dcm_to_4dfp $RELEASE
popd

###############
# 4dfp scripts
###############
chmod 555 $NILSRC/4dfp_scripts/*
cp $NILSRC/4dfp_scripts/* $RELEASE


echo successful 4dfp suite make complete
