#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/blockdesign_to_regressor.awk,v 1.3 2006/01/06 23:17:57 avi Exp $
#$Log: blockdesign_to_regressor.awk,v $
# Revision 1.3  2006/01/06  23:17:57  avi
# read block off time from event file fourth column
#
# Revision 1.2  2005/12/29  02:55:59  avi
# better usage
#
# Revision 1.1  2005/12/29  00:13:17  avi
# Initial revision
#
# fidl event files DO NOT explicitly code the OFF event.
# OFF events occur ontim frames after explicit ON events

BEGIN {
	nevent = 0;	# event counter
#	Rshift = 1;	# HDR right shift relative to event
#	ntrans = 5;	# number of transient frames
#	nframe = 4*116;	# total frames in dataset
}

NR == 1 {
	ntype = NF - 1;	# block type count
}

NR > 1 && NF >= 3 {
	frame[nevent] = $1;
	itype[nevent] = $2;
	ontim[nevent] = $3;
	offtm[nevent] = $4;
	nevent++;
}

END {
	if (nframe == 0) {
		printf ("Usage:\tnawk -f blockdesign_to_regressor.awk nframe=<int> [Rshift=<int>] [ntrans=<int>] <my_event_file>\n");
		printf (" e.g.,\tnawk -f blockdesign_to_regressor.awk nframe=464 Rshift=1 ntrans=5 jf27Feb03.event_file\n");
		printf ("\tVariables:\n");
		printf ("	nframe= total number of frames in conc dataset\n");
		printf ("	Rshift= HDR shift (in frames) relative to events\n");
		printf ("	ntrans= number of frames over which to model block onset/offset transients\n");
		printf ("N.B.:	input event file must be fidl block design type\n");
		printf ("N.B.:	2*(ntrans + 1) columns are output for every event type\n");
		printf ("N.B.:	non-initialized variables default to zero\n");
		printf ("N.B.:	block off times are taken from column 4\n");
		printf ("N.B.:	default block off times are computed but wrong at run transitions\n");
	}
	ncol = ntype*2*(ntrans + 1);	# ntype*(ON+OFF)*(transient + sustained)

	frame[nevent] = nframe;
###################
# compute off times
###################
	for (i = 0; i < nevent; i++) if (offtm[i] < 1) offtm[i] = frame[i + 1] - frame[i] - ontim[i];

##########################
# initialize design matrix
##########################
	for (i = 0; i < nframe; i++) for (icol = 0; icol < ncol; icol++) {
		f[icol,i] = 0;
	}

############
# debug code
############
	for (i = 0; i <= nevent; i++) {
		if (0) printf ("%10d%10d%10d%10d\n", frame[i], itype[i], ontim[i], offtm[i]);
	}

###########################
# fill in the design matrix
###########################
	for (ievent = 0; ievent < nevent; ievent++) {
		icol = 2*(ntrans + 1)*itype[ievent];
		i = frame[ievent] + Rshift;
##############
# on transient
##############
		for (itrans = 0; itrans < ntrans; itrans++) f[icol++,i++] = 1;
##############
# on sustained
##############
		for (ii = 0; ii < ontim[ievent] - ntrans; ii++) f[icol,i++] = 1;
		icol++;
###############
# off transient
###############
		for (itrans = 0; itrans < ntrans; itrans++) f[icol++,i++] = 1;
###############
# off sustained
###############
		for (ii = 0; ii < offtm[ievent] - ntrans; ii++) f[icol,i++] = -1;
		icol++;
	}

######################
# output design matrix
######################
	for (iframe = 0; iframe < nframe; iframe++) {
		for (icol = 0; icol < ncol; icol++) printf ("%4d", f[icol,iframe]);
		printf ("\n");
	}
}
