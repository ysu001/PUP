#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/read_bootstrap_log.awk,v 1.6 2010/12/30 22:42:19 avi Exp $
#$Log: read_bootstrap_log.awk,v $
# Revision 1.6  2010/12/30  22:42:19  avi
# change nrep increment signal from "permutation" to "iter"
#
# Revision 1.5  2010/12/30  03:53:59  avi
# remove obligatory print "savetable="savetable;
#
# Revision 1.4  2010/07/12  06:14:04  avi
# option savetable
# sbig -> 10000
#
# Revision 1.3  2010/02/01  05:23:09  avi
# sbig -> 2500
#
# Revision 1.2  2009/10/26  05:38:40  avi
# dthresh -> 0.1
#
# Revision 1.1  2009/08/15  22:41:13  avi
# Initial revision
#
BEGIN {
	print "$Id: read_bootstrap_log.awk,v 1.6 2010/12/30 22:42:19 avi Exp $";
	freq[1] = 0.05;
	freq[2] = 0.02;
	freq[3] = 0.01;
	for (k = 4;  k  <= 6;  k++) freq[k] = freq[k - 3]/10;
	for (k = 7;  k  <= 9;  k++) freq[k] = freq[k - 3]/10;
	nfreq = 9;
	if (0) for (k = 1; k <= nfreq; k++) print k, freq[k];
	smax = 0;	# maximum encountered ROI voxel count
	ktmin = 1000;
	ktmax = 0;
	nrep = 0;
	dthresh = 0.1	# interval appropriate for t-images
	sbig = 10000;	# maximum tabulated ROI voxel count
}

($1 == "dthresh=") {dthresh = $NF;}

($1 == "iter") {nrep++;}

($1 == "thresh=") {
	thresh = $NF;
	kt = int ((thresh + 0.000001)/dthresh);
	if (kt > ktmax) ktmax = kt;
	if (kt < ktmin) ktmin = kt;
	go = 0;
}

NF == 4 && go > 0 {
	if ($2 > smax) smax = $2
	s = $2 - $2 % 5;
	kt = int ((thresh + 0.000001)/dthresh);
	a[kt,s] += 1;
#	printf ("%3d %6d %6d\n", kt, s, a[kt,s]);
}

($1 == "region") {go++;}

END {
	smax += 5; print "smax="smax;
	print "ktmin="ktmin;
	print "ktmax="ktmax;
	print "dthresh="dthresh;

	if (0) {
		for (kt = ktmin; kt <= ktmax; kt++) {
			for (s = 5; s <= smax; s += 5) {
				printf ("%3d %6d %6d\n", kt, s, a[kt,s]);
			}
		}
		exit;
	}
	for (kt = ktmin; kt <= ktmax; kt++) {
	if (0) print "kt="kt;
		for (s = 5; s <= sbig; s += 5) {
			b[kt,s] = 0;
			for (t = s; t <= smax; t += 5) b[kt,s] += a[kt,t];
			if (0) printf ("%3d %6d %6d\n", kt, s, b[kt,s]);
		}
	}

	if (savetable) {
		printf ("%10s%10s%15s\n", "threshold", "nvoxcrit", "frequency");
		for (kt = ktmin; kt <= ktmax; kt++) {
			thresh = kt*dthresh;
			for (s = 5; s <= sbig; s += 5) {
				if (b[kt,s]) printf ("%10.2f%10d%15.6f\n", thresh, s, b[kt,s]/nrep);
			}
		}
		exit;
	}

	for (ifreq = 1; ifreq <= nfreq; ifreq++) {
		printf ("significance=%.4f\n", freq[ifreq]);
		for (kt = ktmin; kt <= ktmax; kt++) {
			for (s = sbig; s > 0; s -= 5) {
				if (b[kt,s] >= freq[ifreq]*nrep) break;
			}
			if (s > 0 && s < sbig) {
#				printf ("%3d %6d %6d\n", kt, s, b[kt,s]);
				if (b[kt,s+5] > 0) {
					w = (log(b[kt,s]) - log(freq[ifreq]*nrep))/(log(b[kt,s]) - log(b[kt,s+5]));
					logscrit = w*log(s+5) + (1 - w)*log(s);
					scrit = exp(logscrit);

				} else {
					w = (b[kt,s] - freq[ifreq]*nrep)/(b[kt,s] - b[kt,s+5]);
					scrit = s + w*5;
				}
				printf ("%4.2f %6.1f\n", kt*dthresh, scrit);
			}
		}
	}
	exit;
}
