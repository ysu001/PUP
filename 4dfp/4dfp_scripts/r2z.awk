#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/r2z.awk,v 1.3 2013/01/31 23:28:17 avi Exp $
#$Log: r2z.awk,v $
# Revision 1.3  2013/01/31  23:28:17  avi
# report off-diagonal r==1 but keep going
#
# Revision 1.2  2013/01/28  04:06:07  avi
# trap off-diagonal r==1.0
#
# Revision 1.1  2013/01/24  23:23:21  avi
# Initial revision
#

BEGIN {
	nrow = 0;
	ncol = 0;
	err = 0;
	offdiag = 0;
	debug = 1;
}

NR == 1 {
	printf ("#$Id: r2z.awk,v 1.3 2013/01/31 23:28:17 avi Exp $\n");
	printf ("#%s\n", FILENAME);
}

$1!~/#/ {
	if (ncol > 0 && NF != ncol) {
		err++;
		exit;
	}
	ncol = NF;
	for (j = 1; j <= ncol; j++) {
		r[nrow,j-1] = $j;
		if (r[nrow,j-1] > .9999 && nrow != j - 1 && debug) {
			printf ("#row=%3d col=%3d r=%.6f\n", nrow + 1, j, r[nrow,j-1]);
			offdiag++;
		}
	}
	nrow++;
}

END {
	if (err || nrow != ncol) {
		print "format error";
		exit -1;
	}
#	if (offdiag) exit -1;
	for (irow = 0; irow < nrow; irow++) {
		for (icol = 0; icol < ncol; icol++) {
			q = r[irow,icol];
			if (q == 1.0) {
				printf ("%10s", "Inf");
			} else {
				z = 0.5*log ((1 + q)/(1 - q));
				printf ("%10.6f", z);
			}
		}
		printf ("\n");
	}
}
