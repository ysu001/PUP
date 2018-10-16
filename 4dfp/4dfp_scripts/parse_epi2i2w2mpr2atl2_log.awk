#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/parse_epi2i2w2mpr2atl2_log.awk,v 1.1 2009/08/21 00:08:35 avi Exp $
#$Log: parse_epi2i2w2mpr2atl2_log.awk,v $
# Revision 1.1  2009/08/21  00:08:35  avi
# Initial revision
#
BEGIN {
	m = 0;
}

/^Reading image:/ {
	if (m == 0) {
		if ($NF == src) dump();
		tar = $NF;
		m++;
	} else {
		src = $NF;
		m = 0;
	}
}

/^hessian/ {det = $3; cndnum = $6;}

/eta,q/ {eta = $2;}

/search radius/ {nr_radius = NR; np = 0;}
NR == nr_radius + 1 {
	for (k = 1; k <= 6; k++) radius[k] = $k;
	np =+ 6;
}
NR == nr_radius + 2 && $1 !~/hessian/ {
	for (k = 1; k <= NF; k++) radius[6+k] = $k;
	np += NF;
}

/^parameters/   {nr_params = NR;}
NR == nr_params + 1 {
	for (k = 1; k <= 6; k++) params[k] = $k;
}
NR == nr_params + 2 && $1 !~/eta/ {
	for (k = 1; k <= NF; k++) params[6+k] = $k;
}

END {
	dump();
}

function dump() {
	printf ("%s<-%s search radius\n", tar, src)
	for (k = 1; k <= np; k++) printf ("%10.4f", radius[k])
	printf ("\n");
	printf ("parameters\n")
	for (k = 1; k <= np; k++) printf ("%10.4f", params[k])
	printf ("\n");
	printf ("condition number = %.2f\n", cndnum);
	printf ("eta = %.5f\n", eta);
}
