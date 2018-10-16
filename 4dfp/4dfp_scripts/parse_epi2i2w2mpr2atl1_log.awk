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

/search radius/ {nr_radius = NR;}
NR == nr_radius + 1 {
	for (k = 1; k <= 6; k++) radius[k] = $k;
}
NR == nr_radius + 2 {
	for (k = 1; k <= 3; k++) radius[6+k] = $k;
}

/^parameters/   {nr_params = NR;}
NR == nr_params + 1 {
	for (k = 1; k <= 6; k++) params[k] = $k;
}
NR == nr_params + 2 {
	for (k = 1; k <= 3; k++) params[6+k] = $k;
}

END {
	dump();
}

function dump() {
	printf ("%s<-%s search radius\n", tar, src)
	for (k = 1; k <= 9; k++) printf ("%10.4f", radius[k])
	printf ("\n");
	printf ("parameters\n")
	for (k = 1; k <= 9; k++) printf ("%10.4f", params[k])
	printf ("\n");
	printf ("condition number = %.2f\n", cndnum);
	printf ("eta = %.5f\n", eta);
}
