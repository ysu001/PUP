#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/eval_scrubbed_movement.awk,v 1.1 2012/12/13 05:07:15 avi Exp $
#$Log: eval_scrubbed_movement.awk,v $
# Revision 1.1  2012/12/13  05:07:15  avi
# Initial revision
#

$1!~/#/ && NF==4 {
	n++;
	for(i = 0; i < 4; i++) a[n,i] = $(i+1);
}

END {
	for(k = 1; k <= n; k++) {
		printf("%10d%10.4f%10.4f%10.4f\n", a[k,0], a[k,1], a[k,2], a[k,3]);
		ntot  += a[k,0];
		trans += a[k,0]*a[k,1]^2;
		rot   += a[k,0]*a[k,2]^2;
		tot   += a[k,0]*a[k,3]^2;
	}
	printf("%10s%10.4f%10.4f%10.4f\n", "rms", sqrt(trans/ntot), sqrt(rot/ntot), sqrt(tot/ntot));
}
