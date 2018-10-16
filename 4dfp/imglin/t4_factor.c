/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_factor.c,v 1.1 2009/02/23 06:15:22 avi Exp $*/
/*$Log*/
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <t4_io.h>

#define MAXL	256

/*************/
/* externals */
/*************/
void	warp2param_ (float *t4, float *disp, float *ang, float *stretch, float *strang);	/* param12opr.f */
void	bparam18_   (float *par);								/* param12opr.f */

/***********/
/* globals */
/***********/
static char	program[MAXL];
static char	rcsid[] = "$Id: t4_factor.c,v 1.1 2009/02/23 06:15:22 avi Exp $";

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	if (ptr = strchr (string, '#')) *ptr = '\0';
	i = m = 0;
	while (m < maxp) {
		while (!isgraph ((int) string[i]) && string[i]) i++;
		if (!string[i]) break;
		srgv[m++] = string + i;
		while (isgraph ((int) string[i])) i++;
		if (!string[i]) break;
		string[i++] = '\0';
	}
	return m;
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

int main (int argc, char *argv[]) {
	char	t4file[MAXL];
	float	t4[16], d2r;
/*	equivalence (disp,par(1)),(ang,par(4)),(stretch,par(13)),(strang,par(16)) */
	struct	{
		float disp[3];
		float ang[3];
		float unused[6];
		float stretch[3];
		float strang[3];
	} par;

/***********/
/* utility */
/***********/
	char		*str, command[MAXL];
	int 		c, i, k;

/*********/
/* flags */
/*********/
	int		debug = 0;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 'd': debug++;		break;
			}
		} else switch (k) {
			case 0: strcpy (t4file, argv[i]);	k++; break;
			default:				k++; break;
		}
	}
	if (k < 1) {
		printf ("Usage:	%s <t4file>\n", program);
		printf (" e.g.,	%s vm11b_anat_ave_to_vm11b_234-3_t4\n", program);
		if (0) printf ("\toption\n");
		exit (1);
	}

	d2r = 45./atan(1.);
	if (t4_read (t4file, t4)) exit (-1);
	warp2param_ (t4, par.disp, par.ang, par.stretch, par.strang);
	bparam18_ (par.disp);
	if (fabs (par.stretch[0]) < 1.e-6 && fabs (par.stretch[1]) < 1.e-6) par.strang[2] = 0.;
	if (fabs (par.stretch[0]) < 1.e-6 && fabs (par.stretch[2]) < 1.e-6) par.strang[1] = 0.;
	if (fabs (par.stretch[1]) < 1.e-6 && fabs (par.stretch[2]) < 1.e-6) par.strang[0] = 0.;

	printf ("displacement (mm):  %10.4f%10.4f%10.4f\n", par.disp[0], par.disp[1], par.disp[2]);
	printf ("rotation (deg):     %10.4f%10.4f%10.4f\n", d2r*par.ang[0], d2r*par.ang[1], d2r*par.ang[2]);
	printf ("stretch (factor):   %10.6f%10.6f%10.6f\n", exp(par.stretch[0]), exp(par.stretch[1]), exp(par.stretch[2]));
	printf ("stretch ang (deg):  %10.4f%10.4f%10.4f\n", d2r*par.strang[0], d2r*par.strang[1], d2r*par.strang[2]);

	exit (0);
}
