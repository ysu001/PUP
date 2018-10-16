/*$Header: /home/usr/shimonyj/diff4dfp/RCS/t4tolin.c,v 1.2 2008/09/16 02:58:55 avi Exp $*/
/*$Log: t4tolin.c,v $
 * Revision 1.2  2008/09/16  02:58:55  avi
 * #include <JSSutil.h>
 *
 * Revision 1.1  2005/12/08  03:13:59  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <JSSutil.h>

#define	MAXL	256

extern void errr (char* program, char* filespc);	/* endianio.c */

static char rcsid[] = "$Id: t4tolin.c,v 1.2 2008/09/16 02:58:55 avi Exp $";
void	t4tolin (char *t4file, int orient, float **lin) {
	FILE	*matfp;
	char	string[MAXL];
	float	**inv;
	int	a = 1, b = 1, c = 1;

/*******************************************/
/* set orientation specific flip constants */
/*******************************************/
	switch (orient) {
		case 2: a = -1; b = +1; c = -1; break;	/* transverse */
		case 3: a = -1; b = +1; c = +1; break;	/* coronal */
		case 4: a = +1; b = +1; c = +1; break;	/* sagittal */
		default: fprintf (stderr, "t4tolin: illegal 4dfp orientation (%d)\n", orient);
		exit (-1);
	}

	inv = matrix (1, 3, 1, 3);
/****************/
/* read t4 file */
/****************/
	fprintf (stdout, "Reading: %s\n", t4file);
	if (!(matfp = fopen (t4file, "r"))) errr ("t4tolin", t4file);

/**********************************/
/* read until first field is "t4" */
/**********************************/
	while (fgets (string, MAXL, matfp)) {
		if (0) printf ("line:\t%s", string);
		if (!(strncmp (string, "t4", 2))) break;
	}
	inv = matrix (1, 3, 1, 3);
	if (fscanf (matfp, "%f %f %f %*f", &inv[1][1], &inv[1][2], &inv[1][3]) != 3) {
		rewind (matfp);
		if (fscanf (matfp, "%f %f %f %*f", &inv[1][1], &inv[1][2], &inv[1][3]) != 3) errr ("t4tolin", t4file);
	}
	if (fscanf (matfp, "%f %f %f %*f", &inv[2][1], &inv[2][2], &inv[2][3]) != 3) errr ("t4tolin", t4file);
	if (fscanf (matfp, "%f %f %f %*f", &inv[3][1], &inv[3][2], &inv[3][3]) != 3) errr ("t4tolin", t4file);
	fclose (matfp);

	inv_svd (inv, 3, 3, lin);
	free_matrix (inv, 1, 3, 1, 3);
	
	lin[1][2] *= a*b;
	lin[1][3] *= a*c;
	lin[2][1] *= b*a;
	lin[2][3] *= b*c;
	lin[3][1] *= c*a;
	lin[3][2] *= c*b;
	return;
}
