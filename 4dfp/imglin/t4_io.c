/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_io.c,v 1.3 2010/03/03 01:26:02 avi Exp $*/
/*$Log: t4_io.c,v $
 * Revision 1.3  2010/03/03  01:26:02  avi
 * t4_write returns void
 *
 * Revision 1.2  2010/02/19  01:50:19  avi
 * correct bug in response to error return from t4_read()
 *
 * Revision 1.1  2009/02/23  06:14:36  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAXL	256

extern int split (char *string, char *srgv[MAXL], int maxl);
extern void errr (char *program, char *filspc);

int	t4_read1 (FILE *fp, float *t4) {
	int	i, m;

	for (i = 0; i < 4; i++) {
		m = fscanf (fp, " %f %f %f %f", t4 + i + 0,  t4 + i + 4, t4 + i + 8, t4 + i + 12);
		if (m != 4) return -1;
	}
	return 0;
}

void t4_init (float *t4) {
	int	i;

	for (i = 0; i < 16; i++)    t4[i] = 0.;
	for (i = 0; i < 16; i += 5) t4[i] = 1.;
}

int t4_read (char *t4file, float *t4) {
	FILE		*fp;
	char		string[MAXL], *srgv[MAXL];
	int 		m;
	int		done;


	if (!(fp = fopen (t4file, "r"))) {
		fprintf (stderr, "%s read error; transform initialized to I4\n", t4file);
		t4_init (t4);
		return 1;
	}
	printf ("Reading: %s\n", t4file);
	done = 0; while (fgets (string, MAXL, fp)) {
		if (!(m = split (string, srgv, MAXL))) continue;		/* skip blank lines */
		if (m == 1 && !strcmp (srgv[0], "t4")) {
			if (t4_read1 (fp, t4)) errr ("t4_read", t4file);
			done++;
		}
	}
	if (!done) {
		rewind (fp);
		if (t4_read1 (fp, t4)) errr ("t4_read", t4file);
	}
	fclose (fp);
	return 0;
}

int iget_t4_scale (char *t4file, float *scale) {
	FILE		*fp;
	char		string[MAXL], *srgv[MAXL];
	int 		m, err;

	*scale = 1.;
	err = 1;
	if (!(fp = fopen (t4file, "r"))) return -1;
	printf ("Reading: %s\n", t4file);
	while (fgets (string, MAXL, fp)) {
		if (!(m = split (string, srgv, MAXL))) continue;		/* skip blank lines */
		if (!strcmp (srgv[0], "scale:") && m > 1) {
			*scale = atof (srgv[1]);
			err = 0;
		}
	}
	fclose (fp);
	return err;
}

void	t4_list (FILE *fp, float *t4) {
	int	i, m;

	fprintf (fp, "t4\n");
	for (i = 0; i < 4; i++) {
		fprintf (fp, "%10.6f%10.6f%10.6f%10.4f\n", (t4 + i)[0], (t4 + i)[4], (t4 + i)[8], (t4 + i)[12]);
	}
}

void	t4_write (FILE *fp, float *t4, float scale) {
	int	i, m;

	t4_list (fp, t4);
	if (scale != 0) fprintf (fp, "scale:    %10.6f\n", scale);
}
