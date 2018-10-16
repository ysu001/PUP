/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/stretch_out.c,v 1.1 2010/03/03 01:27:10 avi Exp $*/
/*$Log: stretch_out.c,v $
 * Revision 1.1  2010/03/03  01:27:10  avi
 * Initial revision
 **/
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
void	stretchout_ (float *t4);	/* stretchout.f */

/***********/
/* globals */
/***********/
static char	program[MAXL];
static char	rcsid[] = "$Id: stretch_out.c,v 1.1 2010/03/03 01:27:10 avi Exp $";

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

void errw (char* program, char* filespc) {
	fprintf (stderr, "%s: %s write error\n", program, filespc);
	exit (-1);
}

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n#%s\n", rcsid);
}

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
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

int main (int argc, char *argv[]) {
	FILE	*outfp;
	char	t4file[MAXL], outfile[MAXL] = "";
	float	t4[16], scale;

/***********/
/* utility */
/***********/
	char		*str, string[MAXL], command[MAXL];
	int 		c, i, k, m;

/*********/
/* flags */
/*********/
	int		debug = 0;
	int		verbose = 0;

	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 'd': debug++;		break;
				case 'v': verbose++;		break;
			}
		} else switch (k) {
			case 0: strcpy (t4file,  argv[i]);	k++; break;
			case 1: strcpy (outfile, argv[i]);	k++; break;
		}
	}
	if (k < 1) {
		fprintf (stderr, "%s\n", rcsid);
		fprintf (stderr, "Usage:	%s <t4file> [t4file_new]\n", program);
		fprintf (stderr, "N.B.:	default output filename is <t4file>\"r\"\n");
		fprintf (stderr, "\toption\n");
		exit (1);
	}

	if (t4_read       (t4file, t4)) errr (program, t4file);
	if (iget_t4_scale (t4file, &scale)) scale = 0.;

	stretchout_ (t4);

	if (!strlen (outfile)) sprintf (outfile, "%sr", t4file);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);
	write_command_line (outfp, argc, argv);
	t4_write (outfp, t4, scale);
	if (fclose (outfp)) errw (program, outfile);

	exit (0);
}
