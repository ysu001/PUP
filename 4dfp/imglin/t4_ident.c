/*$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4_ident.c,v 1.1 2007/05/01 02:00:12 avi Exp $*/
/*$Log: t4_ident.c,v $
 * Revision 1.1  2007/05/01  02:00:12  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAXL	256

static char	program[MAXL];
static char	rcsid[] = "$Id: t4_ident.c,v 1.1 2007/05/01 02:00:12 avi Exp $";

void errw (char* program, char* filespc) {
	fprintf (stderr, "%s: %s write error\n", program, filespc);
	exit (-1);
}

void	t4list (FILE *fp, float *t4) {
	int	i, m;

	fprintf (fp, "t4\n");
	for (i = 0; i < 4; i++) {
		fprintf (fp, "%10.6f%10.6f%10.6f%10.4f\n", (t4 + i)[0], (t4 + i)[4], (t4 + i)[8], (t4 + i)[12]);
	}
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

int main (int argc, char *argv[]) {
	FILE		*fp;
	char		t4file[1][MAXL];
	float		t4[1][16];

/***********/
/* utility */
/***********/
	char		*str, command[MAXL], string[MAXL];
	int 		c, i, k, m, n;

/*********/
/* flags */
/*********/
	int		debug = 0;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
	t4file[1][0] = '\0';
/************************/
/* process command line */
/************************/
	for (n = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while (c = *str++) switch (c) {
				case 'd': debug++;		break;
			}
		} else if (n < 1) {
			strcpy (t4file[n++], argv[i]);
		}
	}
	if (n < 1) {
		printf ("Usage:	%s <t4file>\n", program);
		printf ("e.g.,	%s vm11b_mpr1_to_711-2B_t4\n", program);
		exit (1);
	}

	for (k = 0; k < 16; k++)    t4[0][k] = 0.;
	for (k = 0; k < 16; k += 5) t4[0][k] = 1.;

	if (!(fp = fopen (t4file[0], "w"))) errw (program, t4file[0]);
	printf ("Writing: %s\n", t4file[0]);
	write_command_line (fp, argc, argv);
	t4list (fp, t4[0]);
	if (fclose (fp)) errw (program, t4file[0]);

	exit (0);
}
