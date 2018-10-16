/*$Header$*/
/*$Log$*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <librms.h>

#define MAXL	256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

int main (int argc, char **argv) {
        int  		imgdim[3] = {64, 80, 32};
        float	 	voxdim[3] = {1., 2., 3.};	
	float		f0 = 0.02, crit = 0.3;
	float		*imgt;
	int		*mask;
	int		vdim;

/***********/
/* utility */
/***********/
	char 		program[MAXL];

	setprog (program, argv);

	vdim = imgdim[0]*imgdim[1]*imgdim[2];
	if (!(imgt = (float *) calloc (vdim, sizeof (float)))) errm (program);
	if (!(mask = (int *)   calloc (vdim, sizeof (int))))   errm (program);

	printf ("before img2lmask\n");
	img2lmask (imgdim+0, imgdim+1, imgdim+2, imgt, mask, voxdim, &f0, &crit);
	printf ("after  img2lmask\n");

	free (imgt);
	free (mask);
	exit (0);
}
