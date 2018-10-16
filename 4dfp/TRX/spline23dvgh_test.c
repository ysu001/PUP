/*$Header: /data/petsun4/data1/src_solaris/TRX/RCS/spline23dvgh_test.c,v 1.3 2007/04/22 02:25:40 avi Exp $*/
/*$Log: spline23dvgh_test.c,v $
 * Revision 1.3  2007/04/22  02:25:40  avi
 * remove splint2dvgh_testf_()
 *
 * Revision 1.2  2007/04/01  00:03:22  avi
 * break out spline23dvgh.h
 * command line control of test
 *
 * Revision 1.1  2005/12/16  04:22:09  avi
 * Initial revision
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <spline23dvgh.h>

#define MAXL	256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[] = "$Id: spline23dvgh_test.c,v 1.3 2007/04/22 02:25:40 avi Exp $";
int main (int argc, char *argv[]) {
/***********/
/* utility */
/***********/
	int		k;
	char		program[MAXL];

	printf ("%s\n", rcsid);
	setprog (program, argv);
	if (argc < 2) {
		printf ("Usage:	%s <(int) k>\n", program);
		printf ("	 k	test\n");
		printf ("	 1	spline2dvgh_rcs\n");
		printf ("	 2	splint2dv_test\n");
		printf ("	 3	splint2dvgh_testh\n");
		printf ("	 4	splint2dvgh_testg\n");
		printf ("	 5	spline3dvgh_rcs\n");
		printf ("	 6	splint3dv_test\n");
		printf ("	 7	splint3dvl_test\n");
		printf ("	 8	splint3dvgh_testh\n");
		printf ("	 9	splint3dvgh_testg\n");
		exit (1);
	}
	k = atoi (argv[1]);
	switch (k) {
		case 1:		spline2dvgh_rcs_ ();	break;
		case 2:		splint2dv_test_	();	break;
		case 3:		splint2dvgh_testh_ ();	break;
		case 4:		splint2dvgh_testg_ ();	break;

		case 5:		spline3dvgh_rcs_ ();	break;
		case 6:		splint3dv_test_ ();	break;
		case 7:		splint3dvl_test_ ();	break;
		case 8:		splint3dvgh_testh_ ();	break;
		case 9:		splint3dvgh_testg_ ();	break;
	}
	exit (0);
}
