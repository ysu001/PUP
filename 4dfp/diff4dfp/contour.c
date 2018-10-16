#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <knee.h>

/* 
** name - con_map
** use  - find all the contours at given grey level and store
** input - ncont: the number of contours up to now
**         z: the level of contour desired
**         src: pointer to source image
** output- returns the number of contours found for this z
** effect- fill the cont array 
** note - assumes the image is 256x256
*/

int
con_map(int ncont, int z, short int *src)
{
    int i,j,nall,ncon,x,y,hv,xbeg,ybeg,hvbeg;
    int xold,yold,area,perim;
    int *shor, *sver;
    shor = (int *) malloc(DZ*sizeof(int));
    sver = (int *) malloc(DZ*sizeof(int));

/* set the borders to zero */
    for (i=0; i<DY; ++i) {
        *(src+i) = 0;
	*(src+i*DY) = 0;
	*(src+255*DY+i) = 0;
	*(src+i*DY+255) = 0;
    }

/* set the horizontal and vertical differences and starts */
    for (i=0; i<255; ++i) {
	for (j=0; j<255; ++j) {
	    if (*(src+i*DY+j+1)>=z && *(src+i*DY+j)<z)
	        *(shor+i*DY+j) = 1;
            else if (*(src+i*DY+j+1)<z && *(src+i*DY+j)>=z)
	        *(shor+i*DY+j) = 2;
	    else *(shor+i*DY+j) = 0;
	    if (*(src+(i+1)*DY+j)>=z && *(src+i*DY+j)<z)
	        *(sver+i*DY+j) = 3;
	    else if (*(src+(i+1)*DY+j)<z && *(src+i*DY+j)>=z)
	        *(sver+i*DY+j) = 4;
	    else *(sver+i*DY+j) = 0;
        }
    }

/* nall - count all contours, ncon - count only large ones */
nall=0; ncon=0;
/* loop till find all contours in image around z */
while(1) {
/* loop till locate a start arrow then break out */
    xbeg = ybeg = 0;
    for (i=0; i<255; ++i) {
	for (j=0; j<255; ++j) {
            if (*(shor+i*DY+j) != 0) {
		xbeg=j; ybeg=i; hvbeg=0; break;
            }
            if (*(sver+i*DY+j) != 0) {
		xbeg=j; ybeg=i; hvbeg=1; break;
            }
        }
	if (xbeg!=0 || ybeg!=0) break;
    }
/* no more contours, so return */
    if (xbeg==0 && ybeg==0) {
	return(ncon);
    }

/* trace the outline of the found curve, erase shor,sver */
    ++nall;
    area = 0; perim = 0;
    xold = 0; yold = 0;
    x = xbeg; y = ybeg; hv = hvbeg;
/* continue looping until return to the beg point */
    do {
	++perim;
	if (perim>5000) {
	    printf("error - too many iterations\n");
	    break;
        } 
	area += yold*x - xold*y;
	xold = x; yold = y;
        if (hv==0) { /* horizontal */
            *(shor+y*DY+x) = 0;
	    if (*(src+y*DY+x+1) > *(src+y*DY+x)) { /* DOWN */
		if (*(src+y*DY+x+1) >= z &&
		    *(src+(y+1)*DY+x+1) < z) { /* left */
		    x +=1; hv=1;
                }
		else if (*(src+y*DY+x) < z &&
			 *(src+(y+1)*DY+x) >= z) { /* right */
		    hv=1;
                }
		else { /* forward */ 
		    y +=1;
                }
            }
            else { /* UP */
		if (*(src+(y-1)*DY+x) < z &&
		    *(src+y*DY+x) >= z) { /* left */
		    y -=1; hv=1;
                }
		else if (*(src+(y-1)*DY+x+1) >= z &&
			 *(src+y*DY+x+1) < z) { /* right */
		    y -=1; x +=1; hv=1;
                }
		else { /* forward */
		    y -=1;
                }
            }
        }
        else { /* vertical */
            *(sver+y*DY+x) = 0;
	    if (*(src+(y+1)*DY+x) > *(src+y*DY+x)) { /* LEFT */
		if (*(src+(y+1)*DY+x-1) < z &&
		    *(src+(y+1)*DY+x) >= z) { /* left */
		    y +=1; x -=1; hv=0;
                }
		else if (*(src+y*DY+x-1) >= z &&
			 *(src+y*DY+x) < z) { /* right */
		    x -=1; hv=0;
                }
		else { /* forward */ 
		    x -=1;
                }
            }
            else { /* RIGHT */
		if (*(src+y*DY+x) >= z &&
		    *(src+y*DY+x+1) < z) { /* left */
		    hv=0;
                }
		else if (*(src+(y+1)*DY+x) < z &&
			 *(src+(y+1)*DY+x+1) >= z) { /* right */
		    y +=1; hv=0;
                }
		else { /* forward */
		    x +=1;
                }
            }
        }
    } while ( x!=xbeg || y!=ybeg || hv!=hvbeg );
/* fix the area */
    area += yold*x - xold*y;
    if (area>0) ++area; else --area;
    area /= 2;
/* check if contour big enough to be saved */
    if (abs(area) > CON_AREA) {
	++ncon;
	if (ncont+ncon >= CON_NMAX) {
	    printf("error - too many contours\n");
	    return(ncon-1);
        }
	cont[ncont].xbeg = xbeg;
	cont[ncont].ybeg = ybeg;
	cont[ncont].hvbeg = hvbeg;
	cont[ncont].z = z;
	cont[ncont].area = area;
	cont[ncont].perim = perim;
/* these will be filled later */
	cont[ncont].parent = 0;
	cont[ncont].depth = 0;
	cont[ncont].nson = 0;
/*
        printf("path %d: beg x=%d y=%d hv=%d area=%d perim=%d\n",
            ncon,xbeg,ybeg,hvbeg,area,perim);
*/
    }
} /* close for the while(1) line */
}

/*
** name - con_parents
** use  - find parent for each contour
** input - ncont: number of contours
**         src: pointer to source image
** output- none
** effect- fill the parent section of contour data
** note - assume the out side contour z=1 is in location 0
*/
void
con_parents(int ncont, short int *src)
{
    int i,j,x,y,hv,zs,zd,hflag;

/* loop over all contours except the first */
    for (i=1; i<ncont; ++i) {
/* determine the possible z's of the parent */
        if (cont[i].area > 0) zd = cont[i].z - CON_STEP;
        else zd = cont[i].z + CON_STEP;
        zs = cont[i].z;
/* find the next location a znew contour would pass */
        hflag = 0;
        hv = cont[i].hvbeg;
        if (hv == 0) {
	    y = cont[i].ybeg;
	    for (x=cont[i].xbeg; x>=0; --x) {
	        if ((*(src+y*DY+x+1)>=zd && *(src+y*DY+x)<zd) ||
                    (*(src+y*DY+x+1)<zd && *(src+y*DY+x)>=zd)) {
                    if (hflag=parent_flag(ncont,i,1,zd,src,x,y,hv)) break;
                }
	        if (x<cont[i].xbeg) {
	        if ((*(src+y*DY+x+1)>=zs && *(src+y*DY+x)<zs) ||
                    (*(src+y*DY+x+1)<zs && *(src+y*DY+x)>=zs)) {
                    if (hflag=parent_flag(ncont,i,2,zs,src,x,y,hv)) break;
                }
                }
            }
        }
	else {
	    x = cont[i].xbeg;
	    for (y=cont[i].ybeg; y>=0; --y) {
	        if ((*(src+(y+1)*DY+x)>=zd && *(src+y*DY+x)<zd) ||
                    (*(src+(y+1)*DY+x)<zd && *(src+y*DY+x)>=zd)) {
                    if (hflag=parent_flag(ncont,i,1,zd,src,x,y,hv)) break;
                }
	        if (y<cont[i].ybeg) {
	        if ((*(src+(y+1)*DY+x)>=zs && *(src+y*DY+x)<zs) ||
                    (*(src+(y+1)*DY+x)<zs && *(src+y*DY+x)>=zs)) {
                    if (hflag=parent_flag(ncont,i,2,zs,src,x,y,hv)) break;
                }
                }
            }
        }
        if (hflag==0) {
            printf("No parent for number %d\n",i);
            continue;
        }
    }
}

/* 
** name - parent_flag
** use  - store found parent and return flag
** inputs -
**     ncont: the number of contours
**     i: the contour who needs a parent
**     hflag: type of possible parent found
**     z: the z for possible parent
**     src: pointer to the original image
**     x,y,hv: the point through which the parent should pass
** output - flag=1 if parent found, else 0
** effect - fill parent number in contour data
** note - none
*/
int
parent_flag(int ncont, int i, int hflag, int z, short int *src, int x, int y, int hv)
{
    int j;

    for (j=0; j<ncont; ++j) {
        if (cont[j].z != z) continue;
        if (hflag==1) {
            if (cont[j].area*cont[i].area < 0) continue;
	    if (trace_flag(j,src,x,y,hv)!=0) {
		cont[i].parent = j;
		return(1);
            }
        }
        else if (hflag==2) {
            if (cont[j].area*cont[i].area > 0) continue;
	    if (trace_flag(j,src,x,y,hv)!=0) {
		cont[i].parent = j;
		return(1);
            }
        }
    }
    return(0);
}

/*
** name - con_trace 
** use  - trace a contour given the start point
** input - icont: contour to trace
**         src: pointer to source image
**         grey: grey value to draw trace
**         dst: pointer to destination image
** output - none
** effect - draw trace in dst with grey value
** note - none
*/
void
con_trace(int icont, short int *src, int grey, short int *dst)
{
    int x,y,z,hv;

    x = cont[icont].xbeg;
    y = cont[icont].ybeg;
    hv = cont[icont].hvbeg;
    z = cont[icont].z;
    do {
	*(dst+y*DY+x) = grey;
        if (hv==0) { /* horizontal */
	    if (*(src+y*DY+x+1) > *(src+y*DY+x)) { /* DOWN */
		if (*(src+y*DY+x+1) >= z &&
		    *(src+(y+1)*DY+x+1) < z) { /* left */
		    x +=1; hv=1;
                }
		else if (*(src+y*DY+x) < z &&
			 *(src+(y+1)*DY+x) >= z) { /* right */
		    hv=1;
                }
		else { /* forward */ 
		    y +=1;
                }
            }
            else { /* UP */
		if (*(src+(y-1)*DY+x) < z &&
		    *(src+y*DY+x) >= z) { /* left */
		    y -=1; hv=1;
                }
		else if (*(src+(y-1)*DY+x+1) >= z &&
			 *(src+y*DY+x+1) < z) { /* right */
		    y -=1; x +=1; hv=1;
                }
		else { /* forward */
		    y -=1;
                }
            }
        }
        else { /* vertical */
	    if (*(src+(y+1)*DY+x) > *(src+y*DY+x)) { /* LEFT */
		if (*(src+(y+1)*DY+x-1) < z &&
		    *(src+(y+1)*DY+x) >= z) { /* left */
		    y +=1; x -=1; hv=0;
                }
		else if (*(src+y*DY+x-1) >= z &&
			 *(src+y*DY+x) < z) { /* right */
		    x -=1; hv=0;
                }
		else { /* forward */ 
		    x -=1;
                }
            }
            else { /* RIGHT */
		if (*(src+y*DY+x) >= z &&
		    *(src+y*DY+x+1) < z) { /* left */
		    hv=0;
                }
		else if (*(src+(y+1)*DY+x) < z &&
			 *(src+(y+1)*DY+x+1) >= z) { /* right */
		    y +=1; hv=0;
                }
		else { /* forward */
		    x +=1;
                }
            }
        }
    } while ( x!=cont[icont].xbeg || y!=cont[icont].ybeg ||
              hv!=cont[icont].hvbeg );
}


/*
** name - trace_flag 
** use  - check if a point is in a contour
** input - icont: contour number to trace
**         src: pointer to source image
**         xif,yif,hvif: point to check if in contour
** output - return 1 if yes, 0 if no
** effect - none
** note - none
*/
int
trace_flag(int icont, short int *src, int xif, int yif, int hvif)
{
    int x,y,z,hv;

    x = cont[icont].xbeg;
    y = cont[icont].ybeg;
    hv = cont[icont].hvbeg;
    z = cont[icont].z;
    do {
	if (x==xif && y==yif && hv==hvif) return(1);
        if (hv==0) { /* horizontal */
	    if (*(src+y*DY+x+1) > *(src+y*DY+x)) { /* DOWN */
		if (*(src+y*DY+x+1) >= z &&
		    *(src+(y+1)*DY+x+1) < z) { /* left */
		    x +=1; hv=1;
                }
		else if (*(src+y*DY+x) < z &&
			 *(src+(y+1)*DY+x) >= z) { /* right */
		    hv=1;
                }
		else { /* forward */ 
		    y +=1;
                }
            }
            else { /* UP */
		if (*(src+(y-1)*DY+x) < z &&
		    *(src+y*DY+x) >= z) { /* left */
		    y -=1; hv=1;
                }
		else if (*(src+(y-1)*DY+x+1) >= z &&
			 *(src+y*DY+x+1) < z) { /* right */
		    y -=1; x +=1; hv=1;
                }
		else { /* forward */
		    y -=1;
                }
            }
        }
        else { /* vertical */
	    if (*(src+(y+1)*DY+x) > *(src+y*DY+x)) { /* LEFT */
		if (*(src+(y+1)*DY+x-1) < z &&
		    *(src+(y+1)*DY+x) >= z) { /* left */
		    y +=1; x -=1; hv=0;
                }
		else if (*(src+y*DY+x-1) >= z &&
			 *(src+y*DY+x) < z) { /* right */
		    x -=1; hv=0;
                }
		else { /* forward */ 
		    x -=1;
                }
            }
            else { /* RIGHT */
		if (*(src+y*DY+x) >= z &&
		    *(src+y*DY+x+1) < z) { /* left */
		    hv=0;
                }
		else if (*(src+(y+1)*DY+x) < z &&
			 *(src+(y+1)*DY+x+1) >= z) { /* right */
		    y +=1; hv=0;
                }
		else { /* forward */
		    x +=1;
                }
            }
        }
    } while ( x!=cont[icont].xbeg || y!=cont[icont].ybeg ||
              hv!=cont[icont].hvbeg );
    return(0);
}

/*
** name - con_depth
** use  - fill in the depth and nson of image tree
** input - ncont: number of contours
**         cont: pointer to contour array
** output- none
** effect- fill the nson and depth data for the full tree
** note - none
*/
void
con_depth(int ncont, CONT cont[])
{
    int i,j,dep;

/* loop over all contours */
    for (i=1; i<ncont; ++i) {
        dep = 1;
/* j is i's parent */
        j = cont[i].parent;
/* increment j's son count */
        ++cont[j].nson;
/* count depth */
        while (j!=0) {
            j = cont[j].parent;
            ++dep;
        }
        cont[i].depth = dep;
    }
}

/*
** con_fill - fill a contour given the start point
** mainly for drawing purposes in dst
** name - con_fill
** use  - fill a contour in dst for display
** input - icont: contour to trace
**         src: pointer to source image
**         grey: grey value to fill
**         dst: pointer to destination image
** output - none
** effect - draw fill in dst with grey value
** note - none
*/
void
con_fill(int icont, short int *src, int grey, short int *dst)
{
    int x,y,z,hv,i,j,jscan;
    short int *hdir;
    hdir = (short int *) malloc(DZ*sizeof(short int));

    for (i=0; i<DZ; ++i) *(hdir+i)=0;
    x = cont[icont].xbeg;
    y = cont[icont].ybeg;
    hv = cont[icont].hvbeg;
    z = cont[icont].z;
    do {
        if (hv==0) { /* horizontal */
	    if (*(src+y*DY+x+1) > *(src+y*DY+x)) { /* DOWN */
	        *(hdir+y*DY+x) = 1;
		if (*(src+y*DY+x+1) >= z &&
		    *(src+(y+1)*DY+x+1) < z) { /* left */
		    x +=1; hv=1;
                }
		else if (*(src+y*DY+x) < z &&
			 *(src+(y+1)*DY+x) >= z) { /* right */
		    hv=1;
                }
		else { /* forward */ 
		    y +=1;
                }
            }
            else { /* UP */
	        *(hdir+y*DY+x) = 2;
		if (*(src+(y-1)*DY+x) < z &&
		    *(src+y*DY+x) >= z) { /* left */
		    y -=1; hv=1;
                }
		else if (*(src+(y-1)*DY+x+1) >= z &&
			 *(src+y*DY+x+1) < z) { /* right */
		    y -=1; x +=1; hv=1;
                }
		else { /* forward */
		    y -=1;
                }
            }
        }
        else { /* vertical */
	    if (*(src+(y+1)*DY+x) > *(src+y*DY+x)) { /* LEFT */
		if (*(src+(y+1)*DY+x-1) < z &&
		    *(src+(y+1)*DY+x) >= z) { /* left */
		    y +=1; x -=1; hv=0;
                }
		else if (*(src+y*DY+x-1) >= z &&
			 *(src+y*DY+x) < z) { /* right */
		    x -=1; hv=0;
                }
		else { /* forward */ 
		    x -=1;
                }
            }
            else { /* RIGHT */
		if (*(src+y*DY+x) >= z &&
		    *(src+y*DY+x+1) < z) { /* left */
		    hv=0;
                }
		else if (*(src+(y+1)*DY+x) < z &&
			 *(src+(y+1)*DY+x+1) >= z) { /* right */
		    y +=1; hv=0;
                }
		else { /* forward */
		    x +=1;
                }
            }
        }
    } while ( x!=cont[icont].xbeg || y!=cont[icont].ybeg ||
              hv!=cont[icont].hvbeg );

/* have list of all ups (2) and downs (1) in hdir */
    if (cont[icont].area>0) {
/* scan for the down arrows to the next up arrow */
        for (i=0; i<DY; ++i) {
            for (j=0; j<DY; ++j) {
                if (*(hdir+i*DY+j)!=1) continue;
                jscan=j;
                do {
                    ++jscan;
                    *(dst+i*DY+jscan)=grey;
                } while(*(hdir+i*DY+jscan)!=2);
            }
        }
    }
    else {
        for (i=0; i<DY; ++i) {
            for (j=0; j<DY; ++j) {
                if (*(hdir+i*DY+j)!=2) continue;
                jscan=j;
                do {
                    ++jscan;
                    *(dst+i*DY+jscan)=grey;
                } while(*(hdir+i*DY+jscan)!=1);
            }
        }
    }
}


/*
** name - tree_prt
** use  - print the contour tree from given ancestor
** input - ind: index to ancestor in cont array
** output- print out of tree structure with indentations
** effect- none
** note - recursive routine
*/
void
tree_prt(int ncont, int ind)
{
    int i;

/* first print */
    for (i=0; i<cont[ind].depth; ++i) printf("  ");
    printf("d%d:#%d a %d z %d\n",cont[ind].depth,ind,
        cont[ind].area,cont[ind].z);

/* if leaf return else recurse */
    if (cont[ind].nson==0) return;
    else for (i=1; i<ncont; ++i) {
            if (cont[i].parent==ind) tree_prt(ncont, i);
    }
}

