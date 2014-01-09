#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <nrutil.h>


#ifndef SMALL
#define SMALL 1.0e-2
#endif

float loganREF(float *REF, float *ROI, float k2, float *t,  float *frd, int nframes, int stframe, int endframe, float *intercept, float *slope, float *Rsq)
{
float BP, r, *x, *y, *intgREF, *intgROI, *tmpx, *tmpy, mx, Sxx, Sx, k, Syy, Sy, my, xt, yt, Sxy,delta, sxx, syy, sxy;
int i;

intgREF=vector(0, nframes-1);
intgROI=vector(0, nframes-1);
tmpx=vector(0, nframes-1);
tmpy=vector(0, nframes-1);
x=vector(0, nframes-1);
y=vector(0, nframes-1);

tmpx[0]=REF[0]*t[0]/2./60.;
tmpy[0]=ROI[0]*t[0]/2./60.;

intgREF[0]=tmpx[0];
intgROI[0]=tmpy[0];
x[0]=(intgREF[0]+REF[0]/k2)/(ROI[0]+SMALL);
y[0]=intgROI[0]/(ROI[0]+SMALL);



for (i=1; i<nframes; i++)
	{
	intgREF[i]=intgREF[i-1]+(REF[i-1]+REF[i])*(t[i]-t[i-1])/2./60.;
	intgROI[i]=intgROI[i-1]+(ROI[i-1]+ROI[i])*(t[i]-t[i-1])/2./60.;
	x[i]=(intgREF[i]+REF[i]/k2)/(ROI[i]+SMALL);
	y[i]=intgROI[i]/(ROI[i]+SMALL);		
	}

k=(float)(endframe-stframe+1);
mx=0.;
Sxx=0.;
my=0.;
Syy=0.;
Sxy=0.;

for (i=stframe-1; i<endframe; i++)
	{
	mx += x[i];
	Sxx += x[i]*x[i];
	
	my += y[i];
	Syy += y[i]*y[i];
	
	Sxy += x[i]*y[i]; 
	}	

Sx = mx;
mx /= k;
sxx = Sxx - k*mx*mx;
delta = k*Sxx -Sx*Sx;

Sy = my;
my /= k;
syy = Syy - k*my*my;

sxy=0;
for (i=stframe-1; i<endframe; i++)
	{
	xt = x[i] - mx;
	yt = y[i] - my;
	sxy += xt*yt;
	}
if (sxx*syy>0) {
	r=sxy/sqrtf(sxx*syy);
	} else {
	r=0;
	}
*intercept=(Sxx*Sy-Sx*Sxy)/delta;
*slope=(k*Sxy-Sx*Sy)/delta;
*Rsq=r*r;

BP=*slope-1.;
return(BP);		
}
