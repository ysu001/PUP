/*$Header: /home/usr/shimonyj/JSSutil/RCS/lin_algebra.c,v 1.4 2012/02/22 23:04:26 avi Exp $*/
/*$Log: lin_algebra.c,v $
 * Revision 1.4  2012/02/22  23:04:26  avi
 * trap errors in gaussj_flg();
 *
 * Revision 1.3  2008/11/10  23:43:24  shimonyj
 * added gaussj_flg
 *
 * Revision 1.2  2007/08/28  03:52:58  avi
 *  remove #include <omp.h>
 *
 * Revision 1.1  2007/08/28  03:33:32  avi
 * Initial revision
 **/

#include <stdio.h>
#include <math.h>
#include <JSSutil.h>

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
#define SIGN(a,b) ((b) > 0 ? fabs(a) : -fabs(a))

/* jacobi method to find eigenvalues of a real symmetric matrix NR pg. 467 */
void jacobi(float **a,int n,float d[],float **v,int *nrot)
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	
	b = vector(1,n);
	z = vector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	*nrot = 0;
	for (i=1; i<=50; i++) {
		sm = 0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				sm += fabs(a[ip][iq]);
			}
		}
		if (sm == 0.0) {
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}
		if (i<4) tresh = 0.2*sm/(n*n);
		else tresh = 0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g = 100.0*fabs(a[ip][iq]);
				if (i>4 && (float) (fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h = d[iq] - d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t = (a[ip][iq])/h;
					else {
						theta = 0.5*h/(a[ip][iq]);
						t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau = s/(1.0+c);
					h = t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}	
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}	
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	free_vector(z,1,n);
	free_vector(b,1,n);
	printf("Error(jacobi): too many iterations\n");
}

/* sort eigenvalues and vectors NR pg. 468 */
void eigsrt(float d[], float **v, int n)
{
	int k,j,i;
	float p;
	
	for (i=1;i<n;i++) {
		p = d[k=i];
		for (j=i+1;j<=n;j++) {
			if (d[j] >= p) p=d[k=j];
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for (j=1;j<=n;j++) {
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}


/* SVD decomposition NR p. 67 */
void svdcmp(float **a, int m, int n, float w[], float **v)
{
	float pythag(float a, float b);
	int flag,i,its,j,jj,k,l,nm;
	float anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1 = vector(1,n);
	g = scale = anorm = 0.0;
	for (i=1; i<=n; i++) {
		l=i+1;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				a[i][i] = f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f = s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i] = scale*g;
		g = s = scale = 0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				a[i][l] = f-g;
				for (k=l;k<=n;k++) rv1[k] = a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g = rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l = i+1;
		g = w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g = 1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f = (s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		}
		else for (j=i;j<=m;j++) a[j][i] = 0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1; its<=30; its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm = l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((float)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0; s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;
					g = w[i];
					h = pythag(f,g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y*c + z*s;
						a[j][i] = z*c - y*s;
					}
				}
			}
			z=w[k];
			if (l==k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) printf("Error(svdcmp): no convergence in 30 iterations\n");
			x = w[l];
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c = s = 1.0;
			for (j=l;j<=nm;j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y*s;
				y *= c;
				for (jj=1; jj<=n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x*c + z*s;
					v[jj][i] = z*c - x*s;
				}
				z = pythag(f,h);
				w[j] = z;
				if (z) {
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g + s*y;
				x = c*y - s*g;
				for (jj=1;jj<=m;jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y*c + z*s;
					a[jj][i] = z*c - y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_vector(rv1,1,n);
}

float pythag(float a, float b)
{
	float absa,absb;
	
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

/* invert a matrix using gauss jordan method NR p. 36 */
void gaussj(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i, icol,irow,j,k,l,ll;
	float big,dum,pivinv,temp;
	float swap;
	
	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for(i=1;i<=n;i++) {
		big=0.0;
		for(j=1;j<=n;j++) {
			if(ipiv[j] !=1) {
				for(k=1;k<=n;k++) {
					if(ipiv[k] == 0) {
						if(fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} 
					else if (ipiv[k] > 1) printf("\nError(gaussj): singular matrix-1\n");
				}
			}
		}
		++(ipiv[icol]);
		if(irow != icol) {
			for(l=1;l<=n;l++) {swap=a[irow][l];a[irow][l]=a[icol][l];a[icol][l]=swap;}
			for(l=1;l<=m;l++) {swap=b[irow][l];b[irow][l]=b[icol][l];b[icol][l]=swap;}
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) printf("\nError(gaussj): singular matrix-2\n");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *=pivinv;
		for (l=1;l<=m;l++) b[icol][l] *=pivinv;
		for(ll=1;ll<=n;ll++) {
			if(ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for(l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for(l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
		}
	}
	for (l=n;l>=1;l--) {
		if(indxr[l] != indxc[l]) {
			for(k=1;k<=n;k++) {
				swap=a[k][indxr[l]];a[k][indxr[l]]=a[k][indxc[l]];a[k][indxc[l]]=swap;
			}
		}
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
	
	

/* invert a matrix using gauss jordan method NR p. 36 */
/* return flag=1 if singular matrix, otherwise return 0 */
int gaussj_flg(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i, icol,irow,j,k,l,ll,status;
	float big,dum,pivinv,temp;
	float swap;
	
	status = 0;
	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for(i=1;i<=n;i++) {
		big=0.0;
		for(j=1;j<=n;j++) {
			if(ipiv[j] !=1) {
				for(k=1;k<=n;k++) {
					if(ipiv[k] == 0) {
						if(fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} 
					else if (ipiv[k] > 1) {status = 1; goto FREE;}
				}
			}
		}
		++(ipiv[icol]);
		if(irow != icol) {
			for(l=1;l<=n;l++) {swap=a[irow][l];a[irow][l]=a[icol][l];a[icol][l]=swap;}
			for(l=1;l<=m;l++) {swap=b[irow][l];b[irow][l]=b[icol][l];b[icol][l]=swap;}
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) return(1);
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *=pivinv;
		for (l=1;l<=m;l++) b[icol][l] *=pivinv;
		for(ll=1;ll<=n;ll++) {
			if(ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for(l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for(l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
		}
	}
	for (l=n;l>=1;l--) {
		if(indxr[l] != indxc[l]) {
			for(k=1;k<=n;k++) {
				swap=a[k][indxr[l]];a[k][indxr[l]]=a[k][indxc[l]];a[k][indxc[l]]=swap;
			}
		}
	}
FREE:	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
	return(status);
}

void covsrt(float **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	float swap;
	
	for(i=mfit+1;i<=ma;i++)
		for(j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for(j=ma;j>=1;j--) {
		if(ia[j]) {
			for(i=1;i<=ma;i++) {swap=covar[i][k];covar[i][k]=covar[i][j];covar[i][j]=swap;}
			for(i=1;i<=ma;i++) {swap=covar[k][i];covar[k][i]=covar[j][i];covar[j][i]=swap;}
			k--;
		}
	}
}
