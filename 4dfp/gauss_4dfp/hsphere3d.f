cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, Washington University,
c Mallinckrodt Institute of Radiology. All Rights Reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of Washington University. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Header: /data/petsun4/data1/src_solaris/gauss_4dfp/RCS/hsphere3d.f,v 1.2 2007/04/17 05:25:58 avi Exp $
c $Log: hsphere3d.f,v $
c Revision 1.2  2007/04/17  05:25:58  avi
c gcc compliant
c
c Revision 1.1  2001/07/11  04:22:13  avi
c Initial revision
c
      subroutine hsphere3d(image,nx,ny,nz,mmppix,radius)
c     frequencey domain filter for hard sphere kernel
c     modified gauss3d.f
c     Author: Avi Snyder 09/Jul/01
c     analyztic expresion for factor derived with help of D. Yablonskiy
c
c     inputs: real*4 image(nx,ny,nz)	! image to be filtered
c             real*4 mmppix(3)		! x, y, z voxel dimensions in mm
c             real*4 radius		! hard shpere radius in mm
c
c     variables modified: image         ! input image overwritten
c
c     subroutines called:  FFT, REALT in FORTRAN source fftsun.f
c                          FFT algorithm by Richard C. Singleton via FFTPE.FTN
c
c     restrictions: mixed radix FFT accepts any nx, ny, nz but 
c                   nx must be divisible by 2 
c                   memory allocation requirement 8*(nx/2+1)*ny*nz bytes 

      real*4 image(nx,ny,nz)
      real*4 mmppix(3)
      pointer (pa,a),(pb,b)
      real*4 a((nx/2+1)*ny*nz),b((nx/2+1)*ny*nz)

      if(mod(nx,2).ne.0)then
        write(*,"('hsphere3d: nx=',i4,' not divisible by 2')")nx
        call exit(-1)
      endif
      twopi=8.0*atan(1.)
      pa=malloc(4*(nx/2+1)*ny*nz)
      pb=malloc(4*(nx/2+1)*ny*nz)

      i=1
      do 21 iz=1,nz
      do 21 iy=1,ny
      do 21 ix=1,nx,2
      a(i)=image(ix,iy,iz)
      b(i)=image(ix+1,iy,iz)
   21 i=i+1

      n1=nx/2
      n=n1*ny*nz
      call FFT(a,b,ny*nz,n1,1,-1)
      call REALT(a,b,ny*nz,n1,1,-1)
      call FFT(a(n+1),b(n+1),nz,ny,1,-1)
      call FFT(a(n+1),b(n+1),1,nz,ny,-1)
      call FFT(a,b,nz,ny,n1,-1)
      call FFT(a,b,1,nz,n1*ny,-1)

      i=1
      do 31 iz=1,nz
      if(iz.le.nz/2+1)then
        fz=float(iz-1)/(float(nz)*mmppix(3))
      else
        fz=float(nz-iz+1)/(float(nz)*mmppix(3))
      endif
      do 31 iy=1,ny
      if(iy.le.ny/2+1)then
        fy=float(iy-1)/(float(ny)*mmppix(2))
      else
        fy=float(ny-iy+1)/(float(ny)*mmppix(2))
      endif
      do 31 ix=1,n1
      fx=float(ix-1)/(float(nx)*mmppix(1))
      f=sqrt(fx**2+fy**2+fz**2)
      if(f.eq.0.0)then
        factor=1.0
      else
        theta=twopi*f*radius
        factor=3.0*(sin(theta)-theta*cos(theta))/theta**3
      endif
      a(i)=factor*a(i)
      b(i)=factor*b(i)
   31 i=i+1
      fx=0.5/mmppix(1)
      do 32 iz=1,nz
      if(iz.le.nz/2+1)then
        fz=float(iz-1)/(float(nz)*mmppix(3))
      else
        fz=float(nz-iz+1)/(float(nz)*mmppix(3))
      endif
      do 32 iy=1,ny
      if(iy.le.ny/2+1)then
        fy=float(iy-1)/(float(ny)*mmppix(2))
      else
        fy=float(ny-iy+1)/(float(ny)*mmppix(2))
      endif
      f=sqrt(fx**2+fy**2+fz**2)
      theta=twopi*f*radius
      factor=3.0*(sin(theta)-theta*cos(theta))/theta**3
      a(i)=factor*a(i)
      b(i)=factor*b(i)
   32 i=i+1

      call FFT(a,b,1,nz,n1*ny,+1)
      call FFT(a,b,nz,ny,n1,+1)
      call FFT(a(n+1),b(n+1),1,nz,ny,+1)
      call FFT(a(n+1),b(n+1),nz,ny,1,+1)
      call REALT(a,b,ny*nz,n1,1,+1)
      call FFT(a,b,ny*nz,n1,1,+1)

      i=1
      do 29 iz=1,nz
      do 29 iy=1,ny
      do 29 ix=1,nx,2
      image(ix,iy,iz)=a(i)
      image(ix+1,iy,iz)=b(i)
   29 i=i+1

      call free(pa)
      call free(pb)
      return
      end
