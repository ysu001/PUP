C$Id: hipass3d.f,v 1.1 1996/04/19 17:10:11 ty7777 Exp $
C$Log: hipass3d.f,v $
c Revision 1.1  1996/04/19  17:10:11  ty7777
c Initial revision
c
c
      subroutine hipass3d (image, nx, ny, nz, cmppix, f0)

c     Author: Avi Snyder 19-Dec-95	! modified gauss3d.f
c
c     inputs: real*4 image(nx,ny,nz)	! image to be filtered
c             real*4 cmppix(3)		! x, y, z voxel dimensions in cm
c             real*4 f0			! half amplitude frequency in cycles/cm
c
c     variables modified: image         ! input image overwritten
c
c     subroutines called:  FFT, REALT in FORTRAN source fftsun.f
c                          FFT algorithm by Richard C. Singleton via FFTPE.FTN
c
c     restrictions: mixed radix FFT accepts any nx, ny, nz but 
c                   nx must be divisible by 2 
c                   memory allocation requirement 8*(nx/2+1)*ny*nz bytes 

      character*256 rcsid/'$Header: /data/petsun4/src_solaris/librms/RCS/hipass3d.f,v 1.1 1996/04/19 17:10:11 ty7777 Exp $'/
      real*4 image(nx,ny,nz)
      real*4 cmppix(3)
      pointer (pa,a),(pb,b)
      real*4 a((nx/2+1)*ny*nz),b((nx/2+1)*ny*nz)

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
        fz=float(iz-1)/(float(nz)*cmppix(3))
      else
        fz=float(nz-iz+1)/(float(nz)*cmppix(3))
      endif
      do 31 iy=1,ny
      if(iy.le.ny/2+1)then
        fy=float(iy-1)/(float(ny)*cmppix(2))
      else
        fy=float(ny-iy+1)/(float(ny)*cmppix(2))
      endif
      do 31 ix=1,n1
      fx=float(ix-1)/(float(nx)*cmppix(1))
      f=sqrt(fx**2+fy**2+fz**2)
      factor=1.-exp(-f/f0)
      a(i)=factor*a(i)
      b(i)=factor*b(i)
   31 i=i+1
      fx=0.5/cmppix(1)
      do 32 iz=1,nz
      if(iz.le.nz/2+1)then
        fz=float(iz-1)/(float(nz)*cmppix(3))
      else
        fz=float(nz-iz+1)/(float(nz)*cmppix(3))
      endif
      do 32 iy=1,ny
      if(iy.le.ny/2+1)then
        fy=float(iy-1)/(float(ny)*cmppix(2))
      else
        fy=float(ny-iy+1)/(float(ny)*cmppix(2))
      endif
      f=sqrt(fx**2+fy**2+fz**2)
      factor=1.-exp(-f/f0)
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
