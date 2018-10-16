c$Header: /data/petsun4/data1/src_solaris/parzen_4dfp/RCS/parzen.f,v 1.1 2007/07/04 02:10:26 avi Exp $
c$Log: parzen.f,v $
c Revision 1.1  2007/07/04  02:10:26  avi
c Initial revision
c
c	2003.02.10 - Designing Parzen 'filters' for 3D images.
c	Definition of arguments: nx,ny,nz - #voxels in the x,y,z directions.
c	voxdim, fwhm are in mm !.

      subroutine parzen(image,nx,ny,nz,voxdim,fwhm)
      real*4 image(nx,ny,nz),voxdim(3)
      real*4 M,pi,sx,sy,sz
      pointer (pa,a),(pb,b)
      real*4 a((nx/2+1)*ny*nz),b((nx/2+1)*ny*nz)
      character*256 rcsid
     &/'$Id: parzen.f,v 1.1 2007/07/04 02:10:26 avi Exp $'/
      data pi/3.141592653589793/

      l=lnblnk(rcsid)
      write(*,"(a)")rcsid(1:l)
      pa=malloc(4*(nx/2+1)*ny*nz)
      pb=malloc(4*(nx/2+1)*ny*nz)

      M=fwhm/0.722351745       ! correct

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
        fz=float(iz-1)/(float(nz)*voxdim(3))
      else
        fz=float(nz-iz+1)/(float(nz)*voxdim(3))
      endif
      if(fz.ne.0.0)then
         sz=(sinc(pi*fz*M/2))**4
      else
         sz=1.0
      endif
      do 31 iy=1,ny
      if(iy.le.ny/2+1)then
        fy=float(iy-1)/(float(ny)*voxdim(2))
      else
        fy=float(ny-iy+1)/(float(ny)*voxdim(2))
      endif
      if(fy.ne.0.0)then
         sy=(sinc(pi*fy*M/2))**4
      else
         sy=1.0
      endif
      do 31 ix=1,n1
        fx=float(ix-1)/(float(nx)*voxdim(1))
      if(fx.ne.0.0)then
         sx=(sinc(pi*fx*M/2))**4
      else
         sx=1.0
      endif
      factor=sx*sy*sz
      a(i)=factor*a(i)
      b(i)=factor*b(i)
   31 i=i+1
      fx=0.5/voxdim(1)
      sx=(sinc(pi*fx*M/2))**4
      do 32 iz=1,nz
      if(iz.le.nz/2+1)then
        fz=float(iz-1)/(float(nz)*voxdim(3))
      else
        fz=float(nz-iz+1)/(float(nz)*voxdim(3))
      endif
      if(fz.ne.0.0)then
         sz=(sinc(pi*fz*M/2))**4
      else
         sz=1.0
      endif
      do 32 iy=1,ny
      if(iy.le.ny/2+1)then
        fy=float(iy-1)/(float(ny)*voxdim(2))
      else
        fy=float(ny-iy+1)/(float(ny)*voxdim(2))
      endif
      if(fy.ne.0.0)then
         sy=(sinc(pi*fy*M/2))**4
      else
         sy=1.0
      endif
      factor=sx*sy*sz
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
        
      real*4 function sinc(x)
      real*4 x
      if(x.ne.0.0) then
         sinc=sin(x)/x
      else 
         sinc=1.0
      endif
C      write(*,"(2e16.6)") x, sinc
      return
      end  
