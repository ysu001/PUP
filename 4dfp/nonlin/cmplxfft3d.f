c$Header: /data/petsun4/data1/src_solaris/nonlin/RCS/cmplxfft3d.f,v 1.1 2004/10/14 23:42:36 rsachs Exp $
c$Log: cmplxfft3d.f,v $
c Revision 1.1  2004/10/14  23:42:36  rsachs
c Initial revision
c
      subroutine cmplxfft3d(imgs,imgf,nx,ny,nz,idir)
      complex*8 imgs(0:nx-1,0:ny-1,0:nz-1),imgf(0:nx-1,0:ny-1,0:nz-1)
      pointer (pa,a),(pb,b)
      real*4 a(0:nx*ny*nz-1),b(0:nx*ny*nz-1)
      character*256 rcsid/'$Id: cmplxfft3d.f,v 1.1 2004/10/14 23:42:36 rsachs Exp $'/

      l=lnblnk(rcsid)
c     write(*,"(a)")rcsid(1:l)
      pa=malloc(4*nx*ny*nz)
      pb=malloc(4*nx*ny*nz)
      if(pa.eq.0.or.pb.eq.0)stop 'cmplxfft3d malloc error'

      i=0
      do 21 iz=0,nz-1
      do 21 iy=0,ny-1
      do 21 ix=0,nx-1
      if(idir.lt.0) then
        a(i)= real(imgs(ix,iy,iz))
        b(i)=aimag(imgs(ix,iy,iz))
      else
        a(i)= real(imgf(ix,iy,iz))
        b(i)=aimag(imgf(ix,iy,iz))
      endif
   21 i=i+1

      call FFT(a,b,ny*nz,nx,1,    idir)
      call FFT(a,b,nz   ,ny,nx,   idir)
      call FFT(a,b,1    ,nz,nx*ny,idir)

      i=0
      do 22 iz=0,nz-1
      do 22 iy=0,ny-1
      do 22 ix=0,nx-1
      if(idir.lt.0) then
        imgf(ix,iy,iz)=cmplx(a(i),b(i))
      else
        imgs(ix,iy,iz)=cmplx(a(i),b(i))
      endif
   22 i=i+1

      call free(pa)
      call free(pb)
      return
      end
