c$Header: /data/petsun4/data1/src_solaris/nonlin/RCS/b2v4m.f,v 1.3 2007/07/04 03:00:53 avi Exp $ 
c$Log: b2v4m.f,v $
c Revision 1.3  2007/07/04  03:00:53  avi
c voxdim initialization gcc compliant
c
c Revision 1.2  2004/11/05  22:12:51  rsachs
c Changed dimensions (nx,ny,nz).
c
c Revision 1.1  2004/10/28  21:49:18  rsachs
c Initial revision
c
c 2004.10.22 - Reintroducing voxdim into the calculation of the derivatives.
c------------------------------------------------------------------------------------      
      subroutine b2v4_test
      parameter (nx=192,ny=224,nz=192)
      real*4 b(0:nx-1,0:ny-1,0:nz-1,3),v(0:nx-1,0:ny-1,0:nz-1,3),g(0:nx-1,0:ny-1,0:nz-1,3)
      real*4 berr(3),voxdim(3)
      real*4 mu/1./
      real*4 lambda/2./ 
      integer*4 iflag/0/
c     data voxdim/1.25,1.50,1.75/
      data voxdim/1.50,1.75,2.50/
      pi=atan2(0.0,-1.0)
      write(*,"('pi=',f10.8)")pi
      pfac=2.0*pi

      do 1 iv=1,3
      u=0.
      do 2 iz=0,nz-1
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
      b(ix,iy,iz,iv)=rand(iflag)
    2 u=u+b(ix,iy,iz,iv)
      u=u/float(nx*ny*nz)
      do 3 iz=0,nz-1
      do 3 iy=0,ny-1
      do 3 ix=0,nx-1
    3 b(ix,iy,iz,iv)=b(ix,iy,iz,iv)-u
    1 continue

      call b2v4(b,nx,ny,nz,voxdim,mu,lambda,v)
      call  v2b(g,nx,ny,nz,voxdim,mu,lambda,v)

      errm=0.
      do 22 iv=1,3
      write(*,"('b vs g volume',i2)")iv
      do 22 jz=1,nz-2
      do 22 jy=1,ny-2
      do 22 jx=1,nx-2
      do 24 i=1,3
   24 berr(i)=b(jx,jy,jz,iv)+g(jx,jy,jz,iv)
      q=berr(1)*berr(1)+berr(2)*berr(2)+berr(3)*berr(3)
      if(q.gt.errm)then
        errm=q
        jxm=jx
        jym=jy
        jzm=jz
        ivm=iv
      endif
c     write(*,"(3i5,6f10.6)")jx,jy,jz,b(jx,jy,jz,iv),g(jx,jy,jz,iv),b(jx,jy,jz,iv)+g(jx,jy,jz,iv)
   22 continue
      write(*,"('greatest b vs g error found in volume',i2)")ivm
      write(*,"(3i5,6f10.6)")jxm,jym,jzm,b(jxm,jym,jzm,ivm),g(jxm,jym,jzm,ivm),b(jxm,jym,jzm,ivm)+g(jxm,jym,jzm,ivm)
      return
      end

      subroutine b2v4(b,nx,ny,nz,voxdim,mu,lambda,v)
c     calculate the velocity field (v) from the 
c     body force field (b) using complex Fourier representation

      real*4 b(0:nx-1,0:ny-1,0:nz-1,3),v(0:nx-1,0:ny-1,0:nz-1,3),voxdim(3),mu,lambda
      complex*8 c(0:nx-1,0:ny-1,0:nz-1),ftv(0:nx-1,0:ny-1,0:nz-1,3),ftb(0:nx-1,0:ny-1,0:nz-1,3)
      pointer (pc,c),(pftv,ftv),(pftb,ftb)

      pc=  malloc(8*nx*ny*nz)
      pftb=malloc(8*nx*ny*nz*3)
      pftv=malloc(8*nx*ny*nz*3)
      if(pftb.eq.0.or.pftv.eq.0.or.pc.eq.0) stop 'b2v4 malloc failure'

      do 21 iv=1,3
      do 22 iz=0,nz-1
      do 22 iy=0,ny-1
      do 22 ix=0,nx-1
   22 c(ix,iy,iz)=cmplx(b(ix,iy,iz,iv),0.0)
c     write(*,"('b2v4 forward FFT',i1)"),iv
      call cmplxfft3d(c,ftb(0,0,0,iv),nx,ny,nz,-1)	! c -> ftb
   21 continue

      call ftb2v(ftv,ftb,nx,ny,nz,voxdim,mu,lambda)

      do 31 iv=1,3
c     write(*,"('b2v4 inverse FFT',i1)"),iv
      call cmplxfft3d(c,ftv(0,0,0,iv),nx,ny,nz,+1)	! c <- ftv
      do 32 iz=0,nz-1
      do 32 iy=0,ny-1
      do 32 ix=0,nx-1
   32 v(ix,iy,iz,iv)=real(c(ix,iy,iz))
   31 continue

      call free(pc)
      call free(pftv)
      call free(pftb)
      return
      end

      subroutine ftv2b(ftv,ftb,nx,ny,nz,voxdim,mu,lambda)
      complex*8 ftb(0:nx-1,0:ny-1,0:nz-1,3),ftv(0:nx-1,0:ny-1,0:nz-1,3)
      real*4 voxdim(3),mu,lambda
      real*4 a(3,3),cfac(3),sfac(3)

      pfac=2.0*atan2(0.0,-1.0)
      do 11 kz=0,nz-1
      f=pfac*float(kz)/float(nz)
      cfac(3)=2.*(cos(f)-1.)/(voxdim(3)*voxdim(3))
      sfac(3)=sin(f)/voxdim(3)
      do 11 ky=0,ny-1
      f=pfac*float(ky)/float(ny)
      cfac(2)=2.*(cos(f)-1.)/(voxdim(2)*voxdim(2))
      sfac(2)=sin(f)/voxdim(2)
      do 11 kx=0,nx-1
      f=pfac*float(kx)/float(nx)
      cfac(1)=2.*(cos(f)-1.)/(voxdim(1)*voxdim(1))
      sfac(1)=sin(f)/voxdim(1)
      delfac=cfac(1)+cfac(2)+cfac(3)
      do 12 i=1,3
      do 12 j=1,3
      if(i.eq.j)then
        a(i,j)=mu*delfac+(mu+lambda)*cfac(i)
      else
        a(i,j)=-(mu+lambda)*sfac(i)*sfac(j)
      endif
   12 continue
      do 14 iv=1,3
      ftb(kx,ky,kz,iv)=cmplx(0.,0.)
      do 14 jv=1,3
   14 ftb(kx,ky,kz,iv)=ftb(kx,ky,kz,iv)-a(iv,jv)*ftv(kx,ky,kz,jv)
   11 continue

      return
      end

      subroutine ftb2v(ftv,ftb,nx,ny,nz,voxdim,mu,lambda)
      complex*8 ftb(0:nx-1,0:ny-1,0:nz-1,3),ftv(0:nx-1,0:ny-1,0:nz-1,3)
      real*4 voxdim(3),mu,lambda
      real*4 a(3,3),ainv(3,3),cfac(3),sfac(3)

      pfac=2.0*atan2(0.0,-1.0)
      do 11 kz=0,nz-1
      f=pfac*float(kz)/float(nz)
      cfac(3)=2.*(cos(f)-1.)/(voxdim(3)*voxdim(3))
      sfac(3)=sin(f)/voxdim(3)
      do 11 ky=0,ny-1
      f=pfac*float(ky)/float(ny)
      cfac(2)=2.*(cos(f)-1.)/(voxdim(2)*voxdim(2))
      sfac(2)=sin(f)/voxdim(2)
      do 11 kx=0,nx-1
      f=pfac*float(kx)/float(nx)
      cfac(1)=2.*(cos(f)-1.)/(voxdim(1)*voxdim(1))
      sfac(1)=sin(f)/voxdim(1)
      delfac=cfac(1)+cfac(2)+cfac(3)
      do 12 i=1,3
      do 12 j=1,3
      if(i.eq.j)then
        a(i,j)=mu*delfac+(mu+lambda)*cfac(i)
      else
        a(i,j)=-(mu+lambda)*sfac(i)*sfac(j)
      endif
   12 continue
      do 13 iv=1,3
   13 ftv(kx,ky,kz,iv)=cmplx(0.,0.)
      if(kz.eq.0.and.ky.eq.0.and.kx.eq.0)goto 11
      call mat3inv(a,ainv,det)
      do 14 iv=1,3
      do 14 jv=1,3
   14 ftv(kx,ky,kz,iv)=ftv(kx,ky,kz,iv)-ainv(iv,jv)*ftb(kx,ky,kz,jv)
   11 continue

      return
      end

      subroutine v2b(b,nx,ny,nz,voxdim,mu,lambda,v)
c     calculate the velocity field (v) from the body force field (b) in the spatial domain

      real*4 b(0:nx-1,0:ny-1,0:nz-1,3),v(0:nx-1,0:ny-1,0:nz-1,3),voxdim(3),mu,lambda
      real*4 ddv(6,3)

      do 1 iz=1,nz-2
      do 1 iy=1,ny-2
      do 1 ix=1,nx-2
      do 2 iv=1,3
      ddv(1,iv)=(v(ix+1,iy,iz,iv)+v(ix-1,iy,iz,iv)-2.*v(ix,iy,iz,iv))/(voxdim(1)*voxdim(1))
      ddv(2,iv)=(v(ix,iy+1,iz,iv)+v(ix,iy-1,iz,iv)-2.*v(ix,iy,iz,iv))/(voxdim(2)*voxdim(2))
      ddv(3,iv)=(v(ix,iy,iz+1,iv)+v(ix,iy,iz-1,iv)-2.*v(ix,iy,iz,iv))/(voxdim(3)*voxdim(3))
      ddv(4,iv)=(v(ix+1,iy+1,iz,iv)+v(ix-1,iy-1,iz,iv)-v(ix+1,iy-1,iz,iv)-v(ix-1,iy+1,iz,iv))*0.25/(voxdim(1)*voxdim(2))
      ddv(5,iv)=(v(ix+1,iy,iz+1,iv)+v(ix-1,iy,iz-1,iv)-v(ix+1,iy,iz-1,iv)-v(ix-1,iy,iz+1,iv))*0.25/(voxdim(1)*voxdim(3))
    2 ddv(6,iv)=(v(ix,iy+1,iz+1,iv)+v(ix,iy-1,iz-1,iv)-v(ix,iy+1,iz-1,iv)-v(ix,iy-1,iz+1,iv))*0.25/(voxdim(2)*voxdim(3))
      b(ix,iy,iz,1)=mu*(ddv(1,1)+ddv(2,1)+ddv(3,1))+(mu+lambda)*(ddv(1,1)+ddv(4,2)+ddv(5,3))
      b(ix,iy,iz,2)=mu*(ddv(1,2)+ddv(2,2)+ddv(3,2))+(mu+lambda)*(ddv(4,1)+ddv(2,2)+ddv(6,3))
      b(ix,iy,iz,3)=mu*(ddv(1,3)+ddv(2,3)+ddv(3,3))+(mu+lambda)*(ddv(5,1)+ddv(6,2)+ddv(3,3))
    1 continue

      return
      end
