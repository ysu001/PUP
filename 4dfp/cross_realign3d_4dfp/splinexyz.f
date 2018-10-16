c$Header: /data/petsun4/data1/src_solaris/cross_realign3d_4dfp/RCS/splinexyz.f,v 1.3 2009/07/03 00:58:10 avi Exp $
c$Log: splinexyz.f,v $
c Revision 1.3  2009/07/03  00:58:10  avi
c increase z slice count limit 64->256
c
c Revision 1.2  2007/08/08  00:29:16  avi
c gcc v3 compliant
c
c Revision 1.1  1997/05/23  01:30:29  yang
c Initial revision
c
c Revision 1.5  1996/05/17  18:53:19  avi
c correct typo in spliney
c
c Revision 1.4  1996/05/14  02:46:51  avi
c put rcshdr on continuation line for FORTRAN compiler
c
c Revision 1.3  1996/05/14  02:38:38  avi
c add Id, Log and Header fields
c
      subroutine splinez_test
      character*256 rcshdr
     &/'$Id: splinexyz.f,v 1.3 2009/07/03 00:58:10 avi Exp $'/
      real*4 f(0:18)/19*0./
      real*4 d2f(0:18)
      real*4 a(19**2)

      do 1 i=0,18,2
    1 f(i)=1.0
      call splinez(f,1,1,19,d2f,a)

      do i=0,18
c       write(*, "(i10,f10.6,f10.4,f10.6)")i,f(i),d2f(i),splintz(f,1,1,19,d2f,0,0,float(i))
      enddo

      do i=-10,190
        z=0.1*float(i)
c       write(*, "(i10,2f10.6)")i,z,splintz(f,1,1,19,d2f,0,0,z)
        write(*, "(2f10.6)")z,splintz(f,1,1,19,d2f,0,0,z)
      enddo

      call exit(0)
      end

      subroutine splinexyz_test
      parameter (nx=16)
      parameter (ny=24)
      parameter (nz=9)
      real*4 f(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2zf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2xf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yf(0:nx-1,0:ny-1,0:nz-1)
      real*4 scratch(nz**2)
      real*4 g(-nx/2:nx/2,-ny/2:ny/2,-nz/2:nz/2)
      real*4 d(-nx/2:nx/2)
      data del/0.2/

      ix0=8
      iy0=8
      iz0=4
      do 2 iz=0,nz-1
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
      f(ix,iy,iz)=0.
c     if(mod(ix+iy+iz,7).eq.0)f(ix,iy,iz)=1.
    2 continue
      f(ix0,iy0,iz0)=1.

      call splinex(f,nx,ny,nz,d2xf)
      call spliney(f,nx,ny,nz,d2yf)
      call splinez(f,nx,ny,nz,d2zf,scratch)

      do 3 iz=-nz/2,nz/2
      do 3 iy=-ny/2,ny/2
      do 3 ix=-nx/2,nx/2
      z=float(iz0)+del*float(iz)
      y=float(iy0)+del*float(iy)
      x=float(ix0)+del*float(ix)
      g(ix,iy,iz)=splintxyz(f,nx,ny,nz,d2xf,d2yf,d2zf,x,y,z)
    3 continue

      do 4 ix=-nx/2,nx/2
    4 d(ix)=float(ix0)+del*float(ix)

      write(*, "('g')")
      do 51 iz=-nz/2,nz/2
      write(*, "('iz=',i2,' z=',f5.2)")iz,float(iz0)+del*float(iz)
      write(*, "(6x,17f5.2)")d
      write(*, "(6x,17('-----'))")
      do iy=-ny/2,ny/2
        write(*, "(f5.2,'|',17f5.2)")float(iy0)+del*float(iy),(g(ix,iy,iz),ix=-nx/2,nx/2)
      enddo
   51 continue
      call exit(0)

      write(*, "('d2xf')")
      do 21 iz=0,nz-1
      write(*, "('iz=',i1)")iz
      do iy=0,ny-1
        write(*, "(i3,17f5.2)")iy,(d2xf(ix,iy,iz),ix=0,nx-1)
      enddo
   21 continue

      write(*, "('d2yf')")
      do 31 iz=0,nz-1
      write(*, "('iz=',i1)")iz
      do iy=0,ny-1
        write(*, "(i3,17f5.2)")iy,(d2yf(ix,iy,iz),ix=0,nx-1)
      enddo
   31 continue

      write(*, "('d2zf')")
      do 41 iz=0,nz-1
      write(*, "('iz=',i1)")iz
      do iy=0,ny-1
        write(*, "(i3,17f5.2)")iy,(d2zf(ix,iy,iz),ix=0,nx-1)
      enddo
   41 continue

      call exit(0)

      end

      subroutine splinez(imag,nx,ny,nz,d2zi,a)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
      real*4 a(0:nz-1,0:nz-1)
      real*4 u(0:255)

      if(nz.gt.256)then
        write(*,"('splinez error: image z dimension exceeds 256')")
        call exit(-2)
      endif

      do 1 iz=0,nz-1
      do 1 jz=0,nz-1
    1 a(iz,jz)=0.
      a(0,0)=1.
      a(nz-1,nz-1)=1.
      do 2 iz=1,nz-2
      a(iz,iz-1)=1./6.
      a(iz,iz+0)=2./3.
    2 a(iz,iz+1)=1./6.
      call matinv(a,nz,det)

      do 11 ix=0,nx-1
      do 11 iy=0,ny-1
      u(0)=0.
      do iz=1,nz-2
        u(iz)=imag(ix,iy,iz+1)-2.*imag(ix,iy,iz)+imag(ix,iy,iz-1)
      enddo
      u(nz-1)=0.
      do iz=0,nz-1
        d2zi(ix,iy,iz)=0.
        do jz=0,nz-1
          d2zi(ix,iy,iz)=d2zi(ix,iy,iz)+a(iz,jz)*u(jz)
        enddo
      enddo
   11 continue

      return
      end

      subroutine splinex(imag,nx,ny,nz,d2xi)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2xi(0:nx-1,0:ny-1,0:nz-1)
      real*4 a(0:255),b(0:255),f(0:128)

      if(nx.gt.256)then
        write(*,"('splinex error: image x dimension exceeds 256')")
        call exit(-2)
      endif
      twopi=8.*atan(1.)
      do ix=0,(nx+1)/2
        c=cos(twopi*float(ix)/float(nx))
        f(ix)=6.*(c-1.)/(c+2.)
      enddo

      do 21 iz=0,nz-1
      do 21 iy=0,ny-1
      do ix=0,nx-1
        a(ix)=imag(ix,iy,iz)
        b(ix)=0.
      enddo
      call fft(a,b,1,nx,1,+1)
      a(0)=0.
      do ix=1,(nx+1)/2-1
        a(ix)=a(ix)*f(ix)
        b(ix)=b(ix)*f(ix)
        a(nx-ix)=a(nx-ix)*f(ix)
        b(nx-ix)=b(nx-ix)*f(ix)
      enddo
      if(mod(nx,2).eq.0)a(nx/2)=a(nx/2)*(-12.)
      call fft(a,b,1,nx,1,-1)
      do ix=0,nx-1
        d2xi(ix,iy,iz)=a(ix)
      enddo
   21 continue

      return
      end

      subroutine spliney(imag,nx,ny,nz,d2yi)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2yi(0:nx-1,0:ny-1,0:nz-1)
      real*4 a(0:255),b(0:255),f(0:128)

      if(ny.gt.256)then
        write(*,"('spliney error: image y dimension exceeds 256')")
        call exit(-2)
      endif
      twopi=8.*atan(1.)
      do iy=0,(ny+1)/2
        c=cos(twopi*float(iy)/float(ny))
        f(iy)=6.*(c-1.)/(c+2.)
      enddo

      do 31 iz=0,nz-1
      do 31 ix=0,nx-1
      do iy=0,ny-1
        a(iy)=imag(ix,iy,iz)
        b(iy)=0.
      enddo
      call fft(a,b,1,ny,1,+1)
      a(0)=0.
      do iy=1,(ny+1)/2-1
        a(iy)=a(iy)*f(iy)
        b(iy)=b(iy)*f(iy)
        a(ny-iy)=a(ny-iy)*f(iy)
        b(ny-iy)=b(ny-iy)*f(iy)
      enddo
      if(mod(ny,2).eq.0)a(ny/2)=a(ny/2)*(-12.)
      call fft(a,b,1,ny,1,-1)
      do iy=0,ny-1
        d2yi(ix,iy,iz)=a(iy)
      enddo
   31 continue

      return
      end

      function splintz(imag,nx,ny,nz,d2zi,ix,iy,z)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)

      iz=nint(z-0.5)
      if(iz.lt.0)then
        splintz=imag(ix,iy,0)
        return
      endif

      if(iz.gt.nz-2)then
        splintz=imag(ix,iy,nz-1)
        return
      endif

      wz=z-float(iz)
      az=1.-wz
      splintz=az*imag(ix,iy,iz)+wz*imag(ix,iy,iz+1)
     1       +((az**3-az)*d2zi(ix,iy,iz)+(wz**3-wz)*d2zi(ix,iy,iz+1))/6.
      return
      end

      function splintxyz(imag,nx,ny,nz,d2xi,d2yi,d2zi,x,y,z)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2xi(0:nx-1,0:ny-1,0:nz-1),d2yi(0:nx-1,0:ny-1,0:nz-1)

      ix=nint(x-.5)
      wx=x-float(ix)
      dowhile(ix.lt.0)
        ix=ix+nx
      enddo
      ix=mod(ix,nx)
      ix1=mod(ix+1,nx)
      ax=1.-wx
      cx=ax*(ax*ax-1.)/6.
      dx=wx*(wx*wx-1.)/6.

      iy=nint(y-.5)
      wy=y-float(iy)
      dowhile(iy.lt.0)
        iy=iy+ny
      enddo
      iy=mod(iy,ny)
      iy1=mod(iy+1,ny)
      ay=1.-wy
      cy=ay*(ay*ay-1.)/6.
      dy=wy*(wy*wy-1.)/6.

      iz=nint(z-.5)
      wz=z-float(iz)
      az=1.-wz
      if(iz.lt.0)then
        iz=0
        az=1.
        wz=0.
      endif
      if(iz.ge.nz-1)then
        iz=nz-2
        az=0.
        wz=1.
      endif
      cz=az*(az*az-1.)/6.
      dz=wz*(wz*wz-1.)/6.

      v00=az*imag(ix ,iy ,iz+0)+wz*imag(ix ,iy ,iz+1)+cz*d2zi(ix ,iy ,iz+0)+dz*d2zi(ix ,iy ,iz+1)
      v10=az*imag(ix1,iy ,iz+0)+wz*imag(ix1,iy ,iz+1)+cz*d2zi(ix1,iy ,iz+0)+dz*d2zi(ix1,iy ,iz+1)
      v01=az*imag(ix ,iy1,iz+0)+wz*imag(ix ,iy1,iz+1)+cz*d2zi(ix ,iy1,iz+0)+dz*d2zi(ix ,iy1,iz+1)
      v11=az*imag(ix1,iy1,iz+0)+wz*imag(ix1,iy1,iz+1)+cz*d2zi(ix1,iy1,iz+0)+dz*d2zi(ix1,iy1,iz+1)
      y00=az*d2yi(ix ,iy ,iz+0)+wz*d2yi(ix ,iy ,iz+1)
      y10=az*d2yi(ix1,iy ,iz+0)+wz*d2yi(ix1,iy ,iz+1)
      y01=az*d2yi(ix ,iy1,iz+0)+wz*d2yi(ix ,iy1,iz+1)
      y11=az*d2yi(ix1,iy1,iz+0)+wz*d2yi(ix1,iy1,iz+1)
      x00=az*d2xi(ix ,iy ,iz+0)+wz*d2xi(ix ,iy ,iz+1)
      x10=az*d2xi(ix1,iy ,iz+0)+wz*d2xi(ix1,iy ,iz+1)
      x01=az*d2xi(ix ,iy1,iz+0)+wz*d2xi(ix ,iy1,iz+1)
      x11=az*d2xi(ix1,iy1,iz+0)+wz*d2xi(ix1,iy1,iz+1)

      v0=ay*v00+wy*v01+cy*y00+dy*y01
      v1=ay*v10+wy*v11+cy*y10+dy*y11
      x0=ay*x00+wy*x01
      x1=ay*x10+wy*x11

      splintxyz=ax*v0+wx*v1+cx*x0+dx*x1

      return
      end
