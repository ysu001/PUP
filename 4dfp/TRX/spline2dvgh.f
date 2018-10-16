cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
c Washington University, Mallinckrodt Institute of Radiology.
c All Rights Reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of Washington University. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$Log: spline2dvgh.f,v $
c Revision 1.6  2007/04/22  03:36:59  avi
c correct error in splint2dvgh_testg
c eliminate splint2dvgh_testf
c
c Revision 1.5  2007/04/01  00:41:44  avi
c eliminate use of rand()
c
c Revision 1.4  2007/03/30  20:42:36  avi
c prepare for f90
c
c Revision 1.3  2000/12/13  02:55:07  avi
c copyright
c
c Revision 1.2  1998/05/25  09:16:11  avi
c splint2dvgh_testf
c
c Revision 1.1  1998/04/20  02:44:53  avi
c Initial revision
c
c$Header: /data/petsun4/data1/src_solaris/TRX/RCS/spline2dvgh.f,v 1.6 2007/04/22 03:36:59 avi Exp $

c     call spline2dvgh_rcs
c     call splint2dv_test
c     call splint2dvgh_testg
c     call splint2dvgh_testh
c     end

      subroutine spline2dvgh_rcs
      write (*,"('$Id: spline2dvgh.f,v 1.6 2007/04/22 03:36:59 avi Exp $')")
      return
      end

      subroutine splint2dv_test
      parameter (nx=15)
      parameter (ny=9)
      real*4 f(0:nx-1,0:ny-1)
      real*4 d2xf(0:nx-1,0:ny-1)
      real*4 d2yf(0:nx-1,0:ny-1)
      real*4 g(-nx/2:nx/2,-ny/2:ny/2)
      real*4 d(-nx/2:nx/2)
      data del/0.2/

      write (*,"('splint2dv_test')")
      ix0=nx/2
      iy0=ny/2
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
    2 f(ix,iy)=0.
      f(ix0,iy0)=1.

      call splinex (f,nx,ny,1,d2xf)
      call spliney (f,nx,ny,1,d2yf)

      do 3 iy=-ny/2,ny/2
      do 3 ix=-nx/2,nx/2
      y=float(iy0)+del*float(iy)
      x=float(ix0)+del*float(ix)
      call splint2dv(f,nx,ny,d2xf,d2yf,x,y,v)
    3 g(ix,iy)=v

      do 4 ix=-nx/2,nx/2
    4 d(ix)=float(ix0)+del*float(ix)

      write (*,"('g')")
      write (*,"(6x,17f5.2)")d
      write (*,"(6x,17('-----'))")
      do iy=-ny/2,ny/2
        write (*,"(f5.2,'|',17f5.2)")float(iy0)+del*float(iy),(g(ix,iy),ix=-nx/2,nx/2)
      enddo

      write (*,"('d2xf')")
      do iy=0,ny-1
        write (*,"(i3,17f5.2)")iy,(d2xf(ix,iy),ix=0,nx-1)
      enddo

      write (*,"('d2yf')")
      do iy=0,ny-1
        write (*,"(i3,17f5.2)")iy,(d2yf(ix,iy),ix=0,nx-1)
      enddo

      return
      end

      subroutine splint2dv(imag,nx,ny,d2xi,d2yi,x,y,v)
      real*4 imag(0:nx-1,0:ny-1),d2xi(0:nx-1,0:ny-1)
      real*4 d2yi(0:nx-1,0:ny-1)

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

      v=
     &+ax*(ay*imag(ix,iy)+wy*imag(ix,iy1))+wx*(ay*imag(ix1,iy)+wy*imag(ix1,iy1))
     &+cx*(ay*d2xi(ix,iy)+wy*d2xi(ix,iy1))+dx*(ay*d2xi(ix1,iy)+wy*d2xi(ix1,iy1))
     &+ax*(cy*d2yi(ix,iy)+dy*d2yi(ix,iy1))+wx*(cy*d2yi(ix1,iy)+dy*d2yi(ix1,iy1))
      return
      end

      subroutine splint2dvgh_testh
      parameter (nx=15)
      parameter (ny=9)
      real*4 f(0:nx-1,0:ny-1)
      real*4 d2xf(0:nx-1,0:ny-1)
      real*4 d2yf(0:nx-1,0:ny-1)
      real*4 x(2),g0(2),gx(2),gy(2),hess(2,2),test(2,2)
      data del/0.001/

      write (*,"('splint2dvgh_testh')")
      ix0=nx/2
      iy0=ny/2
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
    2 f(ix,iy)=0.
      f(ix0,iy0)=1.

      call splinex (f,nx,ny,1,d2xf)
      call spliney (f,nx,ny,1,d2yf)

      do 11 l=-5,5
      x(1)=float(ix0)+float(l)/10.+0.05
      x(2)=float(iy0)+0.5
      call splint2dvgh(f,nx,ny,d2xf,d2yf,x(1),x(2),v,g0,hess)

      call splint2dvg(f,nx,ny,d2xf,d2yf,x(1)+del,x(2),    v,gx)
      call splint2dvg(f,nx,ny,d2xf,d2yf,x(1),    x(2)+del,v,gy)
      do k=1,2
        test(1,k)=(gx(k)-g0(k))/del
        test(2,k)=(gy(k)-g0(k))/del
      enddo
      do i=1,2
        write(*,"(3(2f10.6,5x))")(hess(i,j),j=1,2),(test(i,j),j=1,2),(hess(i,j)-test(i,j),j=1,2)
      enddo
      call splint2dvg(f,nx,ny,d2xf,d2yf,x(1)-del,x(2),    v,gx)
      call splint2dvg(f,nx,ny,d2xf,d2yf,x(1),    x(2)-del,v,gy)
      do k=1,2
        test(1,k)=-(gx(k)-g0(k))/del
        test(2,k)=-(gy(k)-g0(k))/del
      enddo
      do i=1,2
        write(*,"(3(2f10.6,5x))")(hess(i,j),j=1,2),(test(i,j),j=1,2),(hess(i,j)-test(i,j),j=1,2)
      enddo
   11 write(*,"()")

      return
      end

      subroutine splint2dvgh_testg
      parameter (nx=15)
      parameter (ny=9)
      real*4 f(0:nx-1,0:ny-1)
      real*4 d2xf(0:nx-1,0:ny-1)
      real*4 d2yf(0:nx-1,0:ny-1)
      real*4 x(2),grad(2),test(2)
      data del/0.001/

      write(*,"('splint2dvgh_testg')")
      ix0=nx/2
      iy0=ny/2
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
    2 f(ix,iy)=0.
      f(ix0,iy0)=1.

      call splinex (f,nx,ny,1,d2xf)
      call spliney (f,nx,ny,1,d2yf)

      do 11 l=-5,5
      x(1)=float(ix0)+float(l)/10.+0.05
      x(2)=float(iy0)+0.5

      call splint2dv (f,nx,ny,d2xf,d2yf,x(1)+del,x(2),    gx)
      call splint2dv (f,nx,ny,d2xf,d2yf,x(1),    x(2)+del,gy)
      call splint2dv (f,nx,ny,d2xf,d2yf,x(1)-del,x(2),    ux)
      call splint2dv (f,nx,ny,d2xf,d2yf,x(1),    x(2)-del,uy)
      call splint2dvg(f,nx,ny,d2xf,d2yf,x(1),    x(2),     v,grad)

      test(1)=0.5*(gx-ux)/del
      test(2)=0.5*(gy-uy)/del
   11 write (*,"(10f10.4)")x,v,grad,(grad(i)-test(i),i=1,2)

      return
      end

      subroutine splint2dvg(imag,nx,ny,d2xi,d2yi,x,y,v,grad)
      real*4 imag(0:nx-1,0:ny-1),d2xi(0:nx-1,0:ny-1)
      real*4 d2yi(0:nx-1,0:ny-1)
      real*4 grad(3)

      if(nx.lt.2.or.ny.lt.2)then
        write(*,"('splint2dvg: x and y image dimensions must be at least 2')")
        call exit(-2)
      endif

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
      ex=(1.-3.*ax*ax)/6.
      fx=(3.*wx*wx-1.)/6.

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
      ey=(1.-3.*ay*ay)/6.
      fy=(3.*wy*wy-1.)/6.

      v=
     &+ax*(ay*imag(ix,iy)+wy*imag(ix,iy1))+wx*(ay*imag(ix1,iy)+wy*imag(ix1,iy1))
     &+cx*(ay*d2xi(ix,iy)+wy*d2xi(ix,iy1))+dx*(ay*d2xi(ix1,iy)+wy*d2xi(ix1,iy1))
     &+ax*(cy*d2yi(ix,iy)+dy*d2yi(ix,iy1))+wx*(cy*d2yi(ix1,iy)+dy*d2yi(ix1,iy1))

      grad(1)=
     &   -(ay*imag(ix,iy)+wy*imag(ix,iy1))   +(ay*imag(ix1,iy)+wy*imag(ix1,iy1))
     &+ex*(ay*d2xi(ix,iy)+wy*d2xi(ix,iy1))+fx*(ay*d2xi(ix1,iy)+wy*d2xi(ix1,iy1))
     &   -(cy*d2yi(ix,iy)+dy*d2yi(ix,iy1))   +(cy*d2yi(ix1,iy)+dy*d2yi(ix1,iy1))
      grad(2)=
     &+ax*(  -imag(ix,iy)   +imag(ix,iy1))+wx*(  -imag(ix1,iy)   +imag(ix1,iy1))
     &+cx*(  -d2xi(ix,iy)   +d2xi(ix,iy1))+dx*(  -d2xi(ix1,iy)   +d2xi(ix1,iy1))
     &+ax*(ey*d2yi(ix,iy)+fy*d2yi(ix,iy1))+wx*(ey*d2yi(ix1,iy)+fy*d2yi(ix1,iy1))

      return
      end

      subroutine splint2dvgh(imag,nx,ny,d2xi,d2yi,x,y,v,grad,hess)
      real*4 imag(0:nx-1,0:ny-1),d2xi(0:nx-1,0:ny-1)
      real*4 d2yi(0:nx-1,0:ny-1)
      real*4 grad(2),hess(2,2)

      if(nx.lt.2.or.ny.lt.2)then
        write(*,"('splint2dvgh: x and y image dimensions must be at least 2')")
        call exit(-2)
      endif

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
      ex=(1.-3.*ax*ax)/6.
      fx=(3.*wx*wx-1.)/6.

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
      ey=(1.-3.*ay*ay)/6.
      fy=(3.*wy*wy-1.)/6.

      v=
     &+ax*(ay*imag(ix,iy)+wy*imag(ix,iy1))+wx*(ay*imag(ix1,iy)+wy*imag(ix1,iy1))
     &+cx*(ay*d2xi(ix,iy)+wy*d2xi(ix,iy1))+dx*(ay*d2xi(ix1,iy)+wy*d2xi(ix1,iy1))
     &+ax*(cy*d2yi(ix,iy)+dy*d2yi(ix,iy1))+wx*(cy*d2yi(ix1,iy)+dy*d2yi(ix1,iy1))

      grad(1)=
     &   -(ay*imag(ix,iy)+wy*imag(ix,iy1))   +(ay*imag(ix1,iy)+wy*imag(ix1,iy1))
     &+ex*(ay*d2xi(ix,iy)+wy*d2xi(ix,iy1))+fx*(ay*d2xi(ix1,iy)+wy*d2xi(ix1,iy1))
     &   -(cy*d2yi(ix,iy)+dy*d2yi(ix,iy1))   +(cy*d2yi(ix1,iy)+dy*d2yi(ix1,iy1))
      grad(2)=
     &+ax*(  -imag(ix,iy)   +imag(ix,iy1))+wx*(  -imag(ix1,iy)   +imag(ix1,iy1))
     &+cx*(  -d2xi(ix,iy)   +d2xi(ix,iy1))+dx*(  -d2xi(ix1,iy)   +d2xi(ix1,iy1))
     &+ax*(ey*d2yi(ix,iy)+fy*d2yi(ix,iy1))+wx*(ey*d2yi(ix1,iy)+fy*d2yi(ix1,iy1))

      hess(1,1)=
     &+ax*(ay*d2xi(ix,iy)+wy*d2xi(ix,iy1))+wx*(ay*d2xi(ix1,iy)+wy*d2xi(ix1,iy1))
      hess(1,2)=
     &   -(  -imag(ix,iy)   +imag(ix,iy1))   +(  -imag(ix1,iy)   +imag(ix1,iy1))
     &+ex*(  -d2xi(ix,iy)   +d2xi(ix,iy1))+fx*(  -d2xi(ix1,iy)   +d2xi(ix1,iy1))
     &   -(ey*d2yi(ix,iy)+fy*d2yi(ix,iy1))   +(ey*d2yi(ix1,iy)+fy*d2yi(ix1,iy1))
      hess(2,2)=
     &+ax*(ay*d2yi(ix,iy)+wy*d2yi(ix,iy1))+wx*(ay*d2yi(ix1,iy)+wy*d2yi(ix1,iy1))
      hess(2,1)=hess(1,2)

      return
      end

