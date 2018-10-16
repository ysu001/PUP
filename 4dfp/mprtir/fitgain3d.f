c$Header: /data/petsun4/data1/src_solaris/mprtir/RCS/fitgain3d.f,v 1.6 2008/01/02 01:22:55 avi Exp $
c$Log: fitgain3d.f,v $
c Revision 1.6  2008/01/02  01:22:55  avi
c gcc v4 compliant
c
c Revision 1.5  1999/02/18  09:10:00  avi
c again improve speed
c new logic to prevent g -> 0
c
c Revision 1.4  1999/02/18  01:22:54  avi
c improve speed
c Revision 1.3  1999/02/15  10:06:15  avi
c print rcsid
      subroutine fitgain3d_test
      parameter (nx=64,ny=64,nz=17)
      real*4    imgf(0:nx-1,0:ny-1,0:nz-1)
      integer*2 mask(0:nx-1,0:ny-1,0:nz-1)
      parameter (norder=2)
      parameter (nterm=((norder+3)*(norder+2)*(norder+1))/6)
      real*4    a(nterm),p(nterm)
      parameter (nc=10)
      real*4    fv(nc),fu(nc)
      integer*4 nv(nc)

      do i=1,nc
        fv(i)=float(i)*rand(0)
        fu(i)=0.
        nv(i)=0
      enddo

      a(1)=1.
      do i=2,nterm
        a(i)=0.15*(1.-2.*rand(0))
      enddo
      write(*,"(10f10.6)")a

      do 1 iz=0,nz-1
      z=2.*float(iz-nz/2)/float(nz)
      do 1 iy=0,ny-1
      y=2.*float(iy-ny/2)/float(ny)
      do 1 ix=0,nx-1
      x=2.*float(ix-nx/2)/float(nx)
      k=irand(0)
      mask(ix,iy,iz)=mod(k,nc)+1
      k=mask(ix,iy,iz)
      call poly3ev(x,y,z,norder,a,p,g)
      fu(k)=fu(k)+fv(k)*g
      nv(k)=nv(k)+1
    1 imgf(ix,iy,iz)=fv(k)*g
      do i=1,nc
        fu(i)=fu(i)/float(nv(i))
      enddo

      call fitgain3d(imgf,nx,ny,nz,mask,fu,nc,norder,p)
      write(*,"(10f10.6)")(a(i)     ,i=1,nterm)
      write(*,"(10f10.6)")(p(i)/p(1),i=1,nterm)
      write(*,"(10f10.6)")(p(i)     ,i=1,nterm)

      call exit(0)
      end

      subroutine fitgain3d(imgf,nx,ny,nz,mask,fv,nc,norder,a)
      real*4 imgf(0:nx-1,0:ny-1,0:nz-1),fv(nc),a(((norder+3)*(norder+2)*(norder+1))/6)
      integer*2 mask(0:nx-1,0:ny-1,0:nz-1)
      real*4 t(((norder+3)*(norder+2)*(norder+1))/6,((norder+3)*(norder+2)*(norder+1))/6,3)
      pointer (pt,t)
      real*4 b(120),da(120),p(120),at(120)
      logical*4 liter
      logical*4 ldebug/.false./

      write(*,"('$Id: fitgain3d.f,v 1.6 2008/01/02 01:22:55 avi Exp $')")
      if(norder.gt.7)then
        write(*,"('fitgain3d: norder limit 7')")
        call exit(-1)
      endif
      nterm=((norder+3)*(norder+2)*(norder+1))/6
      pt=malloc(4*nterm**2*3)

      if(ldebug)then
        write(*,"('fitgain3d: image dimensions ',3i6)")nx,ny,nz
        do i=1,nc
          write(*,"(i10,f10.4)")i,fv(i)
        enddo
      endif

      err0=0.
      iter=0
   10 do i=1,nterm
        b(i)=0.
        do j=1,nterm
          t(i,j,1)=0.
        enddo
      enddo

      err2=0.
      fmean=0.
      nvox=0
      gmin=1.
      p(1)=1.
      do 1 iz=0,nz-1
      z=2.*float(iz-nz/2)/float(nz)
      p(2)=z
      p(5)=z*z
      do 1 iy=0,ny-1
      y=2.*float(iy-ny/2)/float(ny)
      p(3)=y
      p(6)=z*y
      p(7)=y*y
      g0=a(1)+a(2)*p(2)+a(3)*p(3)+a(5)*p(5)+a(6)*p(6)+a(7)*p(7)
      do 1 ix=0,nx-1
      x=2.*float(ix-nx/2)/float(nx)
      ic=mask(ix,iy,iz)
      if(ic.lt.1.or.ic.gt.nc)goto 1
      nvox=nvox+1
      f=imgf(ix,iy,iz)
      fmean=fmean+f
      if(norder.eq.2)then
        p(4)=x
        p(8)=z*x
        p(9)=y*x
        p(10)=x*x
        g=g0+a(4)*p(4)+a(8)*p(8)+a(9)*p(9)+a(10)*p(10)
      else
        call poly3ev(x,y,z,norder,a,p,g)
      endif
      gmin=amin1(g,gmin)
      g=amax1(g,0.15)
      err=f/g-fv(ic)
      err2=err2+err*err
c     type "(12f8.4)",x,y,z,(p(i),i=1,nterm)
      q=f/g**2
      do 3 i=1,nterm
      b(i)=b(i)+err*q*p(i)
      do 3 j=i,nterm
    3 t(i,j,1)=t(i,j,1)+q*q*p(i)*p(j)
    1 continue
      if(nvox.eq.0)return
      do i=2,nterm
        do j=1,i-1
          t(i,j,1)=t(j,i,1)
        enddo
      enddo

      q=1.0/float(nvox)
      fmean=fmean*q
      err=sqrt(err2*q)
      do 31 i=1,nterm
      b(i)=b(i)*q/fmean**2
      do 31 j=1,nterm
   31 t(i,j,1)=t(i,j,1)*q/fmean**2
      write(*,"('fitgain3d: iter',i3,'  nvox',i10,'  mean',f10.4,'  rmserr',f10.4,'  gmin',f10.6)")
     &      iter,nvox,fmean,err,gmin

      if(ldebug)then
        do i=1,nterm
          if(nterm.eq.4) write(*,"( 4f8.4,5x,f10.6,2x,f10.6)")(t(i,j,1),j=1,nterm),b(i),a(i)
          if(nterm.eq.10)write(*,"(10f8.4,5x,f10.6,2x,f10.6)")(t(i,j,1),j=1,nterm),b(i),a(i)
        enddo
      endif
      call matcop(t(1,1,1),t(1,1,2),nterm)
      call eigen(t(1,1,2),t(1,1,3),nterm)
      call matinv(t,nterm,det)
      write(*,"('fitgain3d: matrix inversion det',e12.4,' condition number',f10.4)")det,t(1,1,2)/t(nterm,nterm,2)
      do 4 i=1,nterm
      da(i)=0.
      do 4 j=1,nterm
    4 da(i)=da(i)+t(i,j,1)*b(j)
      if(ldebug)then
        do i=1,nterm
          if(nterm.eq.4) write(*,"( 4f8.2,5x,f10.6,2x,f10.6)")(t(i,j,1),j=1,nterm),b(i),da(i)
          if(nterm.eq.10)write(*,"(10f8.2,5x,f10.6,2x,f10.6)")(t(i,j,1),j=1,nterm),b(i),da(i)
        enddo
      endif
      write(*,"(10f8.4)")(a(i),i=1,nterm)

      liter=.false.
      q=1.
    6 do k=1,nterm
        if(abs(da(k)).gt.0.0002)liter=.true.
        at(k)=a(k)+q*da(k)
      enddo
      do i=1,nterm
        a(i)=at(i)
      enddo
      iter=iter+1
      if(liter.and.abs(err-err0)/fmean.gt.1.e-6)then
        err0=err
        goto 10
      endif

      call free(pt)
      return
      end

      subroutine evalgain3d(imgr,nx,ny,nz,norder,a)
      real*4 imgr(0:nx-1,0:ny-1,0:nz-1),a((norder+3)*(norder+2)*(norder+1)/6)
      real*4 p(120)
      logical*4 ldebug/.false./

      do 1 iz=0,nz-1
      z=2.*float(iz-nz/2)/float(nz)
      do 1 iy=0,ny-1
      y=2.*float(iy-ny/2)/float(ny)
      do 1 ix=0,nx-1
      x=2.*float(ix-nx/2)/float(nx)
      call poly3ev(x,y,z,norder,a,p,g)
    1 imgr(ix,iy,iz)=g

      return
      end
