c $Id: frmsmri.f,v 1.7 2009/02/25 01:04:00 avi Exp $
c $Log: frmsmri.f,v $
c Revision 1.7  2009/02/25  01:04:00  avi
c recode quiet NaN to accommodate 64 bit architecture
c
c Revision 1.6  2007/08/08  02:33:52  avi
c r_quiet_nan() -> '7fffffff'x for gcc compatibility
c
c Revision 1.5  2007/04/17  05:35:59  avi
c gcc compliant
c
c Revision 1.4  1998/03/20  03:23:07  avi
c correct '7fffffff'x -> r_quiet_nan()
c
c Revision 1.3  1997/10/11  00:38:03  avi
c inhibit parameter initialization ('1000'x) bit
c
c Revision 1.2  1996/04/19  21:53:31  avi
c Revision 1.1  1996/04/19  17:08:34  ty7777
c Initial revision
c
c
      subroutine alignm(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode,err)
c     alignm finds optimal linear transform of img2 to effect realignment for on img1
c     the computations apply in the nonzero part of mask

      character*256 rcsid /'$Header: /data/petsun4/data1/src_solaris/librms/RCS/frmsmri.f,v 1.7 2009/02/25 01:04:00 avi Exp $'/
      real*4 img1(0:nx-1,0:ny-1,0:nz-1),img2(0:nx-1,0:ny-1,0:nz-1)
      integer*2 mask(0:nx-1,0:ny-1,0:nz-1)
      real*4 param(12),mmppix(3),s4(4)
      real*4 taram(12)
      parameter (dd=4.0)		! linear distance search increment in units of mmppix
      parameter (da=0.0349066)		! 2 deg = multiplicative (angle in radians) search increment
      real*4 daram(12,2)/dd,dd,da,da,da,da,0.,0.,0.,0.,0.,0.,
     &                   dd,dd,dd,da,da,da,da,da,da,da,da,da/
      real*4 coef(0:2),array(3,3)
      parameter (nimax=6)
      real*4 x(-nimax:nimax),y(-nimax:nimax)
c     real*4 report(-nimax:nimax,12)
      logical*4 liter,lfar
      logical*4 lenable,l3d,lstretch,lsgrad,lverbose,linit
      logical*4 ldebug/.false./

      if(nz.lt.2)mode=iand(mode,not(2))	! prevent attempt to 3d align
      lenable=iand(mode,1).ne.0
      l3d=iand(mode,2).ne.0
      lstretch=iand(mode,4).ne.0
      lsgrad=iand(mode,8).ne.0
      method=iand(mode,384)/128
      lverbose=iand(mode,512).eq.0
      linit=iand(mode,'1000'x).eq.0
      write(*,"('image alignment mode ',i6,' decimal ',z8,' hex')")mode,mode

      if(linit)then
        do i=1,12				! initialize param
          param(i)=0.
        enddo
      endif

      if(lenable)then
        i=0
        if(lstretch)i=i+1
        if(l3d)then
          i=i+1
          jndex=2
        else
          jndex=1
        endif
        jmax=3*2**i
      else
        jmax=0
      endif

      call alignma(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode,1)	! initialize s4
      if(.not.lenable)goto 90

      if(lverbose)then
        write(*,"('niter, ni, mesh',3i4)")niter,0,1
        call alignme(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode,1,err,q)
        call alignml(param,jmax,s4,err)
      endif
      niter=5
      ni=1
      mesh=4
      nnn=0
   83 liter=.false.
      lfar=.false.
      if(niter.le.3)ni=max0(8-niter*2,2)
      mesh=max0(min0(niter,4),2)
      if(lverbose)write(*,"('niter, ni, mesh',3i4)")niter,ni,mesh
      do 81 j=1,jmax
      taram(j)=param(j)
      do 84 i=-ni,ni
      param(j)=taram(j)+daram(j,jndex)*float(i)/float(ni)
      call alignme(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode,mesh,err,q)
      if(ldebug)call alignml(param,jmax,s4,err)
      x(i)=float(i)/float(ni)
   84 y(i)=err
      call polfit(x(-ni),y(-ni),2*ni+1,coef,array,3,chisqr)
      t=.5*coef(1)/coef(2)
      if(coef(2).lt.0..or.abs(t).gt.0.5)then
        t0=t
        t=sign(.5,coef(1))
        if(lverbose)write(*,"(f10.6,' ->',f10.6,' parameter',i3)")t0,t,j
        lfar=.true.
      endif
      liter=liter.or.abs(t).gt.0.001
   81 param(j)=taram(j)-daram(j,jndex)*t
      call alignme(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode,mesh,err,q)
      if(lsgrad)then
        call alignma(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode,mesh)
      else
        s4(1)=q*s4(1)
      endif
      if(lverbose)call alignml(param,jmax,s4,err)
      nnn=nnn+1
      if(nnn.gt.75)then
        write(*,"('failed realignment convergence')")
        return
      endif
      if(.not.lfar)niter=niter-1
      if(liter.and.niter.gt.0)goto 83

   90 if(lverbose)then
        write(*,"('niter, ni, mesh',3i4)")0,0,1
        call alignme(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode,1,err,q)
        call alignml(param,jmax,s4,err)
      endif
      return
      end

      subroutine alignml(param,jmax,s4,err)
      real*4 param(jmax),s4(4)
      if(jmax.gt.0)write(*,"(6f10.4)")(param(j),j=1,jmax)
      write(*,"(4f10.4,10x,f10.4)")s4,sqrt(err)
      return
      end

      subroutine alignme(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode,mesh,err,q)
      real*4 img1(0:nx-1,0:ny-1,0:nz-1),img2(0:nx-1,0:ny-1,0:nz-1)
      integer*2 mask(0:nx-1,0:ny-1,0:nz-1)
      real*4 param(12),s4(4),mmppix(3)
      real*4 center(3),a(4,4),b(4,4),c(4,4),t4(4,4)
      logical*4 ldebug/.false./

      center(1)=mmppix(1)*float(nx/2)
      center(2)=mmppix(2)*float(ny/2)
      center(3)=mmppix(3)*float(nz/2)
      call img2vrt(mmppix,center,c)
      call param2warp(mode,param,b)
      call matmul(b,c,a,4)
      call vrt2img(mmppix,center,b)
      call matmul(b,a,t4,4)

      n=0
      sum2=0.
      sum1=0.
      method=iand(mode,384)/128
      goto (20,20,40,30),method+1

   20 err=0.			! minimize dif image variance
      do 24 k=0,nz-1
      do 24 j=0,ny-1,mesh
      do 24 i=0,nx-1,mesh
      if(mask(i,j,k).eq.0)goto 24
      call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i,j,k,v,lslice)
      if(lslice.le.0)goto 24
      n=n+1
      sum2=sum2+v
      u=img1(i,j,k)
      sum1=sum1+u
      err=err+(v-u)**2
   24 continue
      q=sum1/sum2
      sum2=sum2/float(n)
      sum1=sum1/float(n)
      err=err/float(n)
      if(ldebug)write(*,"('n,sum1,sum2,rms',i6,3f10.1)")n,sum1,sum2,sqrt(err)
      return

   30 u2p=0.			! image variance accumulator;	zero crossing theory method
      up2=0.			! dif image |grad| variance accumulator
      do 34 k=0,nz-1
      do 34 j=0,ny-1,mesh
      do 34 i=0,nx-1,mesh
      if(mask(i  ,j,  k).eq.0)goto 34
      if(mask(i-1,j,  k).eq.0)goto 34
      if(mask(i  ,j-1,k).eq.0)goto 34
      if(nz.gt.1.and.mask(i,j,k-1).eq.0)goto 34
      call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i  ,j  ,k,v ,lslice)
      if(lslice.le.0)goto 34
      call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i-1,j  ,k,vx,lslice)
      if(lslice.le.0)goto 34
      call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i  ,j-1,k,vy,lslice)
      if(lslice.le.0)goto 34
      if(nz.gt.1)then
        call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i,j,k-1,vz,lslice)
        if(lslice.le.0)goto 34
      endif
      n=n+1
      sum2=sum2+v
      u =img1(i  ,j  ,k)
      sum1=sum1+u
      u2p=u2p+(v-u)**2
      ux=img1(i-1,j  ,k)
      uy=img1(i  ,j-1,k)
      if(nz.gt.1)then
        uz=img1(i,j,k-1)
        up2=up2+(v-u-(vx-ux))**2+(v-u-(vy-uy))**2+(v-u-(vz-uz))**2
      else
        up2=up2+(v-u-(vx-ux))**2+(v-u-(vy-uy))**2
      endif
   34 continue
      sum2=sum2/float(n)
      sum1=sum1/float(n)
      u2p=u2p/float(n)
      up2=up2/float(n)
      err=10000.*u2p/up2
      q=sum1/sum2
      return

   40 u11=0.			! |grad| covariance accumulator;	gradient matching method
      u20=0.			! image1 |grad| variance accumulator
      u02=0.			! image2 |grad| variance accumulator
      do 44 k=0,nz-1
      do 44 j=0,ny-1,mesh
      do 44 i=0,nx-1,mesh
      if(mask(i  ,j,  k).eq.0)goto 44
      if(mask(i-1,j,  k).eq.0)goto 44
      if(mask(i  ,j-1,k).eq.0)goto 44
      if(nz.gt.1.and.mask(i,j,k-1).eq.0)goto 44
      call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i  ,j  ,k,v ,lslice)
      if(lslice.le.0)goto 44
      call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i-1,j  ,k,vx,lslice)
      if(lslice.le.0)goto 44
      call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i  ,j-1,k,vy,lslice)
      if(lslice.le.0)goto 44
      if(nz.gt.1)then
        call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i,j,k-1,vz,lslice)
        if(lslice.le.0)goto 44
      endif
      u =img1(i  ,j  ,k)
      ux=img1(i-1,j  ,k)
      uy=img1(i  ,j-1,k)
      if(nz.gt.1)uz=img1(i,j,k-1)
      if(nz.gt.1)then
        a11=(u-ux)*(u-ux)+(u-uy)*(u-uy)+(u-uz)*(u-uz)
        a12=(u-ux)*(v-vx)+(u-uy)*(v-vy)+(u-uz)*(v-vz)
        a22=(v-vx)*(v-vx)+(v-vy)*(v-vy)+(v-vz)*(v-vz)
      else
        a11=(u-ux)*(u-ux)+(u-uy)*(u-uy)
        a12=(u-ux)*(v-vx)+(u-uy)*(v-vy)
        a22=(v-vx)*(v-vx)+(v-vy)*(v-vy)
      endif
      if(a11*a22.gt.0.)then
        n=n+1
        sum1=sum1+u
        sum2=sum2+v
        u11=u11+a12*a12/sqrt(a11*a22)
        u20=u20+a11
        u02=u02+a22
      endif
   44 continue
      q=sum1/sum2
      sum2=sum2/float(n)
      sum1=sum1/float(n)
      u20=u20/float(n)
      u02=u02/float(n)
      u11=u11/float(n)
      err=10000.*(1.-u11/sqrt(u20*u02))
      if(ldebug)write(*,"('n,sum1,sum2,u20,u02,u11,rms',i6,6f10.1)")n,sum1,sum2,u20,u02,u11,sqrt(err)
      return

      end

      subroutine alignma(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode,mesh)
      real*4 img1(0:nx-1,0:ny-1,0:nz-1),img2(0:nx-1,0:ny-1,0:nz-1)
      integer*2 mask(0:nx-1,0:ny-1,0:nz-1)
      real*4 param(12),s4(4),mmppix(3)
      real*4 center(3),a(4,4),b(4,4),c(4,4),t4(4,4)
      real*4 sarray(4,4),scaleb(4),x4(4)
      logical*4 lsgrad
      logical*4 ldebug/.false./

      lsgrad=iand(mode,8).ne.0

      center(1)=mmppix(1)*float(nx/2)
      center(2)=mmppix(2)*float(ny/2)
      center(3)=mmppix(3)*float(nz/2)
      call img2vrt(mmppix,center,c)
      call param2warp(mode,param,b)
      call matmul(b,c,a,4)
      call vrt2img(mmppix,center,b)
      call matmul(b,a,t4,4)

      do 41 l=1,4
      scaleb(l)=0.
      s4(l)=0.
      do 41 m=1,4
   41 sarray(l,m)=0.
      s4(1)=1.				! start out with scale gradient compensation off

      n=0
      do 24 k=0,nz-1
      do 24 j=0,ny-1,mesh
      do 24 i=0,nx-1,mesh
      if(mask(i,j,k).eq.0)goto 24
      call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i,j,k,v,lslice)
      if(lslice.gt.0)then
        n=n+1
        u=img1(i,j,k)
        x4(1)=1.			! x4 = [1 x y z] (not [x y z 1])
        do l=1,3			! translate target indices to source coordinates
          x4(l+1)=a(l,1)*float(i)+a(l,2)*float(j)+a(l,3)*float(k)+a(l,4)
        enddo
        do l=1,4
          scaleb(l)=scaleb(l)+x4(l)*u
          do m=1,4
            sarray(l,m)=sarray(l,m)+x4(l)*x4(m)*v
          enddo
        enddo
      endif
   24 continue
      if(lsgrad)then
        if(sarray(4,4).eq.0.)sarray(4,4)=sarray(3,3)	! allows matrix inversion when nz = 1
        if(ldebug)then
          write(*,"('gradient normal equations before inversion')")
          do l=1,4
            write(*,"(4e14.4,10x,e14.4)")(sarray(l,m),m=1,4),scaleb(l)
          enddo
        endif
        call matinv(sarray,4,det)
        do l=1,4
          s4(l)=0.
          do m=1,4
            s4(l)=s4(l)+sarray(l,m)*scaleb(m)
          enddo
        enddo
      else
        s4(1)=scaleb(1)/sarray(1,1)
        do l=2,4
          s4(l)=0.
        enddo
      endif
      if(ldebug)write(*,"('alignma scale coefficients',4f10.6)")s4
      return
      end

      subroutine alignmd(nx,ny,nz,img1,img2,mask,param,s4,mmppix,mode)
c     compute (+) img1 (-) [realigned] img2 and leave in img1
      real*4 img1(0:nx-1,0:ny-1,0:nz-1),img2(0:nx-1,0:ny-1,0:nz-1)
      integer*2 mask(0:nx-1,0:ny-1,0:nz-1)
      real*4 param(12),s4(4),mmppix(3)
      real*4 center(3),a(4,4),b(4,4),c(4,4),t4(4,4)
      logical*4 ldebug/.false./

      center(1)=mmppix(1)*float(nx/2)
      center(2)=mmppix(2)*float(ny/2)
      center(3)=mmppix(3)*float(nz/2)
      call img2vrt(mmppix,center,c)
      call param2warp(mode,param,b)
      call matmul(b,c,a,4)
      call vrt2img(mmppix,center,b)
      call matmul(b,a,t4,4)

      n=0
      u=0.
      do 24 k=0,nz-1
      do 24 j=0,ny-1
      do 24 i=0,nx-1
      if(mask(i,j,k).eq.0)then
        img1(i,j,k)=0.
      else
        call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i,j,k,v,lslice)
        if(lslice.le.0)then
          mask(i,j,k)=0
          img1(i,j,k)=0.
        else
          n=n+1
          img1(i,j,k)=img1(i,j,k)-v
          u=u+img1(i,j,k)
        endif
      endif
   24 continue
      u=u/float(n)
      write(*,"('dif image voxels',i8,'  mean',f10.4)")n,u
      return
      end

      subroutine alignmr(nx,ny,nz,img2,imgr,mask,param,s4,mmppix,mode)
      real*4 img2(0:nx-1,0:ny-1,0:nz-1),imgr(0:nx-1,0:ny-1,0:nz-1)
      integer*2 mask(0:nx-1,0:ny-1,0:nz-1)
      real*4 param(12),s4(4),mmppix(3)
      real*4 center(3),a(4,4),b(4,4),c(4,4),t4(4,4)
      real*4 rnan
      integer*4 irnan
      equivalence(rnan,irnan)
      logical*4 ldebug/.false./
      irnan=x'7fffffff'

      center(1)=mmppix(1)*float(nx/2)
      center(2)=mmppix(2)*float(ny/2)
      center(3)=mmppix(3)*float(nz/2)
      call img2vrt(mmppix,center,c)
      call param2warp(mode,param,b)
      call matmul(b,c,a,4)
      call vrt2img(mmppix,center,b)
      call matmul(b,a,t4,4)

      n=0
      u=0.
      do 24 k=0,nz-1
      do 24 j=0,ny-1
      do 24 i=0,nx-1
      call imgvals(nx,ny,nz,img2,t4,s4,a,mode,i,j,k,v,lslice)
      if(lslice.le.0)then
        mask(i,j,k)=0
c       imgr(i,j,k)=0.
        imgr(i,j,k)=rnan	! quiet NaN
      else
        n=n+1
        imgr(i,j,k)=v
        u=u+v
      endif
   24 continue
      u=u/float(n)
      write(*,"('realigned image voxels',i8,'  mean',f10.4)")n,u
      return
      end

      subroutine param2warp(mode,param,a)
      real*4 a(4,4),param(12)
      logical*4 lenable,l3d,lstretch

      lenable=iand(mode,1).ne.0
      l3d=iand(mode,2).ne.0
      lstretch=iand(mode,4).ne.0

      do l=1,4
        do m=1,4
          a(l,m)=0.
        enddo
        a(l,l)=1.
      enddo

      if(lenable)then
        if(lstretch)then
          if(l3d)then
            lmax=3
          else
            lmax=2
          endif
          jj=1
          do l=1,lmax
            a(l,4)=param(jj)
            jj=jj+1
          enddo
          do l=1,lmax
            a(l,l)=param(jj)+1.0
            jj=jj+1
          enddo
          do l=1,lmax
            do m=1,lmax
              if(m.ne.l)then
                a(l,m)=param(jj)
                jj=jj+1
              endif
            enddo
          enddo
        else
          if(l3d)then
            call trotset(param(1),param(4),a)
          else
            a(1,4)=param(1)
            a(2,4)=param(2)
            a(1,1)= cos(param(3))
            a(1,2)=-sin(param(3))
            a(2,1)= sin(param(3))
            a(2,2)= cos(param(3))
          endif
        endif
      endif

      return
      end

      subroutine warp2param(mode,a,param)
      real*4 a(4,4),param(12),rot(3,3)
      logical*4 lenable,l3d,lstretch

      lenable=iand(mode,1).ne.0
      l3d=iand(mode,2).ne.0
      lstretch=iand(mode,4).ne.0

      do m=1,12
        param(m)=0.
      enddo

      if(lenable)then
        if(lstretch)then
          if(l3d)then
            lmax=3
          else
            lmax=2
          endif
          jj=1
          do l=1,lmax
            param(jj)=a(l,4)
            jj=jj+1
          enddo
          do l=1,lmax
            param(jj)=a(l,l)-1.0
            jj=jj+1
          enddo
          do l=1,lmax
            do m=1,lmax
              if(m.ne.l)then
                param(jj)=a(l,m)
                jj=jj+1
              endif
            enddo
          enddo
        else
          if(l3d)then
            do i=1,3
              param(i)=a(i,4)
              do j=1,3
                rot(i,j)=a(i,j)
              enddo
            enddo
            call rot2ang(rot,param(4))
          else
            param(1)=a(1,4)
            param(2)=a(2,4)
            param(3)=atan2(a(2,1),a(2,2))
          endif
        endif
      endif

      return
      end

      real*4 function determ12(a,det)
      real*4 a(4,4)

      det=+a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     &    -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
     &    +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      determ12=det
      return
      end

      subroutine imgvals(nx,ny,nz,imag,t4,s4,a4,mode,i,j,k,v,lslice)
c     returns in v the interpolated value of imag at transformed coordinates (i,j,k)
c     lslice returns the slice from which most of the data is taken.
c     lslice = -1 (and v = 0.) if the image is undefined at transformed (i,j,k).
      real*4 imag(0:nx-1,0:ny-1,0:nz-1)
      real*4 t4(4,4),s4(4),a4(4,4)
      real*4 x(3),xp(3),x4(4)
      logical*4 lwrap,lendslc
      logical*4 ldebug/.false./

      lwrap=iand(mode,1024).ne.0
      lendslc=iand(mode,2048).ne.0
      x(1)=float(i)
      x(2)=float(j)
      x(3)=float(k)
      x4(1)=1.
      do l=1,3
        xp(l)=  t4(l,1)*x(1)+t4(l,2)*x(2)+t4(l,3)*x(3)+t4(l,4)
        x4(l+1)=a4(l,1)*x(1)+a4(l,2)*x(2)+a4(l,3)*x(3)+a4(l,4)
      enddo
      scale=s4(1)*x4(1)+s4(2)*x4(2)+s4(3)*x4(3)+s4(4)*x4(4)

      ix=nint(xp(1)-.5)
      wx=xp(1)-float(ix)
      if(lwrap)then
        ix=mod(ix+nx,nx)
        ix1=mod(ix+1,nx)
      else
        if(ix.eq.-1.and.wx.gt.0.999)then
          ix=0
          wx=0.
        else if(ix.eq.nx-1.and.wx.lt.0.001)then
          ix=nx-2
          wx=1.
        endif
        ix1=ix+1
        if(ix.lt.0.or.ix1.gt.nx-1)goto 9
      endif
      iy=nint(xp(2)-.5)
      wy=xp(2)-float(iy)
      if(lwrap)then
        iy=mod(iy+ny,ny)
        iy1=mod(iy+1,ny)
      else
        if(iy.eq.-1.and.wy.gt.0.999)then
          iy=0
          wy=0.
        else if(iy.eq.ny-1.and.wy.lt.0.001)then
          iy=ny-2
          wy=1.
        endif
        iy1=iy+1
        if(iy.lt.0.or.iy1.gt.ny-1)goto 9
      endif
      if(nz.eq.1)goto 10
      iz=nint(xp(3)-.5)
      wz=xp(3)-float(iz)
      if(lwrap)then
        iz=mod(iz+nz,nz)
        iz1=mod(iz+1,nz)
      else
        if(iz.eq.-1.and.(lendslc.or.wz.gt.0.999))then
          iz=0
          wz=0.
        else if(iz.eq.nz-1.and.(lendslc.or.wz.lt.0.001))then
          iz=nz-2
          wz=1.
        endif
        iz1=iz+1
        if(iz.lt.0.or.iz1.gt.nz-1)goto 9
      endif
      v=scale*(
     & +(1.-wz)*((1.-wy)*(
     &               (1.-wx)*imag(ix ,iy ,iz )
     &                   +wx*imag(ix1,iy ,iz ))
     &               +wy*(
     &               (1.-wx)*imag(ix ,iy1,iz )
     &                   +wx*imag(ix1,iy1,iz )))
     &      +wz*((1.-wy)*(
     &               (1.-wx)*imag(ix ,iy ,iz1)
     &                   +wx*imag(ix1,iy ,iz1))
     &               +wy*(
     &               (1.-wx)*imag(ix ,iy1,iz1)
     &                   +wx*imag(ix1,iy1,iz1)))
     &)
      lslice=iz+1
      if(wz.gt.0.5)lslice=iz1+1
      if((.not.v.ge.0.).and.(.not.v.lt.0.))write(*,"(6i4,f10.4)")ix,ix1,iy,iy1,iz,iz1,v
      return
  10  v=scale*(
     &           (1.-wy)*(
     &               (1.-wx)*imag(ix ,iy ,0)
     &                   +wx*imag(ix1,iy ,0))
     &               +wy*(
     &               (1.-wx)*imag(ix ,iy1,0)
     &                   +wx*imag(ix1,iy1,0))
     &)
      lslice=1
      return
    9 v=0.
      lslice=-1
      return
      end
