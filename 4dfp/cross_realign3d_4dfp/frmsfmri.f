c $Header: /data/petsun4/data1/src_solaris/cross_realign3d_4dfp/RCS/frmsfmri.f,v 1.4 2011/09/01 01:55:03 avi Exp $
c $Log: frmsfmri.f,v $
c Revision 1.4  2011/09/01  01:55:03  avi
c move up initialization of niter in alignmaiz() to prevent reporting value as undefined when lverbose.eq..true.
c
c Revision 1.3  2007/08/08  00:32:53  avi
c gcc compliant
c integer*2 -> integer*4
c
c Revision 1.2  2006/09/28  21:46:06  avi
c determ12 now returns det
c
c Revision 1.1  1997/05/23  01:30:29  yang
c Initial revision
c
c Revision 1.8  1996/05/17  20:54:55  avi
c alignml s4 format 4f10.4 => 4f10.6
c
c Revision 1.7  1996/05/13  07:45:27  avi
c correct do 44 loop in alignmeiz
c
c Revision 1.6  1996/05/13  05:42:32  avi
c correct start indices of do 34 loop in alignmeiz
c
c Revision 1.5  1996/05/12  23:50:48  avi
c add "alignmriz: " to alignmriz stdout message.
c
c Revision 1.4  1996/04/30  20:53:16  avi
c more variable name changes
c
c Revision 1.3  1996/04/30  00:53:31  avi
c variable name changes in preparation for introducing alignmrixyz
c
c Revision 1.2  1996/04/25  22:07:25  ty7777
c Fix bug in lslice computation.
c
c Revision 1.1  1996/04/25  22:04:51  ty7777
c Initial revision
c
      subroutine alignmiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode,err)
c     alignmiz finds optimal linear transform of img2 to effect realignment for on img1
c     the computations apply in the nonzero part of mask
c     img2 is read with interpolation in z using subroutines splinez and splintz

      character*256 rcshdr
     &/'$Id: frmsfmri.f,v 1.4 2011/09/01 01:55:03 avi Exp $'/
      real*4    img1(0:nx-1,0:ny-1,0:nz-1),img2(0:nx-1,0:ny-1,0:nz-1),d2img2(0:nx-1,0:ny-1,0:nz-1)
      integer*4 mask(0:nx-1,0:ny-1,0:nz-1)
      real*4 param(12),mmppix(3),s4(4)
      real*4 taram(12)
      parameter (dd=1.0)		! linear distance search range in mm
      parameter (da=0.0174533)		! 1 deg = multiplicative (angle in radians) search range
      real*4 daram(12,2)/dd,dd,da,da,da,da,0.,0.,0.,0.,0.,0.,
     &                   dd,dd,dd,da,da,da,da,da,da,da,da,da/
      real*4 coef(0:2),array(3,3)
      parameter (nimax=6)
      real*4 x(-nimax:nimax),y(-nimax:nimax)
      real*4 report(-nimax:nimax,12)
      logical*4 liter,lfar
      logical*4 lenable,l3d,lstretch,lsgrad,lverbose
      logical*4 ldebug/.false./

      if(nz.lt.2)mode=iand(mode,not(2))	! prevent attempt to 3d align
      lenable=iand(mode,1).ne.0
      l3d=iand(mode,2).ne.0
      lstretch=iand(mode,4).ne.0
      lsgrad=iand(mode,8).ne.0
      method=iand(mode,384)/128
      lverbose=iand(mode,512).eq.0
      write(*,"('image alignment mode ',i6,' decimal ',z8,' hex')")mode,mode

      do i=1,12				! initialize param
        param(i)=0.
      enddo

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

      call alignmaiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode,1)	! initialize s4
      if(.not.lenable)goto 90

      niter=5
      if(lverbose)then
        write(*,"('niter, ni, mesh',3i4)")niter,0,1
        call alignmeiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode,1,err,q)
        call alignml(param,jmax,s4,err)
      endif
      ni=1
      mesh=4
      nnn=0
   83 liter=.false.
      lfar=.false.
      if(niter.le.2)ni=2
      mesh=max0(min0(niter,5),3)
      if(lverbose)write(*,"('niter, ni, mesh',3i4)")niter,ni,mesh
      do 81 j=1,jmax
      taram(j)=param(j)
      do 84 i=-ni,ni
      param(j)=taram(j)+daram(j,jndex)*float(i)/float(ni)
      call alignmeiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode,mesh,err,q)
      if(ldebug)call alignml(param,jmax,s4,err)
      x(i)=float(i)/float(ni)
      report(i,j)=sqrt(err)
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
      call alignmeiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode,mesh,err,q)
      if(lsgrad)then
        call alignmaiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode,mesh)
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
        write(*,"(6('       +/-'))")
        write(*,"(6f10.4)")(daram(j,jndex),j=1,jmax)
        do i=-ni,ni
          write(*,"(6f10.4)")(report(i,j),j=1,jmax)
        enddo
        write(*,"('niter, ni, mesh',3i4)")0,0,1
        call alignmeiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode,1,err,q)
        call alignml(param,jmax,s4,err)
      endif
      return
      end

      subroutine alignml(param,jmax,s4,err)
      real*4 param(jmax),s4(4)
      if(jmax.gt.0)write(*,"(6f10.4)")(param(j),j=1,jmax)
      write(*,"(4f10.6,10x,f10.4)")s4,sqrt(err)
      return
      end

      subroutine alignmeiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode,mesh,err,q)
      real*4    img1(0:nx-1,0:ny-1,0:nz-1),img2(0:nx-1,0:ny-1,0:nz-1),d2img2(0:nx-1,0:ny-1,0:nz-1)
      integer*4 mask(0:nx-1,0:ny-1,0:nz-1)
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
      call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i,j,k,v,lslice)
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
      if(ldebug)write(*, "('n,sum1,sum2,rms',i6,3f10.1)")n,sum1,sum2,sqrt(err)
      return

   30 u2p=0.			! image variance accumulator;	zero crossing theory method
      up2=0.			! dif image |grad| variance accumulator
      do 34 k=1,nz-1
      do 34 j=1,ny-1,mesh
      do 34 i=1,nx-1,mesh
      if(mask(i  ,j,  k).eq.0)goto 34
      if(mask(i-1,j,  k).eq.0)goto 34
      if(mask(i  ,j-1,k).eq.0)goto 34
      if(nz.gt.1.and.mask(i,j,k-1).eq.0)goto 34
      call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i  ,j  ,k,v ,lslice)
      if(lslice.le.0)goto 34
      call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i-1,j  ,k,vx,lslice)
      if(lslice.le.0)goto 34
      call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i  ,j-1,k,vy,lslice)
      if(lslice.le.0)goto 34
      if(nz.gt.1)then
        call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i,j,k-1,vz,lslice)
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
      do 44 k=1,nz-1
      do 44 j=1,ny-1,mesh
      do 44 i=1,nx-1,mesh
      if(mask(i  ,j,  k).eq.0)goto 44
      if(mask(i-1,j,  k).eq.0)goto 44
      if(mask(i  ,j-1,k).eq.0)goto 44
      if(nz.gt.1.and.mask(i,j,k-1).eq.0)goto 44
      call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i  ,j  ,k,v ,lslice)
      if(lslice.le.0)goto 44
      call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i-1,j  ,k,vx,lslice)
      if(lslice.le.0)goto 44
      call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i  ,j-1,k,vy,lslice)
      if(lslice.le.0)goto 44
      if(nz.gt.1)then
        call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i,j,k-1,vz,lslice)
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
      if(ldebug)write(*, "('n,sum1,sum2,u20,u02,u11,rms',i6,6f10.1)")n,sum1,sum2,u20,u02,u11,sqrt(err)
      return

      end

      subroutine alignmaiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode,mesh)
      real*4    img1(0:nx-1,0:ny-1,0:nz-1),img2(0:nx-1,0:ny-1,0:nz-1),d2img2(0:nx-1,0:ny-1,0:nz-1)
      integer*4 mask(0:nx-1,0:ny-1,0:nz-1)
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
      call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i,j,k,v,lslice)
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
          write (*, "('gradient normal equations before inversion')")
          do l=1,4
            write (*, "(4e14.4,10x,e14.4)")(sarray(l,m),m=1,4),scaleb(l)
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

      subroutine alignmdiz(nx,ny,nz,img1,img2,d2img2,mask,param,s4,mmppix,mode)
c     compute (+) img1 (-) [realigned] img2 and leave in img1
      real*4    img1(0:nx-1,0:ny-1,0:nz-1),img2(0:nx-1,0:ny-1,0:nz-1),d2img2(0:nx-1,0:ny-1,0:nz-1)
      integer*4 mask(0:nx-1,0:ny-1,0:nz-1)
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
        call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i,j,k,v,lslice)
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

      subroutine alignmriz(nx,ny,nz,img2,d2img2,imgr,mask,param,s4,mmppix,mode)
      real*4    img2(0:nx-1,0:ny-1,0:nz-1),imgr(0:nx-1,0:ny-1,0:nz-1),d2img2(0:nx-1,0:ny-1,0:nz-1)
      integer*4 mask(0:nx-1,0:ny-1,0:nz-1)
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
      call imgvalsiz(nx,ny,nz,img2,d2img2,t4,s4,a,mode,i,j,k,v,lslice)
      if(lslice.le.0)then
        mask(i,j,k)=0
        imgr(i,j,k)=0.
c       imgr(i,j,k)='7fffffff'x	! NaN
      else
        n=n+1
        imgr(i,j,k)=v
        u=u+v
      endif
   24 continue
      u=u/float(n)
      write(*,"('alignmriz: realigned image voxels',i8,'  mean',f10.4)")n,u
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

      determ12 = det
      end

      subroutine imgvalsiz(nx,ny,nz,imag,d2zi,t4,s4,a4,mode,i,j,k,v,lslice)
c     returns in v the interpolated value of imag at transformed coordinates (i,j,k)
c     lslice returns the slice from which most of the data is taken.
c     lslice = -1 (and v = 0.) if the image is undefined at transformed (i,j,k).
c     imag is read with interpolation in z using splintz logic
c     d2zi must have been computed by a call to splinez
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
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
        dowhile(ix.lt.0)
          ix=ix+nx
        enddo
        ix=mod(ix,nx)
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
        dowhile(iy.lt.0)
          iy=iy+ny
        enddo
        iy=mod(iy,ny)
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
      az=1.-wz
      if(iz.lt.0)then
        if(.not.lendslc)goto 9
        iz=0
        az=1.
        wz=0.
      endif
      if(iz.ge.nz-1)then
        if(.not.lendslc)goto 9
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
      v=scale*(
     &           (1.-wy)*(
     &               (1.-wx)*v00
     &                   +wx*v10)
     &               +wy*(
     &               (1.-wx)*v01
     &                   +wx*v11)
     &)
      lslice=iz+1+nint(wz)
      if((.not.v.ge.0.).and.(.not.v.lt.0.))then
        write(*,"('i,j,k,x,y,z,xp,yp,zp',3i4,6f10.4)")i,j,k,x,xp
        write(*,"('t4')")
        do l=1,4
          write(*,"(3f10.6,f10.4)")(t4(l,m),m=1,4)
        enddo
        write(*,"('ix,ix1,iy,iy1,iz,iz+1,v',6i4,4x,2f10.4)")ix,ix1,iy,iy1,iz,iz+1,v
        write(*,"('wx,wy,wz,az,cz,dz',6f10.6)")wx,wy,wz,az,cz,dz
        write(*,"(2i4,'|',2i4,4f10.4)")ix ,iy ,iz,iz+1,imag(ix ,iy ,iz+0),imag(ix ,iy ,iz+1),d2zi(ix ,iy ,iz+0),d2zi(ix ,iy ,iz+1)
        write(*,"(2i4,'|',2i4,4f10.4)")ix1,iy ,iz,iz+1,imag(ix1,iy ,iz+0),imag(ix1,iy ,iz+1),d2zi(ix1,iy ,iz+0),d2zi(ix1,iy ,iz+1)
        write(*,"(2i4,'|',2i4,4f10.4)")ix ,iy1,iz,iz+1,imag(ix ,iy1,iz+0),imag(ix ,iy1,iz+1),d2zi(ix ,iy1,iz+0),d2zi(ix ,iy1,iz+1)
        write(*,"(2i4,'|',2i4,4f10.4)")ix1,iy1,iz,iz+1,imag(ix1,iy1,iz+0),imag(ix1,iy1,iz+1),d2zi(ix1,iy1,iz+0),d2zi(ix1,iy1,iz+1)
        call exit(-2)
      endif
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
