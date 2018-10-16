c$Header: /data/petsun4/data1/src_solaris/imglin/RCS/grad_opt.f,v 1.1 2007/09/05 04:05:45 avi Exp $
c$Log: grad_opt.f,v $
c Revision 1.1  2007/09/05  04:05:45  avi
c Initial revision
c
      subroutine grad_opt(imgt,nx,ny,nz,mmppix,center,param,x0,grad,xopt,ierr)
      real*4 imgt(nx,ny,nz)
      real*4 mmppix(3),center(3),x0(3),grad(3),xopt(3)
      real*4 param(16)
c     param(01:03)	normal to osculating plane (used by edge_trace)
c     param(04)		del arclength              (used by edge_trace)
c     param(05)		sigma eval
      real*4 coef(0:2),array(3,3)
      parameter (nk=3)
      parameter (npts=2*nk+1)
      real*4 u(-nk:nk),v(-nk-1:nk+1),dv(-nk:nk)
      real*4 x(3),unorm(3),e(3)
      logical*4 lfar,liter
      logical*4 ldebug/.false./

c     copy x0 into xopt	as inital guess
      do i=1,3
        xopt(i)=x0(i)
      enddo
       
      do 10 ncycle=1,5
      if(ldebug)write (*,"('x0 grad_opt ',3f10.4)")x0

      nnn=0
      niter=5
   83 liter=.false.
      lfar=.false.

c     compute gradient
      do i=1,3
        do k=1,3
          x(k)=xopt(k)
        enddo
        x(i)=xopt(i)+0.5*param(05)
        call imgvalx(imgt,nx,ny,nz,center,mmppix,x,v2,lslice)
        if(lslice.le.0)goto 99
        x(i)=xopt(i)-0.5*param(05)
        call imgvalx(imgt,nx,ny,nz,center,mmppix,x,v1,lslice)
        if(lslice.le.0)goto 99
        grad(i)=(v2-v1)/param(05)
      enddo

c     compute unit normal
      a=sqrt(dot(grad,grad))
      do i=1,3
        unorm(i)=-grad(i)/a
      enddo
c     write (*, "('unit normal',3f10.4)")unorm

      do 9 k=-nk-1,nk+1			! sample image along normal
      do i=1,3
        x(i)=xopt(i)-param(05)*unorm(i)*(float(k)/float(nk))
      enddo
      call imgvalx(imgt,nx,ny,nz,center,mmppix,x,v(k),lslice)
      if(lslice.le.0)goto 99
c     write (*, "('imgvalx x,v',3f7.2,f9.1)",x,v(k)
    9 continue
      do 8 k=-nk,nk			! differentiate
      u(k)=float(k)/float(nk)
    8 dv(k)=v(k+1)-v(k-1)
c     write (*,"('dv',10f9.1)")dv

      call polfit(u,dv,npts,coef,array,3,chisqr)
c     write (*,"('coef ',3g12.4)",coef
      t=-.5*coef(1)/coef(2)
      if(coef(2).gt.0..or.abs(t).gt.0.5)then
        t0=t
        t=sign(.5,coef(1))
        if(ldebug)write(*,"('t ',f10.6,' ->',f10.6)")t0,t
        lfar=.true.
      endif
      liter=liter.or.abs(t).gt.0.02
      do i=1,3
        xopt(i)=xopt(i)-t*param(05)*unorm(i)
      enddo
      nnn=nnn+1
      if(nnn.gt.20)then
        write(*,"('failed convergence')")
        ierr=1
        return
      endif
      if(.not.lfar)then
        niter=niter-1
      endif
      if(liter.and.niter.gt.0)goto 83

c     write (*,"('xopt       ',3f10.4)")xopt
      do i=1,3
        e(i)=x0(i)-xopt(i)
      enddo
      dist=dot(unorm,e)
      do i=1,3
        xopt(i)=x0(i)-dist*unorm(i)
      enddo
c     write (*,"('xopt new   ',3f10.4)")xopt
   10 continue
      ierr=0
      return

   99 write (*,"('grad_opt: image boundary reached')")
      ierr=-1
      return
      end

