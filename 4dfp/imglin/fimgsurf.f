c$Header: /data/petsun4/data1/src_solaris/imglin/RCS/fimgsurf.f,v 1.1 2007/09/05 04:04:43 avi Exp $
c$Log: fimgsurf.f,v $
c Revision 1.1  2007/09/05  04:04:43  avi
c Initial revision
c
      subroutine fimgsurf(imgt,nx,ny,nz,center,mmppix,x0,grad,xopt)
c     param(01:03)	normal to plane of edge trace
c     param(04)		edge trace arclength increment (mm)
c     param(05)		edge trace gradient sampling radius (mm)
      real*4 imgt(nx,ny,nz)
      real*4 center(3),mmppix(3),x0(3),grad(3),xopt(3)
      real*4 unorm(3),e(3)
      real*4 param(16)
      logical*4 ldebug/.false./

      param(05)=4.	! half distance (mm) for gradient optimizatiun

      call grad_opt(imgt,nx,ny,nz,mmppix,center,param,x0,grad,xopt,ierr)

c     compute unit normal
      a=sqrt(dot(grad,grad))
      do i=1,3
        unorm(i)=-grad(i)/a
      enddo

      do i=1,3
        e(i)=x0(i)-xopt(i)
      enddo
      dist=dot(unorm,e)

      if (ldebug) then
      write(*, "('ierr       ',i10)")ierr
      write(*, "('x0         ',3f10.4)")x0
      write(*, "('xopt       ',3f10.4)")xopt
      write(*, "('unit normal',3f10.4)")unorm
      write(*, "('dist       ',3f10.4)")dist
      write(*, "('|grad|     ',3f10.4)")a
      endif

      return
      end
