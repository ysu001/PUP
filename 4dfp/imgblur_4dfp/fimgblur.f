c $Header: /data/petsun4/data1/src_solaris/imgblur_4dfp/RCS/fimgblur.f,v 1.3 2008/01/01 22:39:55 avi Exp $
c $Log: fimgblur.f,v $
c Revision 1.3  2008/01/01  22:39:55  avi
c gcc v4 compliant
c
c Revision 1.2  1999/03/10  02:23:51  avi
c subroutine imgblur3d
c
c Revision 1.1  1998/10/02  19:24:42  mcavoy
c Initial revision
c
c Revision 1.1  1997/10/30  15:38:21  tscull
c Initial revision
c
c Revision 1.1  1997/01/16  04:47:58  avi
c Initial revision
c
      subroutine fimgblur_rcs
      write(*,"('$Id: fimgblur.f,v 1.3 2008/01/01 22:39:55 avi Exp $')")
      return
      end

      subroutine imgblur3d(fwhm,alpha,voxsiz,img0,imgdim,img1)
      real*4 img0(1),img1(1),voxsiz(3)
      integer*4 imgdim(3)

      call fimgblur_rcs
      call imgblurx(fwhm,alpha,voxsiz(1),img0,imgdim(1),imgdim(2),imgdim(3),img1)
      call imgblury(fwhm,alpha,voxsiz(2),img1,imgdim(1),imgdim(2),imgdim(3),img0)
      call imgblurz(fwhm,alpha,voxsiz(3),img0,imgdim(1),imgdim(2),imgdim(3),img1)
      return
      end

      subroutine imgblurx(fwhm,alpha,voxsiz,imag,nx,ny,nz,imgb)
      real*4 imag(nx,ny,nz),imgb(nx,ny,nz)
      parameter (nkmax=256)
      real*4 h(-nkmax:nkmax)
      logical*4 ldebug/.false./

      if(fwhm.le.0..or.alpha.le.0.)then
        write(*,"('imgblurx: fwhm and alpha must be greater than 0')")
        return
      endif
      sigma=0.5*fwhm/sqrt(2.*alog(2.))
      nk=nint(sigma*sqrt(-2.*alog(alpha))/abs(voxsiz))
      if(nk.gt.nkmax)then
        write(*,"('imgblurx: nk reduced to',i5)")nkmax
        nk=nkmax
      endif
      if(ldebug) write(*,"('imgblurx: nk=',i5)")nk

      s=0.
      do 5 k=-nk,nk
      h(k)=exp(-0.5*(float(k)*voxsiz/sigma)**2)
    5 s=s+h(k)
      do 6 k=-nk,nk
    6 h(k)=h(k)/s

      do 1 iz=1,nz
      do 1 iy=1,ny
      do 1 ix=1,nx
      s=0.
      do 2 k=-nk,nk
      kx=min0(max0(ix+k,1),nx)
    2 s=s+h(k)*imag(kx,iy,iz)
    1 imgb(ix,iy,iz)=s

      return
      end

      subroutine imgblury(fwhm,alpha,voxsiz,imag,nx,ny,nz,imgb)
      real*4 imag(nx,ny,nz),imgb(nx,ny,nz)
      parameter (nkmax=256)
      real*4 h(-nkmax:nkmax)
      logical*4 ldebug/.false./

      if(fwhm.le.0..or.alpha.le.0.)then
        write(*,"('imgblurx: fwhm and alpha must be greater than 0')")
        return
      endif
      sigma=0.5*fwhm/sqrt(2.*alog(2.))
      nk=nint(sigma*sqrt(-2.*alog(alpha))/abs(voxsiz))
      if(nk.gt.nkmax)then
        write(*,"('imgblury: nk reduced to',i5)")nkmax
        nk=nkmax
      endif
      if(ldebug) write(*,"('imgblury: nk=',i5)")nk

      s=0.
      do 5 k=-nk,nk
      h(k)=exp(-0.5*(float(k)*voxsiz/sigma)**2)
    5 s=s+h(k)
      do 6 k=-nk,nk
    6 h(k)=h(k)/s

      do 1 iz=1,nz
      do 1 ix=1,nx
      do 1 iy=1,ny
      s=0.
      do 2 k=-nk,nk
      ky=min0(max0(iy+k,1),ny)
    2 s=s+h(k)*imag(ix,ky,iz)
    1 imgb(ix,iy,iz)=s

      return
      end

      subroutine imgblurz(fwhm,alpha,voxsiz,imag,nx,ny,nz,imgb)
      real*4 imag(nx,ny,nz),imgb(nx,ny,nz)
      parameter (nkmax=256)
      real*4 h(-nkmax:nkmax)
      logical*4 ldebug/.false./

      if(fwhm.le.0..or.alpha.le.0.)then
        write(*,"('imgblurz: fwhm and alpha must be greater than 0')")
        return
      endif
      sigma=0.5*fwhm/sqrt(2.*alog(2.))
      nk=nint(sigma*sqrt(-2.*alog(alpha))/abs(voxsiz))
      if(nk.gt.nkmax)then
        write(*,"('imgblurx: nk reduced to',i5)")nkmax
        nk=nkmax
      endif
      if(ldebug) write(*,"('imgblurz: nk=',i5)")nk

      s=0.
      do 5 k=-nk,nk
      h(k)=exp(-0.5*(float(k)*voxsiz/sigma)**2)
    5 s=s+h(k)
      do 6 k=-nk,nk
    6 h(k)=h(k)/s

      do 1 iy=1,ny
      do 1 ix=1,nx
      do 1 iz=1,nz
      s=0.
      do 2 k=-nk,nk
      kz=min0(max0(iz+k,1),nz)
    2 s=s+h(k)*imag(ix,iy,kz)
    1 imgb(ix,iy,iz)=s

      return
      end
