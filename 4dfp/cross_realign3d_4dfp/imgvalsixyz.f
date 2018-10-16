c$Header: /data/petsun4/data1/src_solaris/cross_realign3d_4dfp/RCS/imgvalsixyz.f,v 1.2 2007/08/08 00:34:06 avi Exp $
c$Log: imgvalsixyz.f,v $
c Revision 1.2  2007/08/08  00:34:06  avi
c gcc compliant
c
c Revision 1.1  1997/05/23  01:30:29  yang
c Initial revision
c
c Revision 1.4  1996/06/13  02:49:27  avi
c correct typo (unmasked in wrap disabled mode) in out of bounds test
c
c Revision 1.3  1996/05/14  02:50:17  avi
c put rcshdr on continuation line for FORTRAN compiler
c
c Revision 1.2  1996/05/14  02:36:52  avi
c add Id, Log and Header fields
c
      subroutine imgvalsixyz(nx,ny,nz,imag,d2xi,d2yi,d2zi,t4,s4,a4,mode,i,j,k,v,lslice)
      character*256 rcshdr
     &/'$Id: imgvalsixyz.f,v 1.2 2007/08/08 00:34:06 avi Exp $'/
c     returns in v the interpolated value of imag at transformed coordinates (i,j,k)
c     lslice returns the slice from which most of the data is taken.
c     lslice = -1 (and v = 0.) if the image is undefined at transformed (i,j,k).
c     imag is read with cubic spline interolation in x, y and z.
c     d2xi, d2yi and d2zi must have been computed by calls to spline[xyz].
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2xi(0:nx-1,0:ny-1,0:nz-1),d2yi(0:nx-1,0:ny-1,0:nz-1)
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
      if((.not.lwrap).and.(ix.lt.0.or.ix.ge.nx-1))goto 9
      iy=nint(xp(2)-.5)
      if((.not.lwrap).and.(iy.lt.0.or.iy.ge.ny-1))goto 9
      iz=nint(xp(3)-.5)
      wz=xp(3)-float(iz)
      if(iz.lt.0)then
        if(.not.lendslc)goto 9
        iz=0
        wz=0.
      endif
      if(iz.ge.nz-1)then
        if(.not.lendslc)goto 9
        iz=nz-2
        wz=1.
      endif
      v=scale*splintxyz(imag,nx,ny,nz,d2xi,d2yi,d2zi,xp(1),xp(2),xp(3))
      if((.not.v.ge.0.).and.(.not.v.lt.0.))write(*,"(6i4,f10.4)")ix,ix1,iy,iy1,iz,iz1,v
      lslice=iz+1+nint(wz)
      return
    9 v=0.
      lslice=-1
      return
      end
