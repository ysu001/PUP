c$Header: /data/petsun4/data1/src_solaris/cross_realign3d_4dfp/RCS/alignmrixyz.f,v 1.3 2007/08/08 00:28:07 avi Exp $
c$Log: alignmrixyz.f,v $
c Revision 1.3  2007/08/08  00:28:07  avi
c integer*2 -> integer*4
c
c Revision 1.2  2004/05/09  00:58:17  avi
c output undefined voxels as 1.0e-37
c
c Revision 1.1  1997/05/23  01:30:29  yang
c Initial revision
c
c Revision 1.5  1996/05/14  02:48:18  avi
c put rcshdr on continuation line for FORTRAN compiler
c
c Revision 1.4  1996/05/14  02:34:42  avi
c add Id, Log and Header fields
c
      subroutine alignmrixyz(nx,ny,nz,imag,d2xi,d2yi,d2zi,imgr,mdef,param,s4,mmppix,mode)
      character*256 rcshdr
     &/'$Id: alignmrixyz.f,v 1.3 2007/08/08 00:28:07 avi Exp $'/
      real*4    imag(0:nx-1,0:ny-1,0:nz-1),imgr(0:nx-1,0:ny-1,0:nz-1)
      real*4    d2xi(0:nx-1,0:ny-1,0:nz-1),d2yi(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
      integer*4 mdef(0:nx-1,0:ny-1,0:nz-1)
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
      call imgvalsixyz(nx,ny,nz,imag,d2xi,d2yi,d2zi,t4,s4,a,mode,i,j,k,v,lslice)
      if(lslice.le.0)then
        mdef(i,j,k)=0
        imgr(i,j,k)=1.0e-37
      else
        n=n+1
        imgr(i,j,k)=v
        u=u+v
      endif
   24 continue
      u=u/float(n)
      write(*,"('alignmrixyz: realigned image voxels',i8,'  mean',f10.4)")n,u
      return
      end
