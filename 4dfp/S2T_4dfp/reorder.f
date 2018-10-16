c$Header: /data/petsun4/data1/src_solaris/S2T_4dfp/RCS/reorder.f,v 1.2 2004/03/11 07:05:26 avi Exp $
c$Log: reorder.f,v $
c Revision 1.2  2004/03/11  07:05:26  avi
c subroutine t2sf
c
c Revision 1.1  1998/12/13  04:09:29  avi
c Initial revision
c
      subroutine reorder_rcs
      write(*,"('$Id: reorder.f,v 1.2 2004/03/11 07:05:26 avi Exp $')")
      return
      end

      subroutine c2tf(imgm,nx,nz,ny,imgn)
      real*4 imgm(nx,nz,ny),imgn(nx,ny,nz)
      logical*4 ldebug/.false./

      do 1 k=1,nz
      do 1 j=1,ny
      do 1 i=1,nx
    1 imgn(i,j,k)=imgm(i,k,j)

      return
      end

      subroutine s2tf(imgm,ny,nz,nx,imgn)
      real*4 imgm(ny,nz,nx),imgn(nx,ny,nz)
      logical*4 ldebug/.false./

      do 1 k=1,nz
      do 1 j=1,ny
      do 1 i=1,nx
    1 imgn(i,j,k)=imgm(j,k,i)

      return
      end

      subroutine t2sf(imgt,nx,ny,nz,imgs)
      real*4 imgt(nx,ny,nz),imgs(ny,nz,nx)
      logical*4 ldebug/.false./

      do 1 ix=1,nx
      do 1 iz=1,nz
      do 1 iy=1,ny
    1 imgs(iy,iz,ix)=imgt(ix,iy,iz)

      return
      end
