c$Header: /data/petsun4/data1/src_solaris/crop_4dfp/RCS/fzero_slice.f,v 1.2 2010/08/20 04:43:23 avi Exp $
c$Log: fzero_slice.f,v $
c Revision 1.2  2010/08/20  04:43:23  avi
c eliminate slice limit testing in FORTRAN
c
c Revision 1.1  1999/06/28  22:44:34  avi
c Initial revision
c
      subroutine fzero_slice(dimstr,imgt,nx,ny,nz,istart,iend)
      character*2 dimstr
      real*4 imgt(nx,ny,nz)

      if(dimstr(1:1).eq.'z'.or.dimstr(1:1).eq.'H')goto 10
      if(dimstr(1:1).eq.'y'.or.dimstr(1:1).eq.'C')goto 20
      if(dimstr(1:1).eq.'x'.or.dimstr(1:1).eq.'S')goto 30
      write(*,"('fzero_slice: direction character not one of [xyzCHS]')")
      call exit (-1)

   10 do 11 iz=istart,iend
      do 11 iy=1,ny
      do 11 ix=1,nx
   11 imgt(ix,iy,iz)=0.
      return

   20 do 21 iz=1,nz
      do 21 iy=istart,iend
      do 21 ix=1,nx
   21 imgt(ix,iy,iz)=0.
      return

   30 do 31 iz=1,nz
      do 31 iy=1,ny
      do 31 ix=istart,iend
   31 imgt(ix,iy,iz)=0.
      return

      end
