c$Header: /data/petsun4/data1/src_solaris/img2msk_4dfp/RCS/fimg2msk.f,v 1.2 2007/09/21 22:22:19 avi Exp $
c$Log: fimg2msk.f,v $
c Revision 1.2  2007/09/21  22:22:19  avi
c gcc compliant (type -> write(*)
c
      subroutine fimg2msk_rcs
      write (*,"('$Id: fimg2msk.f,v 1.2 2007/09/21 22:22:19 avi Exp $')")
      return
      end

      subroutine fimg2msk(thresh,imag,nx,ny,nz,mask)
      real*4 imag(nx,ny,nz),mask(nx,ny,nz)
      logical*4 ldebug/.false./

      call fimg2msk_rcs
      ix1=0
      ix2=nx+1
      iy1=0
      iy2=ny+1
      iz1=0
      iz2=nz+1

      do 11 ix=1,nx
      do 12 iy=1,ny
      do 12 iz=1,nz
   12 if(imag(ix,iy,iz).gt.thresh)goto 20
   11 ix1=ix
   20 do 21 ix=nx,1,-1
      do 22 iy=1,ny
      do 22 iz=1,nz
   22 if(imag(ix,iy,iz).gt.thresh)goto 30
   21 ix2=ix

   30 do 31 iy=1,ny
      do 32 ix=1,nx
      do 32 iz=1,nz
   32 if(imag(ix,iy,iz).gt.thresh)goto 40
   31 iy1=iy
   40 do 41 iy=ny,1,-1
      do 42 ix=1,nx
      do 42 iz=1,nz
   42 if(imag(ix,iy,iz).gt.thresh)goto 50
   41 iy2=iy

   50 do 51 iz=1,nz
      do 52 ix=1,nx
      do 52 iy=1,ny
   52 if(imag(ix,iy,iz).gt.thresh)goto 60
   51 iz1=iz
   60 do 61 iz=nz,1,-1
      do 62 ix=1,nx
      do 62 iy=1,ny
   62 if(imag(ix,iy,iz).gt.thresh)goto 70
   61 iz2=iz
   70 continue

      if(ix1.gt.0)ix1=min0(ix1+1,nx)
      if(iy1.gt.0)iy1=min0(iy1+1,ny)
      if(iz1.gt.0)iz1=min0(iz1+1,nz)
      if(ix2.le.nx)ix2=max0(ix2-1,1)
      if(iy2.le.ny)iy2=max0(iy2-1,1)
      if(iz2.le.nz)iz2=max0(iz2-1,1)

      if(ix1.gt.0)then
        if(ldebug)write(*, "('zeroing x ',i3,' to ',i3)")1,ix1
        do iz=1,nz
        do iy=1,ny
        do ix=1,ix1
          mask(ix,iy,iz)=0.
        enddo
        enddo
        enddo
      endif
      if(ix2.le.nx)then
        if(ldebug)write(*, "('zeroing x ',i3,' to ',i3)")ix2,nx
        do iz=1,nz
        do iy=1,ny
        do ix=ix2,nx
          mask(ix,iy,iz)=0.
        enddo
        enddo
        enddo
      endif
        
      if(iy1.gt.0)then
        if(ldebug)write(*, "('zeroing y ',i3,' to ',i3)")1,iy1
        do iz=1,nz
        do iy=1,iy1
        do ix=1,nx
          mask(ix,iy,iz)=0.
        enddo
        enddo
        enddo
      endif
      if(iy2.le.ny)then
        if(ldebug)write(*, "('zeroing y ',i3,' to ',i3)")iy2,ny
        do iz=1,nz
        do iy=iy2,ny
        do ix=1,nx
          mask(ix,iy,iz)=0.
        enddo
        enddo
        enddo
      endif
        
      if(iz1.gt.0)then
        if(ldebug)write(*, "('zeroing z ',i3,' to ',i3)")1,iz1
        do iz=1,iz1
        do iy=1,ny
        do ix=1,nx
          mask(ix,iy,iz)=0.
        enddo
        enddo
        enddo
      endif
      if(iz2.le.nz)then
        if(ldebug)write(*, "('zeroing z ',i3,' to ',i3)")iz2,nz
        do iz=iz2,nz
        do iy=1,ny
        do ix=1,nx
          mask(ix,iy,iz)=0.
        enddo
        enddo
        enddo
      endif

      return
      end
