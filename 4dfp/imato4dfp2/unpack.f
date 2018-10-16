c$Header: /data/petsun4/data1/src_solaris/imato4dfp2/RCS/unpack.f,v 1.2 2003/02/16 05:14:18 avi Exp $
c$Log: unpack.f,v $
c Revision 1.2  2003/02/16  05:14:18  avi
c subroutine unpackf
c
c Revision 1.1  1998/09/16  22:50:43  avi
c Initial revision
c
      subroutine unpack_rcs
      write(*,"('$Id: unpack.f,v 1.2 2003/02/16 05:14:18 avi Exp $')")
      return
      end

      subroutine unpackf(imgm,nxp,nyp,imgu,nx,ny,np)
      real*4 imgm(0:nxp-1,0:nyp-1),imgu(0:nx-1,0:ny-1,0:np-1)

      if((nxp*nyp).ne.(nx*ny*np))then
        write(*,"('unpackf: dimension error')")
        call exit(-2)
      endif

      my=nyp/ny
      mx=nxp/nx
      ip=0
      do 1 j=0,my-1
      do 1 i=0,mx-1
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
    2 imgu(ix,iy,ip)=imgm(i*nx+ix,j*ny+iy)
    1 ip=ip+1

      return
      end

      subroutine unpack(imgm,nxp,nyp,imgu,nx,ny,np)
      integer*2 imgm(0:nxp-1,0:nyp-1),imgu(0:nx-1,0:ny-1,0:np-1)

      if((nxp*nyp).ne.(nx*ny*np))then
        write(*,"('unpack: dimension error')")
        call exit(-2)
      endif

      my=nyp/ny
      mx=nxp/nx
      ip=0
      do 1 j=0,my-1
      do 1 i=0,mx-1
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
    2 imgu(ix,iy,ip)=imgm(i*nx+ix,j*ny+iy)
    1 ip=ip+1

      return
      end

      subroutine resort(imgu,ns,nz,imgr)
      integer*2 imgu(ns,nz)
      real*4    imgr(ns,nz)
      logical*4 ldebug/.false./

      ihalf=(nz+1)/2
      m=mod(nz,2)

      jz=1
      do 1 iz=1,nz-1,2
      if(ldebug)write(*,"('iz=',i2,' jz=',i2)")nz-iz+1,jz
      do 2 is=1,ns
    2 imgr(is,nz-iz+1)=float(imgu(is,jz))
      if(ldebug)write(*,"('iz=',i2,' jz=',i2)")nz-iz,jz+ihalf
      do 3 is=1,ns
    3 imgr(is,nz-iz+0)=float(imgu(is,jz+ihalf))
    1 jz=jz+1
      if(m.eq.0)return
      if(ldebug)write(*,"('iz=',i2,' jz=',i2)")1,jz
      do 4 is=1,ns
    4 imgr(is,1)=float(imgu(is,jz))

      return
      end
