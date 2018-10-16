c$Header: /data/petsun4/data1/src_solaris/parzen_4dfp/RCS/imgpac.f,v 1.2 2007/07/04 02:22:19 avi Exp $
c$Log: imgpac.f,v $
c Revision 1.2  2007/07/04  02:22:19  avi
c type -> write
c
c Revision 1.1  2007/07/04  01:54:05  avi
c Initial revision
c
      subroutine imgpac(imag,nx,ny,nz,imagp,nxp,nyp,nzp)
      real*4 imag(nx,ny,nz),imagp(nxp,nyp,nzp)
      logical*4 ldebug/.false./

      mx=(nxp-nx)/2
      my=(nyp-ny)/2
      mz=(nzp-nz)/2
      do 1 kp=1,nzp
         k=kp-mz
         if(k.lt.1)then
            k=1
         elseif(k.gt.nz)then
            k=nz
         endif
      if(ldebug)write (*, "('kp,k,f ',2i6,f10.6)")kp,k,f
      do 1 jp=1,nyp
         j=jp-my
         if(j.lt.1)then
            j=1
         elseif(j.gt.ny)then
            j=ny
         endif
      do 1 ip=1,nxp
         i=ip-mx
         if(i.lt.1)then
            i=1
         elseif(i.gt.nx)then
            i=nx
         endif
         imagp(ip,jp,kp)=imag(i,j,k)
    1 continue

      return
      end
