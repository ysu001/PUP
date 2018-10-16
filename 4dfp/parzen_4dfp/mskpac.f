c$Header: /data/petsun4/data1/src_solaris/parzen_4dfp/RCS/mskpac.f,v 1.1 2007/07/04 02:00:03 avi Exp $
c$Log: mskpac.f,v $
c Revision 1.1  2007/07/04  02:00:03  avi
c Initial revision
c
      subroutine mskpac(mask,nx,ny,nz,maskp,nxp,nyp,nzp)
      integer*2 mask(nx,ny,nz),maskp(nxp,nyp,nzp)

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
         maskp(ip,jp,kp)=mask(i,j,k)
    1 continue

      return
      end
