cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
c Washington University, Mallinckrodt Institute of Radiology.
c All Rights Reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of Washington University. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$Id: etadeta.f,v 1.9 2008/08/14 23:54:25 avi Exp $
c$Log: etadeta.f,v $
c Revision 1.9  2008/08/14  23:54:25  avi
c linux compliant
c
c Revision 1.8  2008/07/01  02:58:03  avi
c in resample() check validity of ly(3) only in 3D mode (l3d.eq..true.)
c
c Revision 1.7  2000/12/13  02:54:06  avi
c copyright
c
c Revision 1.6  2000/09/10  04:17:36  avi
c in subroutine msktst protect iz computation from float round errors when nz.eq.1
c
c Revision 1.5  2000/08/28  01:13:04  avi
c fix (high end) array limit checking in subroutine msktst
c
c Revision 1.4  2000/06/04  23:17:32  avi
c consistent range and mask testing with subroutine msktst
c
c Revision 1.3  2000/06/01  03:36:11  avi
c correct sampling grid computation to include mmppix1
c
c Revision 1.2  1998/07/10  07:44:43  avi
c mode independence wrappers
c
c Revision 1.1  1998/07/10  07:24:16  avi
c Initial revision
c
      subroutine eta_rcs
      write(*,"('$Id: etadeta.f,v 1.9 2008/08/14 23:54:25 avi Exp $')")
      return
      end

      subroutine etaderr(t4)
      real*4 t4(4,4)
      write(*,"('etadeta: no voxels in register')")
      write(*,"('t4')")
      do 1 i=1,4
    1 write(*,"(3f10.6,f10.4)")(t4(i,j),j=1,4)
      call exit(-2)
      end

      subroutine etagend(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                   img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,deta,q)
      logical*4 lcorrel
c     mode independence wrapper for etafgd and etauvd
      lcorrel=iand(mode,256).eq.0
      if(lcorrel)then
             call etafgd(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                   img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,deta,q)
      else
             call etauvd(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                   img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,deta,q)
      endif
      return
      end

      subroutine etagen(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                  img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,q)
      logical*4 lcorrel
c     mode independence wrapper for etafg and etauv
      lcorrel=iand(mode,256).eq.0
      if(lcorrel)then
             call etafg(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                  img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,q)
      else
             call etauv(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                  img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,q)
      endif
      return
      end

      logical*4 function msktst(mode,fi,mask,nx,ny,nz)
      real*4 fi(3)		! floating index
      integer*2 mask(0:nx-1,0:ny-1,0:nz-1)
      logical*4 lwrap
      logical*4 ldebug/.false./
c     check index ranges and mask state prior to calling splint subroutines

      msktst=.false.
      lwrap=iand(mode,1024).ne.0
      ix=nint(fi(1)-0.5)
      iy=nint(fi(2)-0.5)
      if(nz.gt.1)then
        iz=nint(fi(3)-0.5)
        if(iz.lt.0.or.iz.gt.nz-1)goto 9
      else
        iz=0
      endif
      if(lwrap)then
        dowhile(ix.lt.0)
          ix=ix+nx
        enddo
        dowhile(iy.lt.0)
          iy=iy+ny
        enddo
        ix1=mod(ix+1,nx)
        iy1=mod(iy+1,ny)
      else
        if(ix.lt.0.or.iy.lt.0)goto 9
        ix1=ix+1
        iy1=iy+1
        if(ix1.ge.nx.or.iy1.ge.ny)goto 9
      endif
      if(mask(ix,iy,iz).le.0.or.mask(ix1,iy,iz).le.0.or.mask(ix,iy1,iz).le.0.or.mask(ix1,iy1,iz).le.0)goto 9
      if(nz.eq.1)return
      iz1=iz+1
      if(iz1.ge.nz)goto 9
      if(mask(ix,iy,iz1).le.0.or.mask(ix1,iy,iz1).le.0.or.mask(ix,iy1,iz1).le.0.or.mask(ix1,iy1,iz1).le.0)goto 9
      return
    9 msktst=.true.
      return
      end

      subroutine etafgd(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                  img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,deta,q)
      real*4 img1(0:nx1-1,0:ny1-1,0:nz1-1),d2x1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 d2y1(0:nx1-1,0:ny1-1,0:nz1-1),d2z1(0:nx1-1,0:ny1-1,0:nz1-1)
      integer*2                            msk1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 mmppix1(3),center1(3)
      real*4 img2(0:nx2-1,0:ny2-1,0:nz2-1),d2x2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 d2y2(0:nx2-1,0:ny2-1,0:nz2-1),d2z2(0:nx2-1,0:ny2-1,0:nz2-1)
      integer*2                            msk2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 mmppix2(3),center2(3)
      real*4 omega(12),deta(12)
c     compute common mode eta and deta/domega
      real*4 rot(3,3),drot(3,3,3),sdrot(3,3,3),t(4,4)
      real*4 x(3),xi(3),y(3),yi(3)
      real*4 g,f,gg,gf,ff,grad(3),ggrad(3),fgrad(3)
      real*4 gdfdd(3)
      real*4 gdfdt(3,3)
      real*4 gdfds(3)
      real*4 gdfd2e(3)
      real*4 fdfdd(3)
      real*4 fdfdt(3,3)
      real*4 fdfds(3)
      real*4 fdfd2e(3)
      logical*4 lrigid,lstretch,lwrap
      logical*4 ldebug/.false./
      external  msktst
      logical*4 msktst

      lwrap=iand(mode,1024).ne.0
      lrigid=iand(mode,4).eq.0
      lstretch=.not.lrigid

      call omega_t4(mode,omega,t)
      if(lrigid)then
        call twoe2rotd(omega(4),rot,drot)
        do i=1,3
        do j=1,3
        do k=1,3
          sdrot(k,j,i)=(1.+omega(6+k))*drot(k,j,i)
        enddo
        enddo
        enddo
      endif

      nvox=0
      f=0.
      g=0.
      ff=0.
      gf=0.
      gg=0.
      do i=1,3
        gdfdd(i)=0.
        fdfdd(i)=0.
        gdfds(i)=0.
        fdfds(i)=0.
        gdfd2e(i)=0.
        fdfd2e(i)=0.
        do j=1,3
          gdfdt(i,j)=0.
          fdfdt(i,j)=0.
        enddo
      enddo
      grad(3)=0.						! in case nz2.eq.1
      do 41 iz=0,int(float(nz1)*abs(mmppix1(3))/del)
      do 41 iy=0,int(float(ny1)*abs(mmppix1(2))/del)
      do 41 ix=0,int(float(nx1)*abs(mmppix1(1))/del)
      x(3)=float(iz)*del-center1(3)				! target coordinate
      x(2)=float(iy)*del-center1(2)				! target coordinate
      x(1)=float(ix)*del-center1(1)				! target coordinate
      do k=1,3
        xi(k)=(x(k)+center1(k))/mmppix1(k)			! target fndex
        y(k)=t(k,1)*x(1)+t(k,2)*x(2)+t(k,3)*x(3)+t(k,4)		! source coordinate
        yi(k)=(y(k)+center2(k))/mmppix2(k)			! source fndex
      enddo
      if(msktst(mode,xi,msk1,nx1,ny1,nz1).or.
     &   msktst(mode,yi,msk2,nx2,ny2,nz2))goto 41
      nvox=nvox+1
      if(nz1.gt.1)then
        call splint3dv (img1,nx1,ny1,nz1,d2x1,d2y1,d2z1,xi(1),xi(2),xi(3),v1)		! target
      else
        call splint2dv (img1,nx1,ny1,d2x1,d2y1,xi(1),xi(2),v1)
      endif
      if(nz2.gt.1)then
        call splint3dvg(img2,nx2,ny2,nz2,d2x2,d2y2,d2z2,yi(1),yi(2),yi(3),v2,grad)	! source
      else
        call splint2dvg(img2,nx2,ny2,d2x2,d2y2,yi(1),yi(2),v2,grad)
      endif
      do k=1,3
        grad(k)=grad(k)/mmppix2(k)
        ggrad(k)=v1*grad(k)
        fgrad(k)=v2*grad(k)
      enddo
      g=g+v1
      f=f+v2
      gg=gg+v1*v1
      gf=gf+v1*v2
      ff=ff+v2*v2
      do i=1,3
        gdfdd(i)=gdfdd(i)+ggrad(i)
        fdfdd(i)=fdfdd(i)+fgrad(i)
        do j=1,3
        if(lstretch)then
          gdfdt(i,j)=gdfdt(i,j)+ggrad(i)*x(j)
          fdfdt(i,j)=fdfdt(i,j)+fgrad(i)*x(j)
        else
          gdfds(i)=gdfds(i)+ggrad(i)*rot(i,j)*x(j)
          fdfds(i)=fdfds(i)+fgrad(i)*rot(i,j)*x(j)
          do k=1,3
            gdfd2e(i)=gdfd2e(i)+ggrad(k)*sdrot(k,j,i)*x(j)
            fdfd2e(i)=fdfd2e(i)+fgrad(k)*sdrot(k,j,i)*x(j)
          enddo
        endif
        enddo
      enddo
   41 continue
      if(nvox.eq.0)call etaderr(t)

      eta=gf/sqrt(gg*ff)
      do i=1,3
        deta(i)=eta*(gdfdd(i)/gf-fdfdd(i)/ff)
      enddo
      if(lstretch)then
        k=4
        do i=1,3
        do j=1,3
          deta(k)=eta*(gdfdt(i,j)/gf-fdfdt(i,j)/ff)
          k=k+1
        enddo
        enddo
      else
        do i=1,3
          deta(3+i)=eta*(gdfd2e(i)/gf-fdfd2e(i)/ff)
          deta(6+i)=eta*(gdfds (i)/gf-fdfds (i)/ff)
        enddo
      endif
      q=f/g
      if(ldebug)write(*,"('etafgd: nvox',i8)")nvox
      return
      end

      subroutine etafg(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                 img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,q)
      real*4 img1(0:nx1-1,0:ny1-1,0:nz1-1),d2x1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 d2y1(0:nx1-1,0:ny1-1,0:nz1-1),d2z1(0:nx1-1,0:ny1-1,0:nz1-1)
      integer*2                            msk1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 mmppix1(3),center1(3)
      real*4 img2(0:nx2-1,0:ny2-1,0:nz2-1),d2x2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 d2y2(0:nx2-1,0:ny2-1,0:nz2-1),d2z2(0:nx2-1,0:ny2-1,0:nz2-1)
      integer*2                            msk2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 mmppix2(3),center2(3)
      real*4 omega(12)
c     compute common mode eta
      real*4 t(4,4)
      real*4 x(3),xi(3),y(3),yi(3)
      real*4 g,f,gg,gf,ff
      logical*4 lwrap
      logical*4 ldebug/.false./
      external  msktst
      logical*4 msktst

      lwrap=iand(mode,1024).ne.0

      call omega_t4(mode,omega,t)

      if(ldebug)then
        write(*,"('etafg: mmppix1,center1',3f10.6,3f10.4)")mmppix1,center1
        write(*,"('etafg: mmppix2,center2',3f10.6,3f10.4)")mmppix2,center2
        write(*,"('etafg: omega',12f8.4)")omega
      endif

      nvox=0
      f=0.
      g=0.
      ff=0.
      gf=0.
      gg=0.
      do 41 iz=0,int(float(nz1)*abs(mmppix1(3))/del)
      do 41 iy=0,int(float(ny1)*abs(mmppix1(2))/del)
      do 41 ix=0,int(float(nx1)*abs(mmppix1(1))/del)
      x(3)=float(iz)*del-center1(3)				! target coordinate
      x(2)=float(iy)*del-center1(2)				! target coordinate
      x(1)=float(ix)*del-center1(1)				! target coordinate
      do k=1,3
        xi(k)=(x(k)+center1(k))/mmppix1(k)			! target fndex
        y(k)=t(k,1)*x(1)+t(k,2)*x(2)+t(k,3)*x(3)+t(k,4)		! source coordinate
        yi(k)=(y(k)+center2(k))/mmppix2(k)			! source fndex
      enddo
      if(msktst(mode,xi,msk1,nx1,ny1,nz1).or.
     &   msktst(mode,yi,msk2,nx2,ny2,nz2))goto 41
      nvox=nvox+1
      if(nz1.gt.1)then
        call splint3dv (img1,nx1,ny1,nz1,d2x1,d2y1,d2z1,xi(1),xi(2),xi(3),v1)		! target
      else
        call splint2dv (img1,nx1,ny1,d2x1,d2y1,xi(1),xi(2),v1)
      endif
      if(nz2.gt.1)then
        call splint3dv(img2,nx2,ny2,nz2,d2x2,d2y2,d2z2,yi(1),yi(2),yi(3),v2)		! source
      else
        call splint2dv(img2,nx2,ny2,d2x2,d2y2,yi(1),yi(2),v2)
      endif
      g=g+v1
      f=f+v2
      gg=gg+v1*v1
      gf=gf+v1*v2
      ff=ff+v2*v2
   41 continue
      if(nvox.eq.0)call etaderr(t)

      eta=gf/sqrt(gg*ff)
      q=f/g
      if(ldebug)write(*,"('etafg: nvox',i8)")nvox
      return
      end

      subroutine etauvd(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                  img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,deta,r)
      real*4 img1(0:nx1-1,0:ny1-1,0:nz1-1),d2x1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 d2y1(0:nx1-1,0:ny1-1,0:nz1-1),d2z1(0:nx1-1,0:ny1-1,0:nz1-1)
      integer*2                            msk1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 mmppix1(3),center1(3)
      real*4 img2(0:nx2-1,0:ny2-1,0:nz2-1),d2x2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 d2y2(0:nx2-1,0:ny2-1,0:nz2-1),d2z2(0:nx2-1,0:ny2-1,0:nz2-1)
      integer*2                            msk2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 mmppix2(3),center2(3)
      real*4 omega(12),deta(12)
c     compute cross modal eta and deta/domega
      real*4 rot(3,3),drot(3,3,3),sdrot(3,3,3),t(4,4)
      real*4 x(3),xi(3),y(3),yi(3)
      real*4 v(3),u(3),q(3),d2fdx2(3,3),qd2fdx2(3),ud2fdx2(3),d2t2(2,2)
      real*4 qdudd(3)
      real*4 qdudt(3,3)
      real*4 qduds(3)
      real*4 qdud2e(3)
      real*4 ududd(3)
      real*4 ududt(3,3)
      real*4 ududs(3)
      real*4 udud2e(3)
      logical*4 lrigid,lstretch,lwrap
      logical*4 ldebug/.false./
      external  msktst
      logical*4 msktst

      lwrap=iand(mode,1024).ne.0
      lrigid=iand(mode,4).eq.0
      lstretch=.not.lrigid

      call omega_t4(mode,omega,t)
      if(lrigid)then
        call twoe2rotd(omega(4),rot,drot)
        do i=1,3
        do j=1,3
        do k=1,3
          sdrot(k,j,i)=(1.+omega(6+k))*drot(k,j,i)
        enddo
        enddo
        enddo
      endif

      nvox=0
      f=0.
      g=0.
      uu=0.
      vv=0.
      uvct=0.
      do i=1,3
        qd2fdx2(i)=0.
        ud2fdx2(i)=0.
        qdudd(i)=0.
        ududd(i)=0.
        qduds(i)=0.
        ududs(i)=0.
        qdud2e(i)=0.
        udud2e(i)=0.
        do j=1,3
          qdudt(i,j)=0.
          ududt(i,j)=0.
          d2fdx2(i,j)=0.					! in case nz2.eq.1
        enddo
      enddo
      u(3)=0.							! in case nz2.eq.1
      v(3)=0.							! in case nz1.eq.1
      do 41 iz=0,int(float(nz1)*abs(mmppix1(3))/del)
      do 41 iy=0,int(float(ny1)*abs(mmppix1(2))/del)
      do 41 ix=0,int(float(nx1)*abs(mmppix1(1))/del)
      x(3)=float(iz)*del-center1(3)				! target coordinate
      x(2)=float(iy)*del-center1(2)				! target coordinate
      x(1)=float(ix)*del-center1(1)				! target coordinate
      do k=1,3
        xi(k)=(x(k)+center1(k))/mmppix1(k)			! target fndex
        y(k)=t(k,1)*x(1)+t(k,2)*x(2)+t(k,3)*x(3)+t(k,4)		! source coordinate
        yi(k)=(y(k)+center2(k))/mmppix2(k)			! source fndex
      enddo
      if(msktst(mode,xi,msk1,nx1,ny1,nz1).or.
     &   msktst(mode,yi,msk2,nx2,ny2,nz2))goto 41
      nvox=nvox+1
      if(nz1.gt.1)then
        call splint3dvg (img1,nx1,ny1,nz1,d2x1,d2y1,d2z1,xi(1),xi(2),xi(3),v1,v)	! target
      else
        call splint2dvg (img1,nx1,ny1,d2x1,d2y1,xi(1),xi(2),v1,v)
      endif
      if(nz2.gt.1)then
        call splint3dvgh(img2,nx2,ny2,nz2,d2x2,d2y2,d2z2,yi(1),yi(2),yi(3),v2,u,d2fdx2)	! source
      else
        call splint2dvgh(img2,nx2,ny2,d2x2,d2y2,yi(1),yi(2),v2,u,d2t2)
        d2fdx2(1,1)=d2t2(1,1)
        d2fdx2(2,1)=d2t2(2,1)
        d2fdx2(1,2)=d2t2(1,2)
        d2fdx2(2,2)=d2t2(2,2)
      endif
      do i=1,3
        v(i)=v(i)/mmppix1(i)
        u(i)=u(i)/mmppix2(i)
        do j=1,3
          d2fdx2(i,j)=d2fdx2(i,j)/(mmppix1(i)*mmppix2(j))
        enddo
      enddo
      aa=u(1)*u(1)+u(2)*u(2)+u(3)*u(3)
      bb=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
      ab=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
      if(aa.eq.0.0.or.bb.eq.0.0)goto 41
      ct=ab/sqrt(aa*bb)
      do k=1,3
        q(k)=ct*(2.*v(k)-(ab/aa)*u(k))
      enddo
      f=f+v2
      g=g+v1
      uu=uu+aa
      vv=vv+bb
      uvct=uvct+ab*ct
      do i=1,3
        qd2fdx2(i)=q(1)*d2fdx2(1,i)+q(2)*d2fdx2(2,i)+q(3)*d2fdx2(3,i)
        ud2fdx2(i)=u(1)*d2fdx2(1,i)+u(2)*d2fdx2(2,i)+u(3)*d2fdx2(3,i)
      enddo
      do i=1,3
        qdudd(i)=qdudd(i)+qd2fdx2(i)
        ududd(i)=ududd(i)+ud2fdx2(i)
        do j=1,3
          if(lstretch)then
            qdudt(i,j)=qdudt(i,j)+qd2fdx2(i)*x(j)
            ududt(i,j)=ududt(i,j)+ud2fdx2(i)*x(j)
          else
            qduds(i)=qduds(i)+qd2fdx2(i)*rot(i,j)*x(j)
            ududs(i)=ududs(i)+ud2fdx2(i)*rot(i,j)*x(j)
            do k=1,3
              qdud2e(i)=qdud2e(i)+qd2fdx2(k)*sdrot(k,j,i)*x(j)
              udud2e(i)=udud2e(i)+ud2fdx2(k)*sdrot(k,j,i)*x(j)
            enddo
          endif
        enddo
      enddo
   41 continue
      if(nvox.eq.0)call etaderr(t)

      eta=uvct/sqrt(uu*vv)
      do i=1,3
        deta(i)=qdudd(i)/sqrt(uu*vv)-eta*ududd(i)/uu
      enddo
      if(lstretch)then
        k=4
        do i=1,3
        do j=1,3
          deta(k)=qdudt(i,j)/sqrt(uu*vv)-eta*ududt(i,j)/uu
          k=k+1
        enddo
        enddo
      else
        do i=1,3
          deta(3+i)=qdud2e(i)/sqrt(uu*vv)-eta*udud2e(i)/uu
          deta(6+i)=qduds (i)/sqrt(uu*vv)-eta*ududs (i)/uu
        enddo
      endif
      r=f/g
      if(ldebug)write(*,"('etauvd: nvox',i8)")nvox
      return
      end

      subroutine etauv(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &                 img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,r)
      real*4 img1(0:nx1-1,0:ny1-1,0:nz1-1),d2x1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 d2y1(0:nx1-1,0:ny1-1,0:nz1-1),d2z1(0:nx1-1,0:ny1-1,0:nz1-1)
      integer*2                            msk1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 mmppix1(3),center1(3)
      real*4 img2(0:nx2-1,0:ny2-1,0:nz2-1),d2x2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 d2y2(0:nx2-1,0:ny2-1,0:nz2-1),d2z2(0:nx2-1,0:ny2-1,0:nz2-1)
      integer*2                            msk2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 mmppix2(3),center2(3)
      real*4 omega(12)
c     compute cross modal eta
      real*4 t(4,4)
      real*4 x(3),xi(3),y(3),yi(3)
      real*4 v(3),u(3)
      logical*4 lrigid,lstretch,lwrap
      logical*4 ldebug/.false./
      external  msktst
      logical*4 msktst

      lwrap=iand(mode,1024).ne.0
      lrigid=iand(mode,4).eq.0
      lstretch=.not.lrigid

      call omega_t4(mode,omega,t)

      nvox=0
      f=0.
      g=0.
      uu=0.
      vv=0.
      uvct=0.
      u(3)=0.							! in case nz2.eq.1
      v(3)=0.							! in case nz1.eq.1
      do 41 iz=0,int(float(nz1)*abs(mmppix1(3))/del)
      do 41 iy=0,int(float(ny1)*abs(mmppix1(2))/del)
      do 41 ix=0,int(float(nx1)*abs(mmppix1(1))/del)
      x(3)=float(iz)*del-center1(3)				! target coordinate
      x(2)=float(iy)*del-center1(2)				! target coordinate
      x(1)=float(ix)*del-center1(1)				! target coordinate
      do k=1,3
        xi(k)=(x(k)+center1(k))/mmppix1(k)			! target fndex
        y(k)=t(k,1)*x(1)+t(k,2)*x(2)+t(k,3)*x(3)+t(k,4)		! source coordinate
        yi(k)=(y(k)+center2(k))/mmppix2(k)			! source fndex
      enddo
      if(msktst(mode,xi,msk1,nx1,ny1,nz1).or.
     &   msktst(mode,yi,msk2,nx2,ny2,nz2))goto 41
      nvox=nvox+1
      if(nz1.gt.1)then
        call splint3dvg(img1,nx1,ny1,nz1,d2x1,d2y1,d2z1,xi(1),xi(2),xi(3),v1,v)	! target
      else
        call splint2dvg(img1,nx1,ny1,d2x1,d2y1,xi(1),xi(2),v1,v)
      endif
      if(nz2.gt.1)then
        call splint3dvg(img2,nx2,ny2,nz2,d2x2,d2y2,d2z2,yi(1),yi(2),yi(3),v2,u)	! source
      else
        call splint2dvg(img2,nx2,ny2,d2x2,d2y2,yi(1),yi(2),v2,u)
      endif
      do i=1,3
        v(i)=v(i)/mmppix1(i)
        u(i)=u(i)/mmppix2(i)
      enddo
      aa=u(1)*u(1)+u(2)*u(2)+u(3)*u(3)
      bb=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
      ab=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
      if(aa.eq.0.0.or.bb.eq.0.0)goto 41
      ct=ab/sqrt(aa*bb)
      f=f+v2
      g=g+v1
      uu=uu+aa
      vv=vv+bb
      uvct=uvct+ab*ct
   41 continue
      if(nvox.eq.0)call etaderr(t)

      eta=uvct/sqrt(uu*vv)
      r=f/g
      if(ldebug)write(*,"('etauv: nvox',i8)")nvox
      return
      end

      subroutine resample(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,t4,
     &                    img2,nx2,ny2,nz2,mmppix2,center2)
      real*4 img1(0:nx1-1,0:ny1-1,0:nz1-1),d2x1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 d2y1(0:nx1-1,0:ny1-1,0:nz1-1),d2z1(0:nx1-1,0:ny1-1,0:nz1-1)
      integer*2                            msk1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 mmppix1(3),center1(3),t4(4,4)
      real*4 img2(0:nx2-1,0:ny2-1,0:nz2-1),mmppix2(3),center2(3)
c     resample img1 to img2 under t4 coordinate transform
      real*4 x(3),y(3),yi(3)
      integer*4 ly(3)
      logical*4 lwrap
      logical*4 l3d
      logical*4 ldebug/.false./
      logical*4 invalid

      lwrap=iand(mode,1024).ne.0
      l3d=iand(mode,2).ne.0

      if(ldebug)then
        write(*,"('resample: mmppix1,center1',3f10.6,3f10.4)")mmppix1,center1
        write(*,"('resample: mmppix2,center2',3f10.6,3f10.4)")mmppix2,center2
      endif

      sum=0.
      do 41 iz=0,nz2-1
      do 41 iy=0,ny2-1
      do 41 ix=0,nx2-1
      x(3)=float(iz)*mmppix2(3)-center2(3)			! target coordinate
      x(2)=float(iy)*mmppix2(2)-center2(2)			! target coordinate
      x(1)=float(ix)*mmppix2(1)-center2(1)			! target coordinate
      do k=1,3
        y(k)=t4(k,1)*x(1)+t4(k,2)*x(2)+t4(k,3)*x(3)+t4(k,4)	! source coordinate
        yi(k)=(y(k)+center1(k))/mmppix1(k)			! source fndex
        ly(k)=nint(yi(k)-0.5)					! source index
      enddo
      if(nz1.eq.1)ly(3)=0
      invalid=.not.lwrap.and.(ly(1).lt.0.or.ly(1).ge.nx1.or.ly(2).lt.0.or.ly(2).ge.ny1)
      if(l3d)invalid=invalid.or.ly(3).lt.0.or.ly(3).ge.nz1
      invalid=invalid.or.msk1(ly(1),ly(2),ly(3)).le.0
c     if(iz.eq.0)then
c       write(*,"('l3d=',l1,' ix iy iz invalid',4i6)")l3d,ix,iy,iz,invalid
c     endif
      if(invalid)then
        img2(ix,iy,iz)=0.
        goto 41
      endif
      nvox=nvox+1
      if(nz1.gt.1)then
        call splint3dv(img1,nx1,ny1,nz1,d2x1,d2y1,d2z1,yi(1),yi(2),yi(3),v)
      else
        call splint2dv(img1,nx1,ny1,d2x1,d2y1,yi(1),yi(2),v)
      endif
      sum=sum+v
      img2(ix,iy,iz)=v
   41 continue

      if(ldebug)write(*,"('resample: nvox',i8,' mean resampled value',f10.4)")nvox,sum/float(nvox)
      return
      end
