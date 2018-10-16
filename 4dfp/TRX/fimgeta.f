ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007			c
c Washington University, Mallinckrodt Institute of Radiology.				c
c All Rights Reserved.									c
c This software may not be reproduced, copied, or distributed without written		c
c permission of Washington University. For further information contact A. Z. Snyder.	c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$Header: /data/petsun4/data1/src_solaris/TRX/RCS/fimgeta.f,v 1.19 2007/07/18 22:15:05 avi Exp $
c$Log: fimgeta.f,v $
c Revision 1.19  2007/07/18  22:15:05  avi
c gcc compliant
c
c Revision 1.18  2005/12/29  03:14:21  avi
c update copyright
c
c Revision 1.17  2001/01/04  05:36:52  avi
c mode 2048 (parameter search only)
c
c Revision 1.16  2000/12/13  02:53:11  avi
c copyright
c
c Revision 1.15  2000/06/05  03:19:56  avi
c precede discrete Newton method with parameter search
c
c Revision 1.14  2000/06/05  01:49:04  avi
c rscale = 0.1
c improved test mode
c
c Revision 1.13  1998/07/29  04:59:22  avi
c write to screen only free elements of omega
c
c Revision 1.12  1998/07/29  02:41:53  avi
c make termination more dependent on stepacc
c
c Revision 1.11  1998/07/28  04:20:32  avi
c Revision 1.10  1998/07/28  04:12:18  avi
c terminate by disc
c remove if(lunstable) code from subroutine fimgeta
c
c Revision 1.9  1998/07/28  02:50:04  avi
c hessinv prints after each iteration
c
c Revision 1.8  1998/07/28  02:34:09  avi
c correct absent eigv(144) declaration in subroutine fimgeta
c force positive definite hessian uwing knowledge of least eigenvalue
c eliminate call to save_omega
c
c Revision 1.7  1998/07/27  20:51:50  avi
c basically works but Hessian becomes indeterminate too easily
c
c Revision 1.6  1998/07/27  03:46:17  avi
c correct hessinv
c do not adjust eps using eigenvalues
c
c Revision 1.5  1998/07/26  00:19:52  avi
c Revision 1.4  1998/07/25  22:36:21  avi
c after using etagen
c
c Revision 1.3  1998/07/25  22:21:28  avi
c before using etagen
c
c Revision 1.2  1998/06/22  04:56:06  avi
c image regsitration subroutines
c
      subroutine fimgeta_rcs
      character*256 rcsid/'$Id: fimgeta.f,v 1.19 2007/07/18 22:15:05 avi Exp $'/
      write(*,"(a)")rcsid
      return
      end

      subroutine mode_lfree(mode,lfree,nterm)
      integer*4 mode,lfree(12),nterm
      logical*4 lenable,l2d,lrigid,lstretch,l3d,ltrans_only
      logical*4 ldebug/.false./

      lenable=iand(mode,1).ne.0
      l2d=iand(mode,2).eq.0
      lrigid=iand(mode,4).eq.0
      ltrans_only=iand(mode,8).ne.0
      lstretch=.not.lrigid
      l3d=.not.l2d
      
      nterm=0
      do i=1,12
        lfree(i)=0
      enddo
      if(.not.lenable)return
      if(l3d.and.lstretch)then
        do i=1,12
          lfree(i)=1
        enddo
      else if(l3d.and.lrigid)then
        do i=1,6
          lfree(i)=1
        enddo
        mskbit=16
        do i=7,9
          if(iand(mode,mskbit).ne.0)lfree(i)=1
          mskbit=mskbit*2
        enddo
      else if(l2d.and.lstretch)then
        do j=0,2
        do i=1,2
          lfree(i+3*j)=1
        enddo
        enddo
      else if(l2d.and.lrigid)then
        lfree(1)=1
        lfree(2)=1
        lfree(6)=1
        if(iand(mode,16).ne.0)lfree(7)=1
        if(iand(mode,32).ne.0)lfree(8)=1
      endif
      if(ltrans_only)then
        do i=4,6
          lfree(i)=0
        enddo
      endif
      do 1 i=1,12
    1 if(lfree(i).ne.0)nterm=nterm+1

      if(ldebug)write(*,"('lfree ',12i2)")lfree
      return
      end

      subroutine fimgeta(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,del,
     &                   img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,t4)
      real*4 img1(0:nx1-1,0:ny1-1,0:nz1-1),d2x1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 d2y1(0:nx1-1,0:ny1-1,0:nz1-1),d2z1(0:nx1-1,0:ny1-1,0:nz1-1)
      integer*2                            msk1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 mmppix1(3),center1(3)
      real*4 img2(0:nx2-1,0:ny2-1,0:nz2-1),d2x2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 d2y2(0:nx2-1,0:ny2-1,0:nz2-1),d2z2(0:nx2-1,0:ny2-1,0:nz2-1)
      integer*2                            msk2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 mmppix2(3),center2(3)
      real*4 t4(4,4)
c     find optimal translation/rotation to achieve registration of img2 on img1
      real*4 deta(12),deta1(12),deta2(12)
      real*4 omega(12),omega1(12),omega2(12),domega(12),omegam(12)
      real*4 eps(12)/3*3.0,9*0.075/	! 0.075 ~ 3 / (characteristic_radius=40)
      real*4 etat(-64:63)
      real*4 heta(144),eigv(144),eigw(144),hinv(144)
      logical*4 lfree(12)
      logical*4 lenable,l3d,lstretch,lcorrel,ltest,lsearch
      logical*4 lunstable
      logical*4 ldebug/.false./
      parameter (kfac=2)
c     static eps,eta,omega,omegam

      if(ldebug)then
        write (*, "(6i6)")nx1,ny1,nz1,nx2,ny2,nz2
        write (*, "(6f10.4)")mmppix1,center1
        write (*, "(6f10.4)")mmppix2,center2
      endif
      lenable=iand(mode,1).ne.0
      l3d=iand(mode,2).ne.0
      lstretch=iand(mode,4).ne.0
      lcorrel=iand(mode,256).eq.0
      ltest=iand(mode,128).ne.0
      lsearch=iand(mode,2048).ne.0
      call mode_lfree(mode,lfree,nterm)
      write(*,"('image alignment mode ',i6,' decimal ',z8,' hex')")mode,mode
      escale=-0.1*float(nterm)
      if(lcorrel)escale=escale*25.
      rscale=0.1

      call t4_omega(mode,t4,omega)
      if(ldebug)then
        do i=1,4
          write(*,"(3f10.6,f10.4)")(t4(i,k),k=1,4)
        enddo
      endif
      if(ltest)goto 60

      call etagen(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &            img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,q)
      write(*,"('omega'/3f8.2,9f8.4)")(omega(k),k=1,nterm)
      write(*,"('eta,q   ',f8.6,f8.4)")eta,q
      if(nterm.eq.0)return

      niter=1
   40 lunstable=.false.
      write(*,"('fimgeta: parameter search iteration',i3)")niter
      do 41 j=1,12
      if(.not.lfree(j))goto 41
      call etagen(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &            img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta0,q)
      do 42 i=1,12
      omega1(i)=omega(i)
   42 omega2(i)=omega(i)
      omega1(j)=omega(j)-eps(j)
      omega2(j)=omega(j)+eps(j)
      call etagen(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega1,del,
     &            img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta1,q1)
      call etagen(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega2,del,
     &            img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta2,q2)
      if(eta1.gt.eta0.and.eta1.gt.eta2)then
        omega(j)=omega1(j)
        write(*,"('omega'/3f8.2,9f8.4)")(omega(k),k=1,nterm)
        write(*,"('eta,q   ',f8.6,f8.4)")eta1,q1
        lunstable=.true.
      endif
      if(eta2.gt.eta0.and.eta2.gt.eta1)then
        omega(j)=omega2(j)
        write(*,"('omega'/3f8.2,9f8.4)")(omega(k),k=1,nterm)
        write(*,"('eta,q   ',f8.6,f8.4)")eta2,q2
        lunstable=.true.
      endif
   41 continue
      niter=niter+1
      if(niter.gt.32)then
        write(*,"('fimgeta: iteration limit reached')")
        goto 90
      endif
      if(lunstable)goto 40
      if(lsearch)goto 90

      stepacc=0.
      nhess=0
   50 nhess=nhess+1
      nunstable=0
      write(*,"('fimgeta: discrete Newton method iteration',i3)")nhess
      if(nhess.gt.32)then
        write(*,"('fimgeta: iteration limit reached')")
        goto 90
      endif
      call etagend(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &             img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta0,deta,q)
      jj=0
      do 51 j=1,12
      if(.not.lfree(j))goto 51
      do 52 i=1,12
      omega1(i)=omega(i)
   52 omega2(i)=omega(i)
      omega1(j)=omega(j)-eps(j)*rscale
      omega2(j)=omega(j)+eps(j)*rscale
      call etagend(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega1,del,
     &             img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta1,deta1,q1)
      call etagend(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega2,del,
     &             img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta2,deta2,q2)
      ii=0
      do 57 i=1,12
      if(.not.lfree(i))goto 57
      ii=ii+1
      heta(ii+nterm*jj)=0.5*eps(i)*(deta2(i)-deta1(i))*escale/rscale
   57 continue
      jj=jj+1
   51 continue
      call hessinv(nterm,heta,eigv,eigw,hinv,lunstable)
      step=0.
      ii=0
      do 55 i=1,12
      domega(i)=0.
      if(.not.lfree(i))goto 55
      ii=ii+1
      jj=0
      do j=1,12
        if(lfree(j))then
          domega(i)=domega(i)-eps(i)*hinv(ii+nterm*jj)*deta(j)*eps(j)*escale
          jj=jj+1
        endif
      enddo
      step=step+(domega(i)/eps(i))**2
   55 continue
      step=sqrt(step/float(nterm))
      if(step.gt.2.)then
        do i=1,12
          domega(i)=domega(i)/step
        enddo
        step=1.
      endif

      do k=-64,63
        etat(k)=0.
      enddo
      k=kfac/2
      loop=0
   58 do 59 l=k-1,k+1
      loop=loop+1
      if(loop.gt.30)goto 90
      if(etat(l).ne.0.)goto 59
      do i=1,12
        omega1(i)=omega(i)+float(l)*domega(i)/float(kfac)
      enddo      
      call etagen(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega1,del,
     &            img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,etat(l),q)
   59 continue
      disc=etat(k+1)-2.*etat(k)+etat(k-1)
      if(ldebug)write(*,"(i10,4f10.6)")k,etat(k-1),etat(k),etat(k+1),disc
      if(abs(disc).lt.1.e-6)goto 90
      if(etat(k).gt.etat(k-1).and.etat(k).gt.etat(k+1))then
        a=float(k)-0.5*(etat(k+1)-etat(k-1))/disc
        goto 56
      endif
      if(etat(k-1).gt.etat(k).and.etat(k-1).gt.etat(k+1))then
        k=k-1
        if(k.ge.-2*kfac)goto 58
      else
        k=k+1
        if(k.le.+2*kfac)goto 58
      endif
      a=float(k)
   56 do i=1,12
        omega(i)=omega(i)+a*domega(i)/float(kfac)
      enddo      

      if(stepacc.eq.0.)stepacc=step*abs(a)/float(kfac)
      stepacc=    0.5*(stepacc+step*abs(a)/float(kfac))
      call etagen(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &            img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,q)
      write(*,"('omega post step ',2f8.4/3f8.2,9f8.4)")step*a/float(kfac),stepacc,(omega(k),k=1,nterm)
      write(*,"('eta,eta0,q      ',3f8.4)")eta,eta0,q
      if(stepacc.gt.0.01)goto 50

   90 call omega_t4(mode,omega,t4)
   70 call etagen(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega,del,
     &            img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta,q)
      write(*,"('omega'/3f8.2,9f8.4)")(omega(k),k=1,nterm)
      write(*,"('eta,q   ',f8.6,f8.4)")eta,q
      return

   60 write(*,"('fimgeta range check')")
      do 61 j=1,12
      if(.not.lfree(j))goto 61
      do 62 i=1,12
   62 omega1(i)=omega(i)
      write(*,"('   epsfract      domega         eta deta/domega      direct    analytic')")
      do 63 k=-6,5
      t=float(k)*0.1
      omega1(j)=omega(j)+t*eps(j)
      call etagend(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,omega1,del,
     &             img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,eta1,deta1,q1)
      q=(eta1-eta)/(0.1*eps(j))
      if(k.gt.-6)write(*,"(i3,f8.4,7f12.6)")j,t,omega1(j),eta1,deta1(j),q,0.5*(deta(j)+deta1(j))
      eta=eta1
   63 deta(j)=deta1(j)
   61 continue
      return

      end

      subroutine hessinv_test
      real*4 hess(144),eigv(144),eigw(144),hinv(144)

      do 9 n=1,20
      nterm=int(12.*rand(0))+1
      write (*, "(i6)")nterm
      do 8 i=1,nterm**2
    8 hess(i)=rand(0)-0.5
      call hessinv(nterm,hess,eigv,eigw,hinv,lunstable)
      call matmul(hess,hinv,eigv,nterm)
    9 write(*,"(12f8.4)")(eigv(i+nterm*(i-1)),i=1,nterm)

      call exit(0)
      end

      subroutine hessinv(nterm,hess,eigv,eigw,hinv,lunstable)
      real*4 hess(nterm,nterm),eigv(nterm,nterm),hinv(nterm,nterm),eigw(nterm,nterm)
      logical*4 lunstable
      real*4 htmp(144),eigt(144)
      logical*4 ldebug/.false./

      small=0.
    5 if(ldebug)then
        write(*,"('previous least eigenvalue',f10.4/'hessian')")small
        do i=1,nterm
          write(*,"(12f8.4)")(hess(i,j),j=1,nterm)
        enddo
      endif
      do 1 i=1,nterm
      do 1 j=1,nterm
    1 eigv(i,j)=0.5*(hess(i,j)+hess(j,i))
      call matcop(eigv,hess,nterm)
      call eigen(eigv,eigw,nterm)
      call transpos(eigw,eigt,nterm)
      big=eigv(1,1)
      if(.not.(big.gt.0.))then
        write(*,"('hessinv: numerical disaster')")
        call exit(-1)
      endif
      small=eigv(nterm,nterm)
      cndnum=big/small
      lunstable=cndnum.lt.0.
      det=1.
      do 7 i=1,nterm
    7 det=det*eigv(i,i)
      write(*,"('hessian condition number',f8.2,' det',e12.4)")cndnum,det
      if(lunstable)then
        do i=1,nterm
          hess(i,i)=hess(i,i)-2.*small
        enddo
        goto 5
      endif
      do 6 i=1,nterm
    6 eigv(i,i)=1./eigv(i,i)
      call matmul(eigw,eigv,htmp,nterm)
      call matmul(htmp,eigt,hinv,nterm)

      return
      end
