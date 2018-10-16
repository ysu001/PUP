c$Header: /data/petsun4/data1/src_solaris/TRX/RCS/fimgetae.f,v 1.2 2008/08/15 03:08:50 avi Exp $
c$Log: fimgetae.f,v $
c Revision 1.2  2008/08/15  03:08:50  avi
c linux compliant
c
c Revision 1.1  2005/12/29  07:32:21  avi
c Initial revision
c
c     subroutine essentially identical to fimgeta() except that
c     eps(12) is passed as an argument
      subroutine fimgetae(img1,d2x1,d2y1,d2z1,nx1,ny1,nz1,mmppix1,center1,msk1,mode,del,eps,
     &                    img2,d2x2,d2y2,d2z2,nx2,ny2,nz2,mmppix2,center2,msk2,t4)
      real*4 img1(0:nx1-1,0:ny1-1,0:nz1-1),d2x1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 d2y1(0:nx1-1,0:ny1-1,0:nz1-1),d2z1(0:nx1-1,0:ny1-1,0:nz1-1)
      integer*2                            msk1(0:nx1-1,0:ny1-1,0:nz1-1)
      real*4 mmppix1(3),center1(3)
      real*4 eps(12)		!/3*3.0,9*0.075/	!0.075 ~3/(characteristic_radius=40)
      real*4 img2(0:nx2-1,0:ny2-1,0:nz2-1),d2x2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 d2y2(0:nx2-1,0:ny2-1,0:nz2-1),d2z2(0:nx2-1,0:ny2-1,0:nz2-1)
      integer*2                            msk2(0:nx2-1,0:ny2-1,0:nz2-1)
      real*4 mmppix2(3),center2(3)
      real*4 t4(4,4)
c     find optimal translation/rotation to achieve registration of img2 on img1
      real*4 deta(12),deta1(12),deta2(12)
      real*4 omega(12),omega1(12),omega2(12),domega(12),omegam(12)
      real*4 etat(-64:63)
      real*4 heta(144),eigv(144),eigw(144),hinv(144)
      logical*4 lfree(12)
      logical*4 lenable,l3d,lstretch,lcorrel,ltest,lsearch
      logical*4 lunstable
      logical*4 ldebug/.false./
      parameter (kfac=2)
c     static eta,omega,omegam

      if(ldebug)then
        write(*, "(6i6)")nx1,ny1,nz1,nx2,ny2,nz2
        write(*, "(6f10.4)")mmppix1,center1
        write(*, "(6f10.4)")mmppix2,center2
        write(*, "(6f10.4)")eps
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
      write(*,"('fimgetae: parameter search iteration',i3)")niter
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
        write(*,"('fimgetae: iteration limit reached')")
        goto 90
      endif
      if(lunstable)goto 40
      if(lsearch)goto 90

      stepacc=0.
      nhess=0
   50 nhess=nhess+1
      nunstable=0
      write(*,"('fimgetae: discrete Newton method iteration',i3)")nhess
      if(nhess.gt.32)then
        write(*,"('fimgetae: iteration limit reached')")
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

   60 write(*,"('fimgetae range check')")
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
