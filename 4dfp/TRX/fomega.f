cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008
c Washington University, Mallinckrodt Institute of Radiology.
c All Rights Reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of Washington University. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$Header: /data/petsun4/data1/src_solaris/TRX/RCS/fomega.f,v 1.5 2008/08/15 03:02:31 avi Exp $
c$Log: fomega.f,v $
c Revision 1.5  2008/08/15  03:02:31  avi
c correct memory allocation error in sr_factor() due to not allocation array rs(3,3)
c
c Revision 1.4  2007/07/18  22:55:46  avi
c typo (inot->not)
c
c Revision 1.3  2007/07/18  22:44:30  avi
c gcc compliant
c
c Revision 1.2  2000/12/13  02:57:03  avi
c copyright
c
c Revision 1.1  2000/05/22  19:27:22  avi
c Initial revision
c
      subroutine t4_omega_test
      real*4 t4(4,4),omega(12),t(4,4)
      real*4 twoe(3),s(3,3),r(3,3),sr(3,3)
      integer*4 mode/4/
      logical*4 lstretch

      lstretch=iand(mode,4).ne.0

      do 9 m=1,20
      call t4_init(t4)
      if(m.eq.1)goto 7
      do 1 i=1,3
    1 t4(i,4)=20.*rand(0)-10.
      if(lstretch)then
        do i=1,3
          do j=1,3
            t4(i,j)=t4(i,j)+2.*rand(0)-1.
          enddo
        enddo
      else
        do i=1,3
          do j=1,3
            s(i,j)=0.
          enddo
          s(i,i)=rand(0)+0.5
          twoe(i)=2.*rand(0)-1.
        enddo
        call twoe2rot(twoe,r)
        call matmul(s,r,sr,3)
        if(rand(0).gt.0.5)then
          do j=1,3
            f      =sr(2,j)
            sr(2,j)=sr(3,j)
            sr(3,j)=f
          enddo
        endif
        do i=1,3
          do j=1,3
            t4(i,j)=sr(i,j)
          enddo
        enddo
      endif
    7 call t4_omega(mode,t4,omega)
      write(*,"(3f10.4,3f10.6/6f10.6)")omega
      call omega_t4(mode,omega,t)
      do 8 i=1,4
    8 write(*,"(3(3f10.6,f10.4,10x))")(t4(i,j),j=1,4),(t(i,j)-t4(i,j),j=1,4)
    9 write(*,"()")

      call exit(0)
      end

      subroutine t4_omega(mode,t4,omega)
      real*4 t4(4,4),omega(12)
      real*4 sr(3,3),s(3,3),r(3,3)
      logical*4 lstretch,lreflect
      logical*4 ldebug/.true./

      lstretch=iand(mode,4).ne.0

      do 1 i=1,3
      do 2 j=1,3
    2 sr(i,j)=t4(i,j)
    1 omega(i)=t4(i,4)

      det=   sr(1,1)*(sr(2,2)*sr(3,3)-sr(3,2)*sr(2,3))
     &      -sr(1,2)*(sr(2,1)*sr(3,3)-sr(3,1)*sr(2,3))
     &      +sr(1,3)*(sr(2,1)*sr(3,2)-sr(3,1)*sr(2,2))
      write(*,"('t4_omega: determinant',f10.6)")det
      lreflect=det.lt.0.
      mode=iand(mode,not('40000'x))		! clear negative determinant bit
      if(lreflect)then
        do j=1,3
          f      =sr(2,j)
          sr(2,j)=sr(3,j)
          sr(3,j)=f
        enddo
        mode=ior(mode,'40000'x)			! set negative determinant bit
      endif

      if(lstretch)then				! stretch mode
        i=4
        do l=1,3
          do m=1,3
            omega(i)=sr(l,m)
            if(m.eq.l)omega(i)=omega(i)-1.
            i=i+1
          enddo
        enddo
      else					! rigid body mode
        call sr_factor(sr,s,r)
        call rot2twoe(r,omega(4))
        mskbit=16
        do i=1,3
          if(iand(mode,mskbit).ne.0)then
            omega(6+i)=s(i,i)-1.
          else
            omega(6+i)=0.
          endif
          mskbit=mskbit*2
          omega(9+i)=0.
        enddo
      endif

      return
      end

      subroutine omega_t4(mode,omega,t4)
      real*4 omega(12),t4(4,4)
      real*4 sr(3,3)
      logical*4 lstretch,lreflect
      logical*4 ldebug/.false./

      lstretch=iand(mode,4).ne.0
      lreflect=iand(mode,'40000'x).ne.0
      call t4_init(t4)

      do 1 i=1,3				! translation
    1 t4(i,4)=omega(i)

      if(lstretch)then				! stretch mode
        i=4
        do l=1,3
          do m=1,3
            t4(l,m)=t4(l,m)+omega(i)
            i=i+1
          enddo
        enddo
      else					! rigid body mode
        call twoe2rot(omega(4),sr)
        do i=1,3
          do j=1,3
            t4(i,j)=sr(i,j)
          enddo
        enddo
        do i=1,3
          do j=1,3
            t4(i,j)=t4(i,j)*(1.+omega(6+i))
          enddo
        enddo
      endif

      if(lreflect)then				! exchange rows 2 and 3
        if(ldebug)write(*,"('exchanging')")
        do j=1,3
          f      =t4(2,j)
          t4(2,j)=t4(3,j)
          t4(3,j)=f
        enddo
      endif

      return
      end

      subroutine t4_init(t4)
      real*4 t4(4,4)
      do i=1,4
        do j=1,4
          t4(i,j)=0.
        enddo
        t4(i,i)=1.
      enddo
      return
      end

      subroutine sr_factor(sr,s,r)
      real*4 sr(3,3),s(3,3),r(3,3),rs(3,3)
      real*4 ss(3,3),sinv(3,3),w(3,3),wt(3,3),rt(3,3),q(3,3),a(3,3)
      logical*4 ldebug/.false./

      do 1 i=1,3
      do 1 j=1,3
    1 q(i,j)=0.

      call transpos(sr,rs,3)
      call matmul(sr,rs,ss,3)
      call eigen(ss,w,3)
      call transpos(w,wt,3)
      do 3 i=1,3
    3 q(i,i)=1./sqrt(ss(i,i))
      call matmul(w,q,a,3)
      call matmul(a,wt,sinv,3)
      call matmul(sinv,sr,r,3)
      call transpos(r,rt,3)
      call matmul(r,rt,q,3)
      if(ldebug)then
        write(*,"('r,i')")
        do i=1,3
          write(*,"(2(3f10.6,2x))")(r(i,j),j=1,3),(q(i,j),j=1,3)
        enddo
      endif
      call matmul(sr,rt,s,3)
      call matmul(s,r,q,3)
      if(ldebug)then
        write(*,"('s,sr,0')")
        do i=1,3
          write(*,"(3(3f10.6,2x))")(s(i,j),j=1,3),(sr(i,j),j=1,3),(sr(i,j)-q(i,j),j=1,3)
        enddo
      endif
      return
      end

