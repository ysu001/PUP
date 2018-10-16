cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011
c Washington University, Mallinckrodt Institute of Radiology.
c All rights reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of washington university. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Header: /data/petsun4/data1/src_solaris/librms/RCS/dmatinv.f,v 1.6 2011/01/19 04:54:43 avi Exp $
c $Log: dmatinv.f,v $
c Revision 1.6  2011/01/19  04:54:43  avi
c restore malloc functionality
c
c Revision 1.5  2010/02/12  07:14:43  avi
c hard-coded matrix dimenssion limit 512->2048
c
c Revision 1.4  2007/05/05  20:17:03  avi
c typo
c
c Revision 1.3  2007/05/05  20:14:08  avi
c Revision 1.2  2007/05/05  20:05:19  avi
c gcc v3 compliant (no pointers, drand()->rand)
c
c Revision 1.1  2005/09/04  06:24:51  avi
c Initial revision
c
c     program dmatinvtst
      subroutine dmatinvtst
      implicit real*8 (a-h,o-z)
      real*4 rand
      parameter (n=384)
      real*8 q(n,n),t(n,n),r(n,n)
      do 1 i=1,n
      do 1 j=1,n
    1 t(i,j)=0.5*(rand(0)-.1)
      do 2 i=1,n
      do 2 j=1,n
      q(i,j)=0.0
      do 2 k=1,n
      q(i,j)=q(i,j)+t(i,k)*t(j,k)
    2 r(i,j)=q(i,j)
      call dmatinv(q,n,det)
      write(*,"('dimension=',i4,' det=',e12.6)")n,det
  101 format(3(3f16.12,4x))
      do 12 i=1,n
      do 12 j=1,n
      t(i,j)=0.0
      do 12 k=1,n
   12 t(i,j)=t(i,j)+q(i,k)*r(k,j)
      write(*,"()")
      if (n.eq.3)then
        do i=1,n
          write(*,101)(r(i,k),k=1,n),(q(i,k),k=1,n),(t(i,k),k=1,n)
         enddo
      endif
      derr=0.0
      ierr=0
      jerr=0
      do 6 i=1,n
      do 6 j=1,n
      x=t(i,j)
      if(i.eq.j)x=x-1.0
      if(dabs(x).gt.derr)then
        ierr=i
        jerr=j
        derr=dabs(x)
      endif
    6 continue
      write(*,"('maximum error = ',e12.6)")derr
      end
      subroutine dmatinv(a,n,det)
c     algorithm follows Bevington, Phillip R.,
c     Data Reduction and Error Analysis for the Physical Sciences,
c     Mcgraw-Hill, New York, 1969.
      implicit real*8 (a-h,o-z)
      real*8 a(n,n),det
      pointer (pik,ik),(pjk,jk)
      integer*4 ik(n),jk(n)
c     integer*4 ik(2048),jk(2048)
c     if(n.gt.2048)then
c       write (*,"('dmatinv matrix dimension limit (',i4,') exeeded')")2048
c       call exit (-1)
c     endif
      pik=malloc(4*n)
      pjk=malloc(4*n)
      det=1.
      do 100 k=1,n
      vmax=0.
   21 do 30 i=k,n
      do 30 j=k,n
      if(dabs(vmax).gt.dabs(a(i,j)))goto 30
      vmax=a(i,j)
      ik(k)=i
      jk(k)=j
   30 continue
      i=ik(k)
c     if(i-k)21,51,43
      if(i.eq.k)goto 51
      if(i.lt.k)goto 21
   43 do 50 j=1,n
      t=a(k,j)
      a(k,j)=a(i,j)
   50 a(i,j)=-t
   51 j=jk(k)
c     if(j-k)21,61,53
      if(j.eq.k)goto 61
      if(j.lt.k)goto 21
   53 do 60 i=1,n
      t=a(i,k)
      a(i,k)=a(i,j)
   60 a(i,j)=-t
   61 do 70 i=1,n
      if(i.ne.k)a(i,k)=-a(i,k)/vmax
   70 continue
   71 do 80 i=1,n
      do 80 j=1,n
      if(i.ne.k.and.j.ne.k)a(i,j)=a(i,j)+a(i,k)*a(k,j)
   80 continue
   81 do 90 j=1,n
      if(j.ne.k)a(k,j)=a(k,j)/vmax
   90 continue
      a(k,k)=1./vmax
  100 det=det*vmax
      do 130 l=1,n
      k=n-l+1
      j=ik(k)
c     if(j-k)111,111,105
      if(j.le.k)goto 111
  105 do 110 i=1,n
      t=a(i,k)
      a(i,k)=-a(i,j)
  110 a(i,j)=t
  111 i=jk(k)
c     if(i-k)130,130,113
      if(i.le.k)goto 130
  113 do 120 j=1,n
      t=a(k,j)
      a(k,j)=-a(i,j)
  120 a(i,j)=t
  130 continue
      call free(pik)
      call free(pjk)
      return
      end
