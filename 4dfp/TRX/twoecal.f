cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
c Washington University, Mallinckrodt Institute of Radiology.
c All Rights Reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of Washington University. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$Header: /data/petsun4/data1/src_solaris/TRX/RCS/twoecal.f,v 1.4 2007/07/18 22:32:26 avi Exp $
c$Log: twoecal.f,v $
c Revision 1.4  2007/07/18  22:32:26  avi
c gcc compliant
c
c Revision 1.3  2000/12/13  02:58:04  avi
c copyright
c
c Revision 1.2  1999/09/20  03:50:25  avi
c subroutines twoecal_testd1 and twoecal_testd2
c
c Revision 1.1  1999/09/20  01:37:38  avi
c Initial revision
c
      subroutine twoecal_rcs
      write(*,"('$Id: twoecal.f,v 1.4 2007/07/18 22:32:26 avi Exp $')")
      return
      end

      subroutine twoecal_testl

      do 1 k=1,3
      do 1 i=1,3
      do 1 j=1,3
    1 write (*, "(3i1,7x,i10)")i,j,k,levi_civita(i,j,k)

      call exit(0)
      end

      subroutine twoecal_testd2
      real*4 rot(3,3),twoe(3),twof(3),drot(3,3,3)
      real*4 a(3,3),rotinv(3,3),drotinv(3,3,3)

      do k=1,10
        do i=1,3
           twoe(i)=2.*rand(0)-1.
           twof(i)=-twoe(i)
        enddo
        write(*,"('k=',i2,'  twoe',3f10.6)")k,twoe
        call twoe2rotd(twoe,rot,drot)
        call twoe2rotd(twof,rotinv,drotinv)
	call matmul(rot,rotinv,a,3)
        write(*,"('[R][Rt]')")
        call rotlst(a)
        do i=1,3
          write(*,"('dR/d2e',i1)")i
          call rotlst(drot(1,1,i))
          call transpos(drotinv(1,1,i),a,3)
          write(*,"('-dRt/d2e',i1)")i
          call rotlst(a)
          write(*,"()")
        enddo
      enddo
      call exit(0)
      end

      subroutine twoecal_testd1
      real*4 rot(3,3),twoe(3),drot(3,3,3)
      real*4 a(3,3),b(3,3),rotinv(3,3)

      do k=1,10
        do i=1,3
           twoe(i)=2.*rand(0)-1.
        enddo
        write(*,"('k=',i2,'  twoe',3f10.6)")k,twoe
        call twoe2rotd(twoe,rot,drot)
        do i=1,3
          write(*,"('dRt/d2e',i1)")i
          call transpos(drot(1,1,i),a,3)
          call rotlst(a)
          call transpos(rot,rotinv,3)
          call matmul(rotinv,drot(1,1,i),a,3)
          call matmul(a,rotinv,b,3)
          write(*,"('[Rt][dR/d2e',i1,'][Rt]')")i
          call rotlst(b)
          write(*,"()")
        enddo
      enddo
      call exit(0)
      end

      subroutine twoecal_testd
      real*4 rot(3,3),twoe(3),drot(3,3,3)
      real*4 a(3,3),b(3,3),c(3,3)
      data del/0.01/

      do k=1,10
        do i=1,3
           twoe(i)=2.*rand(0)-1.
        enddo
        write(*,"('k=',i2,'  twoe',3f10.6)")k,twoe
        call twoe2rotd(twoe,rot,drot)
        do i=1,3
          twoe(i)=twoe(i)+del
          call twoe2rot(twoe,b)
          twoe(i)=twoe(i)-2.*del
          call twoe2rot(twoe,a)
          twoe(i)=twoe(i)+del
          call rotdif(del,a,b,c)
          call rotlst(c)
          call rotlst(drot(1,1,i))
          write(*,"()")
        enddo
      enddo
      call exit(0)
      end

      subroutine rotdif(del,a,b,c)
      real*4 a(3,3),b(3,3),c(3,3)
      do 1 i=1,3
      do 1 j=1,3
    1 c(i,j)=(b(i,j)-a(i,j))/(2.*del)
      return
      end

      subroutine twoecal_test1
      real*4 rot(3,3),twoe(3)
c     real*4 a(3,3),b(3,3)

      do k=1,10
        do i=1,3
           twoe(i)=2.*rand(0)-1.
        enddo
        write(*,"('twoe      ',3f10.6)")twoe
        call twoe2rot(twoe,rot)
        call rotlst(rot)
c       call transpos(rot,a,3)
c       call matmul(rot,a,b,3)
c       call rotlst(b)
        call rot2twoe(rot,twoe)
        write(*,"('twoe      ',3f10.6)")twoe
        write(*,"()")
      enddo
      call exit(0)
      end

      subroutine rot2twoe(rot,twoe)
      real*4 rot(3,3),twoe(3)
      logical*4 ldebug/.false./

      twoe0=sqrt(rot(1,1)+rot(2,2)+rot(3,3)+1.)
      twoe(1)=(rot(2,3)-rot(3,2))/twoe0
      twoe(2)=(rot(3,1)-rot(1,3))/twoe0
      twoe(3)=(rot(1,2)-rot(2,1))/twoe0
      return
      end

      subroutine twoe2rot(twoe,rot)
      real*4 rot(3,3),twoe(3)

      twoe0=sqrt(4.-twoe(1)**2-twoe(2)**2-twoe(3)**2)
      rot(1,1)=1.-0.5*(twoe(2)**2+twoe(3)**2)
      rot(2,2)=1.-0.5*(twoe(1)**2+twoe(3)**2)
      rot(3,3)=1.-0.5*(twoe(1)**2+twoe(2)**2)
      rot(1,2)=0.5*(twoe(1)*twoe(2)+twoe0*twoe(3))
      rot(2,1)=0.5*(twoe(1)*twoe(2)-twoe0*twoe(3))
      rot(1,3)=0.5*(twoe(1)*twoe(3)-twoe0*twoe(2))
      rot(3,1)=0.5*(twoe(1)*twoe(3)+twoe0*twoe(2))
      rot(2,3)=0.5*(twoe(2)*twoe(3)+twoe0*twoe(1))
      rot(3,2)=0.5*(twoe(2)*twoe(3)-twoe0*twoe(1))
      return
      end

      subroutine twoe2rotd(twoe,rot,drot)
      real*4 rot(3,3),twoe(3)
      real*4 drot(3,3,3)
      logical*4 ldebug/.false./

      call twoe2rot(twoe,rot)
      twoe0inv=1./sqrt(4.-twoe(1)**2-twoe(2)**2-twoe(3)**2)
      
      if(ldebug)then
      drot(1,1,1)=0.
      drot(1,2,1)=-twoe0inv*rot(1,3)
      drot(1,3,1)= twoe0inv*rot(1,2)
      drot(2,1,1)= twoe0inv*rot(3,1)
      drot(2,2,1)=-twoe(1)
      drot(2,3,1)= twoe0inv*(rot(2,2)+rot(3,3))
      drot(3,1,1)=-twoe0inv*rot(2,1)
      drot(3,2,1)=-drot(2,3,1)
      drot(3,3,1)=-twoe(1)

      drot(1,1,2)=-twoe(2)
      drot(1,2,2)=-twoe0inv*rot(3,2)
      drot(1,3,2)=-twoe0inv*(rot(1,1)+rot(3,3))
      drot(2,1,2)= twoe0inv*rot(2,3)
      drot(2,2,2)=0.
      drot(2,3,2)=-twoe0inv*rot(2,1)
      drot(3,1,2)=-drot(1,3,2)
      drot(3,2,2)= twoe0inv*rot(1,2)
      drot(3,3,2)=-twoe(2)

      drot(1,1,3)=-twoe(3)
      drot(1,2,3)= twoe0inv*(rot(1,1)+rot(2,2))
      drot(1,3,3)= twoe0inv*rot(2,3)
      drot(2,1,3)=-drot(1,2,3)
      drot(2,2,3)=-twoe(3)
      drot(2,3,3)=-twoe0inv*rot(1,3)
      drot(3,1,3)=-twoe0inv*rot(3,2)
      drot(3,2,3)= twoe0inv*rot(3,1)
      drot(3,3,3)=0.

      write (*, "('in-line')")
      do k=1,3
        call rotlst(drot(1,1,k))
      enddo
      endif

      do 1 k=1,3
      do 1 i=1,3
      do 1 j=1,3
      m=6/(i*j)
      emji=float(levi_civita(m,j,i))
      eijk=float(levi_civita(i,j,k))
      if(i.eq.k.and.j.eq.k)then
        drot(i,j,k)=0.
      elseif(i.eq.j)then
        drot(i,j,k)=-twoe(k)
      elseif(i.eq.k)then
        drot(i,j,k)=emji*twoe0inv*rot(k,m)
      elseif(j.eq.k)then
        drot(i,j,k)=emji*twoe0inv*rot(m,k)
      else
        drot(i,j,k)=eijk*twoe0inv*(rot(i,i)+rot(j,j))
      endif
    1 continue

      if(ldebug)then
      write (*, "('do-loop')")
      do k=1,3
        call rotlst(drot(1,1,k))
      enddo
      endif

      return
      end

      integer*4 function levi_civita(i,j,k)
      integer*4 n(3)/2,3,1/
c     static n

      if(j.eq.n(i).and.k.eq.n(j).and.i.eq.n(k))then
        levi_civita= 1
      elseif(j.eq.n(k).and.k.eq.n(i).and.i.eq.n(j))then
        levi_civita=-1
      else
        levi_civita= 0
      endif

      return
      end

      subroutine rotlst(rot)
      real*4 rot(3,3)
      det=+rot(1,1)*(rot(2,2)*rot(3,3)-rot(2,3)*rot(3,2))
     &    -rot(1,2)*(rot(2,1)*rot(3,3)-rot(2,3)*rot(3,1))
     &    +rot(1,3)*(rot(2,1)*rot(3,2)-rot(2,2)*rot(3,1))
      write(*,"('det=      ',f10.6)")det
      do 1 i=1,3
    1 write(*,"(3f10.6)")(rot(i,j),j=1,3)
      return
      end
