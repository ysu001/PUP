c$Header: /data/petsun4/data1/src_solaris/imglin/RCS/param12opr.f,v 1.2 2009/02/23 06:14:16 avi Exp $
c$Log: param12opr.f,v $
c Revision 1.2  2009/02/23  06:14:16  avi
c gcc compliant
c
c     param12opr.f
c     Avi Snyder 13-Nov-93
c                31-Oct-94 ustretchset and uparam subroutines

C     PARAM(01:03) NUDGE
C     PARAM(04:06) ANGLE
C     PARAM(07:09) CENTER
C     PARAM(10:12) MMPPIX
C     PARAM(13:15) STRETCH
C     PARAM(16:18) STRANG

c     uparam(01:03) nudge
c     uparam(04:06) angle
c     uparam(07:09) center
c     uparam(10:12) mmppix
c     uparam(13:18) ustretch

      subroutine param18tst
      real*4 t(4,4),t1(4,4)
      real*4 param0(18),param1(18),param2(18)

      do 9 n=1,30

      do 2 i=1,3
    2 param0(i)=5.0*(2.*rand(0)-1.)			! nudge
      do 3 i=4,6
    3 param0(i)=0.2*(2.*rand(0)-1.)			! angle

      if(.true.)then
        do i=13,18
          param0(i)=0.1*(2.*rand(0)-1.)			! ustretch
        enddo
        write(*, "(6f12.6)")(param0(i),i=1,6),(param0(i),i=13,18)
        call uparam2warp(param0(01),param0(04),param0(13),t)
        write(*, "('test warp matrix via uparam2warp')")
        call t4lst(t)
        call warp2uparam(t,param2(01),param2(04),param2(13))
        write(*, "(6f12.6)")(param2(i),i=1,6),(param2(i),i=13,18)
        call uparam2warp(param2(01),param2(04),param2(13),t1)
        write(*, "('test warp matrix after passage through warp2uparam and back')")
        call t4lst(t1)
        do i=1,4
          do j=1,4
            t(i,j)=t(i,j)-t1(i,j)
          enddo
        enddo
        write(*, "('difference')")
        call t4lst(t)
      endif

      if(.false.)then
        do i=13,15
          param0(i)=0.1*(2.*rand(0)-1.)			! stretch
        enddo
        do i=16,18
          param0(i)=1.0*(2.*rand(0)-1.)			! strang
        enddo
        write(*, "(6f12.6)")(param0(i),i=1,6),(param0(i),i=13,18)
        call param2warp(param0(01),param0(04),param0(13),param0(16),t)
        call t4lst(t)
        call warp2param(t,param2(01),param2(04),param2(13),param2(16))
        write(*, "(6f12.6)")(param2(i),i=1,6),(param2(i),i=13,18)
        call param2warp(param2(01),param2(04),param2(13),param2(16),t1)
        call t4lst(t1)
        call warp2param(t1,param2(01),param2(04),param2(13),param2(16))
        write(*, "(6f12.6)")(param2(i),i=1,6),(param2(i),i=13,18)
        do 5 i=1,4
        do 5 j=1,4
    5   t(i,j)=t(i,j)-t1(i,j)
        call t4lst(t)
      endif
      if(.false.)then
        write(*, "(6f12.6)")(param0(i),i=1,6),(param0(i),i=13,18)
        call param18inv(param0,param1)
        write(*, "(6f12.6)")(param1(i),i=1,6),(param1(i),i=13,18)
        call param18mul(param0,param1,param2)
        write(*, "(6f12.6)")(param2(i),i=1,6),(param2(i),i=13,18)
      endif
    9 continue
      return
      end

      subroutine t4lst(a)
      real*4 a(4,4)
      do 1 i=1,4
    1 write(*, "(4f10.6)")(a(i,j),j=1,4)
      return
      end

      subroutine param2warp(nudge,angle,stretch,strang,t)
      real*4 nudge(3),angle(3),stretch(3),strang(3),t(4,4)
      real*4 a(4,4),b(4,4)
      real*4 z(3)/3*0./
      call trotset(z,angle,a)
      call stretchset(stretch,strang,b)
      call matmul(b,a,t,4)
      do 1 i=1,3
    1 t(i,4)=nudge(i)
      return
      end

      subroutine param18inv(param1,param2)
      real*4 param1(18),param2(18)
      real*4 t(4,4),tinv(4,4)
      call param2warp(param1(01),param1(04),param1(13),param1(16),t)
      call t4inv(t,tinv)
      call warp2param(tinv,param2(01),param2(04),param2(13),param2(16))
      return
      end

      subroutine param18mul(param1,param2,param3)
      real*4 param1(18),param2(18),param3(18)
      real*4 a(4,4),b(4,4),c(4,4)
      call param2warp(param1(01),param1(04),param1(13),param1(16),a)
      call param2warp(param2(01),param2(04),param2(13),param2(16),b)
      call matmul(a,b,c,4)
      call warp2param(c,param3(01),param3(04),param3(13),param3(16))
      return
      end

      subroutine t4inv(t,tinv)
      real*4 t(4,4),tinv(4,4)
      real*4 sr(3,3),d(3),g(3,3),q(3,3)
      do 1 i=1,3
      d(i)=t(i,4)
      do 1 j=1,3
    1 sr(i,j)=t(i,j)
      call geninv(sr,g,q,3,det)
      do 2 i=1,3
      do 2 j=1,3
    2 tinv(i,j)=q(i,j)
      do 3 i=1,3
      tinv(4,i)=0.
      tinv(i,4)=0.
      do 3 k=1,3
    3 tinv(i,4)=tinv(i,4)-tinv(i,k)*d(k)
      tinv(4,4)=1.
      return
      end

      subroutine warp2param(t,nudge,angle,stretch,strang)
      real*4 nudge(3),angle(3),stretch(3),strang(3),t(4,4)
      real*4 sr(3,3),r(3,3),sinv(3,3),w(3,3),q(3,3),wt(3,3)
      do 1 i=1,3
      nudge(i)=t(i,4)
      do 1 j=1,3
    1 sr(i,j)=t(i,j)
      do 2 i=1,3
      do 2 j=1,3
      q(i,j)=0.
      do 2 k=1,3
    2 q(i,j)=q(i,j)+sr(i,k)*sr(j,k)
      call eigen(q,w,3)
c     write(*, "('warp2parm stretch eigenvalues',3f10.6)")(q(i,i),i=1,3)
      do 3 i=1,3
      stretch(i)=0.5*alog(q(i,i))
    3 q(i,i)=1./sqrt(q(i,i))
      call transpos(w,wt,3)
      call rot2ang(wt,strang)
      call matmul(w,q,r,3)
      call matmul(r,wt,sinv,3)
      call matmul(sinv,sr,r,3)
      call rot2ang(r,angle)
      return
      end

      subroutine bparam18(param)
c     beautifies param by minimizing stretch angles
      real*4 param(19)
      real*4 t(4,4),q(3,3),w(3,3),wt(3,3)
      real*4 strang(3)
      real*4 sum(6)
      integer*4 iperm(3,6)/1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/
      call param2warp(param(01),param(04),param(13),param(16),t)
      do 2 i=1,3
      do 2 j=1,3
      q(i,j)=0.
      do 2 k=1,3
    2 q(i,j)=q(i,j)+t(i,k)*t(j,k)
      call eigen(q,w,3)
      do 9 m=1,6
      do 8 i=1,3
      ii=iperm(i,m)
      do 8 j=1,3
    8 wt(i,j)=w(j,ii)
      call rot2ang(wt,strang)
    9 sum(m)=abs(strang(1))+abs(strang(2))+abs(strang(3))
      mm=1
      do 11 m=2,6
   11 if(sum(m).lt.sum(mm))mm=m
      do 12 i=1,3
      ii=iperm(i,mm)
      param(12+i)=0.5*alog(q(ii,ii))
      do 12 j=1,3
   12 wt(i,j)=w(j,ii)
      call rot2ang(wt,param(16))
      return
      end

c     uparam(01:03) nudge
c     uparam(04:06) angle
c     uparam(07:09) center
c     uparam(10:12) mmppix
c     uparam(13:18) ustretch

      subroutine ustretchset(ustretch,t)
      real*4 ustretch(6),t(4,4)
      do 1 i=1,4
      do 2 j=1,4
    2 t(i,j)=0.
    1 t(i,i)=1.
      k=1
      do 3 i=1,3
      do 3 j=i,3
      t(i,j)=t(i,j)+ustretch(k)
      t(j,i)=t(i,j)
    3 k=k+1
      return
      end

      subroutine uparam2warp(nudge,angle,ustretch,t)
      real*4 nudge(3),angle(3),ustretch(6),t(4,4)
      real*4 a(4,4),b(4,4)
      real*4 z(3)/3*0./
      call trotset(z,angle,a)
      call ustretchset(ustretch,b)
      call matmul(b,a,t,4)
      do 1 i=1,3
    1 t(i,4)=nudge(i)
      return
      end

      subroutine warp2uparam(t,nudge,angle,ustretch)
      real*4 nudge(3),angle(3),ustretch(6),t(4,4)
      real*4 r(3,3),s(3,3),w(3,3),q(3,3),wt(3,3),b(3,3)
      do 1 i=1,3
    1 nudge(i)=t(i,4)
      do 2 i=1,3
      do 2 j=1,3
      q(i,j)=0.
      do 2 k=1,3
    2 q(i,j)=q(i,j)+t(i,k)*t(j,k)
      call eigen(q,w,3)
      call transpos(w,wt,3)
      do 4 i=1,3
      do 3 j=1,3
    3 b(i,j)=0.
    4 b(i,i)=sqrt(q(i,i))
      call matmul(w,b,q,3)
      call matmul(q,wt,s,3)
      k=1
      do 5 i=1,3
      do 5 j=i,3
      if(j.eq.i)then
        ustretch(k)=s(i,j)-1.
      else
        ustretch(k)=s(i,j)
      endif
    5 k=k+1
      do 6 i=1,3
      do 6 j=1,3
    6 q(i,j)=t(i,j)
      call matinv(s,3,det)
      call matmul(s,q,r,3)
      call rot2ang(r,angle)
      return
      end
