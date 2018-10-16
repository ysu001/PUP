c$Header: /data/petsun4/data1/src_solaris/librms/RCS/deconv.f,v 1.2 2007/09/18 22:57:39 avi Exp $
c$Log: deconv.f,v $
c Revision 1.2  2007/09/18  22:57:39  avi
c gcc compliant (type -> write(*,)
c
c Revision 1.1  2007/09/18  22:54:47  avi
c Initial revision
c
c     program convtst
      subroutine convtst
      parameter (n=128)
      parameter (scale=40.)
      real*4 h(n),g(n),f(n)
      do 1 i=1,n/2
    1 f(i)=1.
      do 2 i=n/2+1,n
    2 f(i)=-1.
      sum=0.0
      do 3 i=1,n
      t=float(i-1)*5./float(n)
      g(i)=exp(-t)
    3 sum=sum+g(i)
      do 13 i=1,n
   13 g(i)=g(i)/sum
      call conv(f,g,h,n)
      do 4 i=1,n
    4 write(*,101)float(i-1)/scale,f(i),scale*g(i),h(i)
  101 format(4f10.6)
      call deconv(h,g,f,n)
      do 5 i=1,n
    5 write(*,101)float(i-1)/scale,f(i),scale*g(i),h(i)
      stop
      end
      subroutine conv(f,g,h,n)
c     calculate h=f*g given f,g
      real*4 h(n),g(n),f(n)
      call gconv(f,g,h,n,.true.)
      return
      end
      subroutine deconv(h,g,f,n)
c     calculate f where h=f*g given h,g
      real*4 h(n),g(n),f(n)
      call gconv(h,g,f,n,.false.)
      return
      end
      subroutine gconv(h,g,f,n,l)
      real*4 h(n),g(n),f(n),a(2*n+4)
      logical*4 l
      complex*8 ch,cg
      pointer (pa,a)
      pa=malloc(4*(2*n+4))
      if(pa.eq.0) stop "gconv: memory allocation error"
      j=n/2+1
      ihr=1
      ihi=ihr+j
      igr=ihi+j
      igi=igr+j
      j=1
      do 1 i=0,n/2-1
      a(ihr+i)=h(j)
      a(igr+i)=g(j)
      j=j+1
      a(ihi+i)=h(j)
      a(igi+i)=g(j)
      j=j+1
    1 continue
      call fft(a(ihr),a(ihi),1,n/2,1,-1)
      call reals(a(ihr),a(ihi),n/2,-1)
      call fft(a(igr),a(igi),1,n/2,1,-1)
      call reals(a(igr),a(igi),n/2,-1)
      do 2 i=0,n/2
      ch=cmplx(a(ihr+i),a(ihi+i))
      cg=cmplx(a(igr+i),a(igi+i))
      if(l)then
      ch=ch*cg
      else
      ch=ch/cg
      endif
      a(ihr+i)=real(ch)
      a(ihi+i)=aimag(ch)
    2 continue
      call reals(a(ihr),a(ihi),n/2,+1)
      call fft(a(ihr),a(ihi),1,n/2,1,+1)
      j=1
      do 3 i=0,n/2-1
      f(j)=a(ihr+i)
      j=j+1
      f(j)=a(ihi+i)
      j=j+1
    3 continue
      call free(pa)
      return
      end
