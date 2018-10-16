c$Header: /data/petsun4/data1/src_solaris/interp_4dfp/RCS/spline.f,v 1.3 2007/11/20 07:06:36 avi Exp $
c$Log: spline.f,v $
c Revision 1.3  2007/11/20  07:06:36  avi
c gcc v4 compliant (type -> write)
c
c Revision 1.2  2004/05/17  22:33:06  avi
c correct bug in splintv() for out of bounds evaluation
c
c Revision 1.1  2002/04/04  03:53:01  avi
c Initial revision
c
      subroutine spline_test
      real*4 f(0:18)/19*0./
      real*4 d2f(0:18)

      do 1 i=8,12
    1 f(i)=1.0
      call spline(f,19,d2f)

      do i=-10,190
        z=0.1*float(i)
        write(*,"(3f10.6)")z,splintd(f,19,d2f,z,df),df
      enddo

      call exit(0)
      end

      subroutine spline(data,n,d2da)
      real*4 data(0:n-1),d2da(0:n-1)
      real*4 a(0:0),b(0:0),f(0:0)
      pointer (pa,a),(pb,b),(pf,f)

      pa=malloc(4*n)
      pb=malloc(4*n)
      pf=malloc(4*((n+1)/2+1))
      twopi=8.*atan(1.)
      do i=0,(n+1)/2
        c=cos(twopi*float(i)/float(n))
        f(i)=6.*(c-1.)/(c+2.)
      enddo

      do i=0,n-1
        a(i)=data(i)
        b(i)=0.
      enddo
      call fft(a,b,1,n,1,+1)
      a(0)=0.
      do i=1,(n+1)/2-1
        a(i)=a(i)*f(i)
        b(i)=b(i)*f(i)
        a(n-i)=a(n-i)*f(i)
        b(n-i)=b(n-i)*f(i)
      enddo
      if(mod(n,2).eq.0)a(n/2)=a(n/2)*(-12.)
      call fft(a,b,1,n,1,-1)
      do i=0,n-1
        d2da(i)=a(i)
      enddo

      call free(pa)
      call free(pb)
      call free(pf)
      return
      end

      function splint(data,n,d2da,z)
      real*4 data(0:n-1),d2da(0:n-1)

      i=nint(z-0.5)
      if(i.lt.0)then
        splint=data(0)
        return
      endif

      if(i.gt.n-2)then
        splint=data(n-1)
        return
      endif

      w=z-float(i)
      a=1.-w
      splint=a*data(i)+w*data(i+1)+((a**3-a)*d2da(i)+(w**3-w)*d2da(i+1))/6.

      return
      end

      subroutine splintv(data,n,d2da,z,v)
      real*4 data(0:n-1),d2da(0:n-1)

      i=nint(z-0.5)
      if(i.lt.0)then
        v=data(0)
        return
      endif

      if(i.gt.n-2)then
        v=data(n-1)
        return
      endif

      w=z-float(i)
      a=1.-w
      v=a*data(i)+w*data(i+1)+((a**3-a)*d2da(i)+(w**3-w)*d2da(i+1))/6.

      return
      end

      function splintd(data,n,d2da,z,dz)
      real*4 data(0:n-1),d2da(0:n-1)

      i=nint(z-0.5)
      w=z-float(i)
      dowhile(i.lt.0)
        i=i+n
      enddo
      i=mod(i,n)
      ii=mod(i+1,n)

      a=1.-w
      c=a*(a*a-1.)/6.
      d=w*(w*w-1.)/6.
      e=(1.-3.*a*a)/6.
      f=(3.*w*w-1.)/6.
      splintd= a*data(i)+w*data(ii)+c*d2da(i)+d*d2da(ii)
      dz=       -data(i)  +data(ii)+e*d2da(i)+f*d2da(ii)

      return
      end
