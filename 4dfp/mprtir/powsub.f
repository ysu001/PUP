c$Header: /data/petsun4/data1/src_solaris/mprtir/RCS/powsub.f,v 1.5 2010/11/29 06:35:32 avi Exp $
c$Log: powsub.f,v $
c Revision 1.5  2010/11/29  06:35:32  avi
c correct termination of powstr
c
c Revision 1.4  2008/01/02  01:24:01  avi
c gcc v4 compliant
c
c Revision 1.3  1999/12/08  22:59:51  avi
c right justify powstring
c
c Revision 1.2  1999/02/18  01:24:48  avi
c improve speed
c
c Revision 1.1  1999/02/15  10:07:28  avi
c Initial revision
      subroutine evalpow_test
      parameter(n=2)
      parameter(nterm=((n+3)*(n+2)*(n+1))/6)
      real*4 a(nterm),b(nterm),t(nterm,nterm)
      pointer (pa,a),(pb,b),(pt,t)
      character*256 string
      data x,y,z/1.,2.,3./
      
      pa=malloc(4*nterm)
      pb=malloc(4*nterm)
      pt=malloc(4*nterm**2)

      write(*,"('nterm=',i5)")nterm
      if(8*nterm.lt.255)then
        call powstring(n,8,string)
        write(*,"(a)")string(1:8*nterm)
      endif
      call evalpo1(x,y,z,n,a)
      call evalpow(x,y,z,n,b)
      do i=1,nterm
        write(*,"(i10,2f10.4)")i,a(i),b(i)
      enddo

      call free(pa)
      call free(pb)
      call free(pt)
      call exit(0)
      end

      subroutine evalpow(x,y,z,norder,p)
      real*4 p(0:((norder+3)*(norder+2)*(norder+1))/6-1)

      p(0)=1.
      if(norder.lt.1)return
      p(1)=z
      p(2)=y
      p(3)=x
      if(norder.lt.2)return
      p(4)=z*z
      p(5)=y*z
      p(6)=y*y
      p(7)=x*z
      p(8)=x*y
      p(9)=x*x
      if(norder.lt.3)return
      j=10
      do 7 i=3,norder
      fx=1.
      do 4 k=0,norder
      fy=fx
      do 5 l=0,norder
      fz=fy
      do 6 m=0,norder
      if(k+l+m.eq.i)then
        p(j)=fz
        j=j+1
      endif
    6 fz=fz*z
    5 fy=fy*y
    4 fx=fx*x
    7 continue

      return
      end

      subroutine poly3ev(x,y,z,norder,a,p,g)
      real*4 a(0:((norder+3)*(norder+2)*(norder+1))/6-1)
      real*4 p(0:((norder+3)*(norder+2)*(norder+1))/6-1)

      p(0)=1.
      g=a(0)
      if(norder.lt.1)return
      p(1)=z
      p(2)=y
      p(3)=x
      g=g+z*a(1)+y*a(2)+x*a(3)
      if(norder.lt.2)return
      p(4)=z*z
      p(5)=y*z
      p(6)=y*y
      p(7)=x*z
      p(8)=x*y
      p(9)=x*x
      g=g+a(4)*p(4)+a(5)*p(5)+a(6)*p(6)+a(7)*p(7)+a(8)*p(8)+a(9)*p(9)
      if(norder.lt.3)return
      j=10
      do 7 i=3,norder
      fx=1.
      do 4 k=0,norder
      fy=fx
      do 5 l=0,norder
      fz=fy
      do 6 m=0,norder
      if(k+l+m.eq.i)then
        p(j)=fz
        g=g+a(j)*fz
        j=j+1
      endif
    6 fz=fz*z
    5 fy=fy*y
    4 fx=fx*x
    7 continue

      return
      end

      subroutine evalpo1(x,y,z,norder,a)
      real*4 a(((norder+3)*(norder+2)*(norder+1))/6)

      j=0
      do 7 i=0,norder
      do 7 k=0,norder
      do 7 l=0,norder
      do 7 m=0,norder
      if(k+l+m.ne.i)goto 7
      j=j+1
      a(j)=x**k*y**l*z**m
    7 continue

      return
      end

      subroutine powstring(norder,lenfield,string)
      character*255 string
      logical*4 ldebug/.false./

      nterm=((norder+3)*(norder+2)*(norder+1))/6
      string=' '
      j=lenfield-5
      do 7 i=0,norder
      do 7 k=0,norder
      do 7 l=0,norder
      do 7 m=0,norder
      if(k+l+m.ne.i)goto 7
      write(string(j:j+lenfield-1),"('x',i1,'y',i1,'z',i1)")k,l,m
      j=j+lenfield
    7 continue
      if(ldebug)write(*,"(a)")string(1:lenfield*nterm)

      string(lenfield*nterm+1:)=char(0)
      return
      end

