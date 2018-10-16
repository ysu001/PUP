c$Header: /data/petsun4/data1/src_solaris/librms/RCS/determ12.f,v 1.1 2007/06/19 22:34:50 avi Exp $
c$Log: determ12.f,v $
c Revision 1.1  2007/06/19  22:34:50  avi
c Initial revision
c
      real*4 function determ12(a,det)
      real*4 a(4,4)

      det=+a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     &    -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
     &    +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      determ12 = det
      return
      end
