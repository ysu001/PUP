c$Header: /data/petsun4/data1/src_solaris/imglin/RCS/stretchout.f,v 1.1 2010/03/03 01:28:53 avi Exp $
c$Log: stretchout.f,v $
c Revision 1.1  2010/03/03  01:28:53  avi
c Initial revision
c
      subroutine stretchout(t)
      real*4 t(4,4)
      real*4 g(3,3),w(3,3),q(3,3),wt(3,3),a(3,3),b(3,3)

      do 1 i=1,3
      do 1 j=1,3
    1 g(i,j)=t(i,j)
      call transpos(g,b,3)
      call matmul(b,g,q,3)
      call eigen(q,w,3)
      call transpos(w,wt,3)
      do 3 i=1,3
    3 q(i,i)=1./sqrt(q(i,i))
      call matmul(w,q,a,3)
      call matmul(a,wt,b,3)
      call matmul(g,b,q,3)
      do 4 i=1,3
      do 4 j=1,3
    4 t(i,j)=q(i,j)

      return
      end
