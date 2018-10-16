c$Header: /data/petsun4/data1/src_solaris/imglin/RCS/t4inv.f,v 1.1 2007/05/01 01:18:57 avi Exp $
c$Log: t4inv.f,v $
c Revision 1.1  2007/05/01  01:18:57  avi
c Initial revision
c
      subroutine t4inv(t,tinv)
c     extracted from param12opr.f
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
