c$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/dglm_4dfp.f,v 1.3 2008/11/26 02:06:51 avi Exp $
c$Log: dglm_4dfp.f,v $
c Revision 1.3  2008/11/26  02:06:51  avi
c move cndnum test to main c program
c
c Revision 1.2  2008/03/17  18:09:05  avi
c trap ill-conditioned design matrix
c
c Revision 1.1  2005/09/05  00:52:46  avi
c Initial revision
c
c Revision 1.1  2004/05/26  05:31:05  avi
c Initial revision
c
      subroutine df2finvt(f,npts,ncol,a,finvt,nnez)
      implicit real*8 (a-h,o-z)
      integer*4 npts,ncol
      real*8 f(npts,ncol),finvt(npts,ncol),a(ncol,ncol)
      real*8 ainv(ncol,ncol),e(ncol,ncol),w(ncol,ncol)
      pointer (painv,ainv)
      pointer (pe,e)
      pointer (pw,w)

c     write(*,"('f')")
c     call matlst(f,npts,ncol)
c     write(*,"('a')")
c     call matlst(a,ncol,ncol)

      painv=malloc(8*ncol*ncol)
      pe   =malloc(8*ncol*ncol)
      pw   =malloc(8*ncol*ncol)
      if(painv.eq.0.or.pe.eq.0.or.pw.eq.0) stop 'df2finvt memory allocation error'

      do 23 i=1,ncol
      do 23 j=1,ncol
   23 ainv(i,j)=a(i,j)
      call dmatinv(ainv,ncol,det)
      write(*,"('det=',e12.6)")det

      do 22 i=1,npts
      do 22 j=1,ncol
      finvt(i,j)=0.0
      do 22 k=1,ncol
   22 finvt(i,j)=finvt(i,j)+f(i,k)*ainv(k,j)

      do 24 i=1,ncol
      do 24 j=1,ncol
      e(i,j)=0.0
      do 25 k=1,npts
   25 e(i,j)=e(i,j)+finvt(k,i)*f(k,j)
   24 e(i,j)=e(i,j)/float(nnez)
c     write(*,"('identity matrix')")
c     call matlst(e,ncol,ncol)
      call errlist(e,ncol)

      call free(painv)
      call free(pe)
      call free(pw)
      return
      end

      subroutine errlist(t,n)
      implicit real*8 (a-h,o-z)
      real*8 t(n,n)
      derr=0.0
      do 6 i=1,n
      do 6 j=1,n
      x=t(i,j)
      if(i.eq.j)x=x-1.0
      if(dabs(x).gt.derr)derr=dabs(x)
    6 continue
      write(*,"('maximum_error=',e12.6)")derr
      end

      subroutine matlst(a,ni,nj)
      real*8 a(ni,nj)
      do 1 i=1,ni
    1 write(*,"(10f10.6)")(a(i,j),j=1,nj)
      return
      end
