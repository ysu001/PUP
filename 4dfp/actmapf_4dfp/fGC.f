c$Header: /data/petsun4/data1/src_solaris/actmapf_4dfp/RCS/fGC.f,v 1.8 2008/09/14 02:45:54 avi Exp $
c$Log: fGC.f,v $
c Revision 1.8  2008/09/14  02:45:54  avi
c correct format statement in debug code
c
c Revision 1.7  2006/08/08  03:40:06  avi
c make compute_gc1() return Fx.y and Fx,y correctly
c
c Revision 1.6  2006/08/06  03:35:03  avi
c more general compute_gc1()
c
c Revision 1.5  2006/08/05  02:41:44  avi
c subroutine compute_gc1()
c
c Revision 1.4  2006/07/30  04:54:56  avi
c include computation of sigma in solve_var()
c
c Revision 1.3  2006/07/30  03:25:29  avi
c rewrite in line with corrected AR equations
c
c Revision 1.2  2006/07/27  01:29:46  avi
c subroutine compute_gc
c
c Revision 1.1  2006/07/26  02:24:20  avi
c Initial revision
c
      subroutine solve_var(ncol,iorder,g,cc,ca,cb,a,det,sigma)
      real*4 g(ncol,ncol,0:iorder)
      real*4 cc(ncol*iorder,ncol*iorder),ca(ncol*iorder,ncol),cb(ncol*iorder,ncol)
      real*4 a(ncol,ncol,iorder),det,sigma(ncol,ncol)
      logical*4 ldebug/.false./

      do 1 i = 0,iorder-1
      do 1 j = i,iorder-1
      do 1 ii = 1,ncol
      do 1 jj = 1,ncol
      cc(ncol*i+ii,ncol*j+jj)=g(ii,jj,j-i)
    1 cc(ncol*j+jj,ncol*i+ii)=g(ii,jj,j-i)
      if(ldebug)then
      write(*,"('cc')")
      call matlist(cc,ncol*iorder,ncol*iorder)
      endif

      do 2 i = 1,iorder
      do 2 ii = 1,ncol
      do 2 jj = 1,ncol
    2 cb(ncol*(i-1)+ii,jj)=-g(jj,ii,i)
      if(ldebug)then
      write(*,"('cb')")
      call matlist(cb,ncol*iorder,ncol)
      endif

      call matinv(cc,ncol*iorder,det)
      if(ldebug)then
      write(*,"('cc inverse det=',e12.4)")det
      call matlist(cc,ncol*iorder,ncol*iorder)
      endif

      do 11 j = 1,ncol
      do 11 i = 1,ncol*iorder
      ca(i,j)=0.0
      do 11 k = 1,ncol*iorder
   11 ca(i,j)=ca(i,j)+cc(i,k)*cb(k,j)
      if(ldebug)then
      write(*,"('ca')")
      call matlist(ca,ncol*iorder,ncol)
      endif

      do 12 i = 1,iorder
      do 12 ii = 1,ncol
      do 12 jj = 1,ncol
   12 a(jj,ii,i)=ca(ncol*(i-1)+ii,jj)

      call matcop(g,sigma,ncol)
      do 21 l = 1,iorder
      do 21 ii=1,ncol
      do 21 jj=1,ncol
      do 21 kk=1,ncol
   21 sigma(ii,jj)=sigma(ii,jj)+a(ii,kk,l)*g(jj,kk,l)
      if(ldebug)then
      write(*,"('sigma')")
      call matlist(sigma,ncol,ncol)
      endif

      return
      end

      subroutine verify_var(ncol,iorder,g,a)
      real*4 g(ncol,ncol,0:iorder)
      real*4 a(ncol,ncol,iorder)
      real*4 carra(ncol,ncol,iorder,iorder)
      real*4 test1(ncol,ncol),test2(ncol,ncol),test3(ncol,ncol)
      logical*4 ldebug/.false./
      pointer (ptest1,test1),(ptest2,test2),(ptest3,test3),(pcarra,carra)

      ptest1=malloc(4*ncol**2)
      ptest2=malloc(4*ncol**2)
      ptest3=malloc(4*ncol**2)
      pcarra=malloc(4*ncol**2*iorder**2)

      do 1 i = 1,iorder
      do 2 j = i,iorder
    2 call matcop  (g(1,1,j-i),carra(1,1,i,j),ncol)
      do 3 j = 1,i-1
    3 call transpos(g(1,1,i-j),carra(1,1,i,j),ncol)
    1 continue

      if(ldebug)then
      do 21 i = 1,iorder
      do 21 j = 1,iorder
      write(*,"('carra(',i1,',',i1,')')")i,j
   21 call matlist(carra(1,1,i,j),ncol,ncol)
      endif

      do 11 i = 1,iorder
      do 12 ii = 1,ncol
      do 12 jj = 1,ncol
   12 test1(ii,jj)=0.0
      do 13 j = 1,iorder
      call transpos(a(1,1,j),test3,ncol)
      call matmul(carra(1,1,i,j),test3,test2,ncol)
      do 13 ii = 1,ncol
      do 13 jj = 1,ncol
   13 test1(ii,jj)=test1(ii,jj)+test2(ii,jj)
      call transpos(g(1,1,i),test2,ncol)
      do 14 ii = 1,ncol
      write(*,"('verify VAR')")
   14 write(*,"(16f8.4)")(test1(ii,jj),-test2(ii,jj),jj=1,ncol)
   11 continue

      call free(ptest1)
      call free(ptest2)
      call free(ptest3)
      call free(pcarra)
      return
      end

      subroutine apply_var(ncol,iorder,npts,a,format,f,z)
      character*1 format(0:npts-1)
      real*4 a(ncol,ncol,iorder),f(0:npts-1,ncol),z(0:npts-1,ncol)
      logical*4 ldebug/.false./

      if(ldebug)then
      do 13 m = 1,iorder
      write(*,"('a(',i1,')')")m
   13 call matlist(a(1,1,m),ncol,ncol)
      endif

      do 1 l = 0,npts-1
      do 2 i = 1,ncol
    2 z(l,i)=f(l,i)
      do 3 m=1,iorder
      if (l.ge.m.and.format(l-m).ne."x")then
        do i=1,ncol
        do j=1,ncol
          z(l,i)=z(l,i)+a(i,j,m)*f(l-m,j)
        enddo
        enddo
      endif
    3 continue
    1 continue

      return
      end

      subroutine format_test(format,npts)
      character*1 format(0:npts-1)
      do 1 i = 0,npts-1
    1 write(*,"(i10,4x,a1,l1)")i,format(i),format(i).eq."x"
      return
      end

      subroutine matlist(a,n,m)
      real*4 a(n,m)
      do 1 i=1,n
    1 write(*,"(16f8.4)")(a(i,j),j=1,m)
      return
      end

      subroutine compute_gc(dimx,dimy,cx,cy,cq)
      integer*4 dimx,dimy
      real*4 cx(dimx,dimx),cy(dimy,dimy),cq(dimx+dimy,dimx+dimy)
      real*4 tx(dimx,dimx),ty(dimy,dimy),tq(dimx+dimy,dimx+dimy)
      pointer (ptx,tx),(pty,ty),(ptq,tq)

      ncol=dimx+dimy
      ptx=malloc(4*dimx**2)
      pty=malloc(4*dimy**2)
      ptq=malloc(4*ncol**2)

      call matcop(cx,tx,dimx)
      call matcop(cy,ty,dimy)
      call matinv(tx,dimx,sigma1)
      call matinv(ty,dimy,tigma1)

      do 1 i=1,dimx
      do 1 j=1,dimx
    1 tx(i,j)=cq(i,j)
      do 2 i=1,dimy
      do 2 j=1,dimy
    2 ty(i,j)=cq(dimx+i,dimx+j)
      call matinv(tx,dimx,sigma2)
      call matinv(ty,dimy,tigma2)

      call matcop(cq,tq,ncol)
      call matinv(tq,ncol,detq)

c     write(*,"('sigma1,tigma1,sigma2,tigma2,detq='5f10.4)")sigma1,tigma1,sigma2,tigma2,detq
      write(*,"('#Fx,y     ',f10.6)")alog(sigma1*tigma1/detq)
      write(*,"('#Fx->y    ',f10.6)")alog(tigma1/tigma2)
      write(*,"('#Fy->x    ',f10.6)")alog(sigma1/sigma2)
      write(*,"('#Fx.y     ',f10.6)")alog(sigma2*tigma2/detq)

      call free(ptx)
      call free(pty)
      call free(ptq)
      return
      end

      subroutine compute_gc1(dimx,dimy,cx,cy,cq,tx,ty,tq,f)
      integer*4 dimx,dimy
      real*4 cx(dimx,dimx),cy(dimy,dimy),cq(dimx+dimy,dimx+dimy)
      real*4 tx(dimx,dimx),ty(dimy,dimy),tq(dimx+dimy,dimx+dimy),f(0:3)

      call matcop(cx,tx,dimx)
      call matcop(cy,ty,dimy)
      call matinv(tx,dimx,sigma1)
      call matinv(ty,dimy,tigma1)

      do 1 i=1,dimx
      do 1 j=1,dimx
    1 tx(i,j)=cq(i,j)
      do 2 i=1,dimy
      do 2 j=1,dimy
    2 ty(i,j)=cq(dimx+i,dimx+j)
      call matinv(tx,dimx,sigma2)
      call matinv(ty,dimy,tigma2)

      ncol=dimx+dimy
      call matcop(cq,tq,ncol)
      call matinv(tq,ncol,detq)

      f(0)=alog(sigma1*tigma1/detq)
      f(1)=alog(tigma1/tigma2)
      f(2)=alog(sigma1/sigma2)
      f(3)=alog(sigma2*tigma2/detq)

      return
      end

