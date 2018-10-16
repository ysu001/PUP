c$Header: /data/petsun4/data1/src_solaris/imglin/RCS/ft4_resolve.f,v 1.5 2009/09/09 05:51:36 avi Exp $
c$Log: ft4_resolve.f,v $
c Revision 1.5  2009/09/09  05:51:36  avi
c correct minor format statement errors
c
c Revision 1.4  1999/10/18  04:59:04  avi
c subroutine symerr()
c
c Revision 1.3  1999/10/14  03:50:10  avi
c Revision 1.2  1999/10/10  21:45:59  avi
c simplified resolve_scale()
c
c Revision 1.1  1999/10/10  07:20:39  avi
c Initial revision
c
      subroutine ft4_resolve_rcs
      write(*,"('$Id: ft4_resolve.f,v 1.5 2009/09/09 05:51:36 avi Exp $')")
      return
      end

      subroutine symerr(n,ck6)
c     compute symmetric rigid body transform error
      real*4 ck6(6,n*(n-1))
      real*4 ck6a(6),ck6b(6),ck6v(6)
      pointer (pjndex,jndex)
      integer*4 jndex(n,n)		! look up table of index corresponding to pair ij in input list

      degperck=90./atan(1.)
      pjndex=malloc(4*n**2)
      l=0
      do i=1,n
        do j=1,n
        if(j.ne.i)then
            l=l+1
            jndex(i,j)=l
        endif
        enddo
      enddo
      do 11 m=1,6						! initialize variance accumulator
   11 ck6v(m)=0.

      do 12 i=1,n
      do 12 j=1,n
      if(j.eq.i)goto 12
      call ck6mul(ck6(1,jndex(j,i)),ck6(1,jndex(i,j)),ck6a)	! Tij * Tji
      do 14 m=1,6
   14 ck6v(m)=ck6v(m)+ck6a(m)**2				! accumulate variance
   12 continue
      do 13 m=1,6
   13 ck6v(m)=ck6v(m)/float(2*n*(n-1))				! 2 for Var[] in pair
      write(*,"('rms error averaged over',i4,' pairs')")n*(n-1)
      write(*,"(6f10.4)")(degperck*sqrt(ck6v(m)),m=1,3),(sqrt(ck6v(m)),m=4,6)
      rsum=ck6v(1)+ck6v(2)+ck6v(3)
      tsum=ck6v(4)+ck6v(5)+ck6v(6)
      write(*,"('pairs total rotation error    ',f10.4,' (rms deg)')")degperck*sqrt(rsum)
      write(*,"('pairs total translation error ',f10.4,' (rms mm)')")sqrt(tsum)

      if(n.lt.3)goto 90
      do 21 m=1,6						! initialize variance accumulator
   21 ck6v(m)=0.
      do 22 i=1,n
      do 22 j=1,n
      do 22 k=1,n
      if(j.eq.i.or.k.eq.i.or.k.eq.j)goto 22
      call ck6mul(ck6(1,jndex(j,i)),ck6(1,jndex(k,j)),ck6a)	! Tij * Tjk
      call ck6mul(ck6a,ck6(1,jndex(i,k)),ck6b)			! Tij * Tjk * Tki
      do 24 m=1,6
   24 ck6v(m)=ck6v(m)+ck6b(m)**2				! accumulate variance
   22 continue
      do 23 m=1,6
   23 ck6v(m)=ck6v(m)/float(3*n*(n-1)*(n-2))			! 3 for Var[] in triple
      write(*,"('rms error averaged over',i4,' triples')")n*(n-1)*(n-2)
      write(*,"(6f10.4)")(degperck*sqrt(ck6v(m)),m=1,3),(sqrt(ck6v(m)),m=4,6)
      rsum=ck6v(1)+ck6v(2)+ck6v(3)
      tsum=ck6v(4)+ck6v(5)+ck6v(6)
      write(*,"('triples total rotation error  ',f10.4,' (rms deg)')")degperck*sqrt(rsum)
      write(*,"('triples total translation error',f9.4,' (rms mm)')")sqrt(tsum)

      if(n.lt.4)goto 90
      do 31 m=1,6						! initialize variance accumulator
   31 ck6v(m)=0.
      do 32 i=1,n
      do 32 j=1,n
      do 32 k=1,n
      do 32 l=1,n
      if(j.eq.i.or.k.eq.i.or.k.eq.j)goto 32
      if(l.eq.i.or.l.eq.j.or.l.eq.k)goto 32
      call ck6mul(ck6(1,jndex(j,i)),ck6(1,jndex(k,j)),ck6a)	! Tij * Tjk
      call ck6mul(ck6a,ck6(1,jndex(l,k)),ck6b)			! Tij * Tjk * Tkl
      call ck6mul(ck6b,ck6(1,jndex(i,l)),ck6a)			! Tij * Tjk * Tkl * Tli
      do 34 m=1,6
   34 ck6v(m)=ck6v(m)+ck6a(m)**2				! accumulate variance
   32 continue
      do 33 m=1,6
   33 ck6v(m)=ck6v(m)/float(4*n*(n-1)*(n-2)*(n-3))		! 4 for Var[] in quadruple
      write(*,"('rms error averaged over',i5,' quadruples')")n*(n-1)*(n-2)*(n-3)
      write(*,"(6f10.4)")(degperck*sqrt(ck6v(m)),m=1,3),(sqrt(ck6v(m)),m=4,6)
      rsum=ck6v(1)+ck6v(2)+ck6v(3)
      tsum=ck6v(4)+ck6v(5)+ck6v(6)
      write(*,"('quadrpls total rotation error ',f10.4,' (rms deg)')")degperck*sqrt(rsum)
      write(*,"('quadrpls total translation error',f8.4,' (rms mm)')")sqrt(tsum)

   90 call free(pjndex)
      return
      end

      subroutine resolve_rigid(radius,n,jac,err,domega,errsum)
      real*4 jac(6*n*(n-1),6*(n-1)),err(6*n*(n-1)),domega(6*(n-1))
      real*4 jtj(6*(n-1),6*(n-1)),jinv(6*(n-1),6*n*(n-1))
      real*4 w(6*(n-1),6*(n-1)),wt(6*(n-1),6*(n-1))
      real*4 ck6v(6)
      pointer (pjtj,jtj),(pjinv,jinv)
      pointer (pwt,wt),(pw,w)
      logical*4 ldebug/.false./

      m=n*(n-1)
      degperck=90./atan(1.)

      do 11 k=1,6
   11 ck6v(k)=0.
      do 12 i=0,m-1		! accumulate variance over m observations
      do 12 k=1,6
   12 ck6v(k)=ck6v(k)+err(6*i+k)**2
      do 13 k=1,6
   13 ck6v(k)=ck6v(k)/float((n-1)*(n-1))	! (n/(n-1)*(1/m)
      write(*,"('estimated observational error based on',i4,' observations')")m
      write(*,"(6f10.4)")(degperck*sqrt(ck6v(k)),k=1,3),(sqrt(ck6v(k)),k=4,6)
      rsum=ck6v(1)+ck6v(2)+ck6v(3)
      tsum=ck6v(4)+ck6v(5)+ck6v(6)
      write(*,"('estimate total rotation error   ',f8.4,' (rms deg)')")degperck*sqrt(rsum)
      write(*,"('estimate total translation error',f8.4,' (rms mm)')")sqrt(tsum)
c     write(*,"('quadrpls total translation error',f8.4,' (rms mm)')")sqrt(tsum)

      do 1 i=0,m-1
      err(6*i+1)=err(6*i+1)*2.*radius
      err(6*i+2)=err(6*i+2)*2.*radius
      err(6*i+3)=err(6*i+3)*2.*radius
      do 1 j=1,6*(n-1)
      jac(6*i+1,j)=jac(6*i+1,j)*2.*radius
      jac(6*i+2,j)=jac(6*i+2,j)*2.*radius
    1 jac(6*i+3,j)=jac(6*i+3,j)*2.*radius

      n1=min0(n-1,2)
      if(ldebug)then
        do i=1,6*m
          write(*,"(12f8.4,2x,f8.4)")(jac(i,j),j=1,6*n1),err(i)
        enddo
      endif

      errsum=0.
      do i=1,6*m
        errsum=errsum+err(i)**2
      enddo
      write(*,"('rigid body rmserr   ',f10.4)")sqrt(errsum/float(6*m))

      pjtj= malloc(4*36*(n-1)*(n-1))
      pw=   malloc(4*36*(n-1)*(n-1))
      pwt=  malloc(4*36*(n-1)*(n-1))
      pjinv=malloc(4*36*(n-1)*m)

      do 3 i=1,6*(n-1)
      do 3 j=1,6*(n-1)
      jtj(i,j)=0.
      do 3 k=1,6*m
    3 jtj(i,j)=jtj(i,j)+jac(k,i)*jac(k,j)
c     call matinv(jtj,6*(n-1),det)
c     write(*,"('JtJ determinant',e12.4)")det
      call eigen(jtj,w,6*(n-1))
      write(*,"('JtJ condition number',e10.4)")jtj(1,1)/jtj(6*(n-1),6*(n-1))

      call transpos(w,wt,6*(n-1))
      do 5 i=1,6*(n-1)
      do 5 j=1,6*(n-1)
    5 w(i,j)=w(i,j)/jtj(j,j)
      call matmul(w,wt,jtj,6*(n-1))

      do 2 i=1,6*(n-1)
      do 2 j=1,6*m
      jinv(i,j)=0.
      do 2 k=1,6*(n-1)
    2 jinv(i,j)=jinv(i,j)+jtj(i,k)*jac(j,k)
      do 4 i=1,6*(n-1)
      domega(i)=0.
      do 4 j=1,6*m
    4 domega(i)=domega(i)+jinv(i,j)*err(j)
      if(ldebug)then
        do i=1,6*n1
          write(*,"(12f8.4,2x,f8.4)")(jtj(i,j),j=1,6*n1),domega(i)
        enddo
      endif

      call free(pjtj)
      call free(pw)
      call free(pwt)
      call free(pjinv)
      return
      end

      subroutine resolve_scale(n,err,dscale)
      real*4 err(n*(n-1)),dscale(n*(n-1))
      real*4 jinv(n-1,n*(n-1))
      pointer (pjinv,jinv)
      logical*4 ldebug/.false./

      m=n*(n-1)
      t=0.
      do i=1,m
        t=t+err(i)**2
      enddo
c     write(*,"('01234567890123456789',f10.4)")sqrt(t/float(m))
      write(*,"('scale rmserr        ',f10.4)")sqrt(t/float(m))

      pjinv=malloc(4*m*(n-1))
      do 11 j=1,m
      do 11 i=1,n-1
   11 jinv(i,j)=0.
      k=0
      do 12 i=1,n
      do 12 j=1,n
      if(i.eq.j)goto 12
      k=k+1
      if(j.gt.1)jinv(j-1,k)=+0.5/float(n)
      if(i.gt.1)jinv(i-1,k)=-0.5/float(n)
      do 13 l=1,n-1
      if(i.eq.1)jinv(l,k)=jinv(l,k)+0.5/float(n)
      if(j.eq.1)jinv(l,k)=jinv(l,k)-0.5/float(n)
   13 continue
   12 continue
      if(ldebug)then
        write(*,"('jinv transposed')")
        do i=1,m
          write(*,"(8f8.4)")(jinv(j,i),j=1,n-1)
        enddo
      endif

      do 6 i=1,n-1
      dscale(i)=0.
      do 6 j=1,m
    6 dscale(i)=dscale(i)+jinv(i,j)*err(j)

      call free(pjinv)
      return
      end

      subroutine resolve_scale_test
      parameter (n=5)
      parameter (m=n*(n-1))
      real*4 jac(m,n-1)
      real*4 jtj(n-1,n-1),jinv(n-1,m)

      do 1 i=1,m
      do 1 j=1,n-1
    1 jac(i,j)=0.

      k=0
      do 2 i=1,n
      do 2 j=1,n
      if(i.eq.j)goto 2
      k=k+1
      if(j.gt.1)jac(k,j-1)=+1.
      if(i.gt.1)jac(k,i-1)=-1.
    2 continue
      if(k.ne.m)stop 'k error'
      do i=1,m
        write(*,"(8f8.4)")(jac(i,j),j=1,n-1)
      enddo

      do 3 i=1,n-1
      do 3 j=1,n-1
      jtj(i,j)=0.
      do 3 k=1,m
    3 jtj(i,j)=jtj(i,j)+jac(k,i)*jac(k,j)
      write(*,"('jtj before inversion')")
      do i=1,n-1
        write(*,"(8f8.4)")(jtj(i,j),j=1,n-1)
      enddo
      call matinv(jtj,n-1,det)
      write(*,"('Jacobian determinant',e12.4)")det
      write(*,"('jtj after inversion')")
      do i=1,n-1
        write(*,"(8f8.4)")(jtj(i,j),j=1,n-1)
      enddo
      do 5 i=1,n-1
      jtj(i,i)=1./float(n)
      do 5 j=i+1,n-1
      jtj(i,j)=0.5*jtj(1,1)
    5 jtj(j,i)=jtj(i,j)
      write(*,"('jtj after inversion')")
      do i=1,n-1
        write(*,"(8f8.4)")(jtj(i,j),j=1,n-1)
      enddo
      do 4 i=1,n-1
      do 4 j=1,m
      jinv(i,j)=0.
      do 4 k=1,n-1
    4 jinv(i,j)=jinv(i,j)+jtj(i,k)*jac(j,k)
      write(*,"('jinv transposed')")
      do i=1,m
        write(*,"(8f8.4)")(jinv(j,i),j=1,n-1)
      enddo

      do 11 j=1,m
      do 11 i=1,n-1
   11 jinv(i,j)=0.
      k=0
      do 12 i=1,n
      do 12 j=1,n
      if(i.eq.j)goto 12
      k=k+1
      if(j.gt.1)jinv(j-1,k)=+0.5/float(n)
      if(i.gt.1)jinv(i-1,k)=-0.5/float(n)
      do 13 l=1,n-1
      if(i.eq.1)jinv(l,k)=jinv(l,k)+0.5/float(n)
      if(j.eq.1)jinv(l,k)=jinv(l,k)-0.5/float(n)
   13 continue
   12 continue
      write(*,"('jinv transposed')")
      do i=1,m
        write(*,"(8f8.4)")(jinv(j,i),j=1,n-1)
      enddo

      call exit(0)
      end
