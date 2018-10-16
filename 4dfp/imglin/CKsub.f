c$Header: /data/petsun4/data1/src_solaris/imglin/RCS/CKsub.f,v 1.5 1999/10/18 02:20:51 avi Exp $
c$Log: CKsub.f,v $
c Revision 1.5  1999/10/18  02:20:51  avi
c additional subroutines
c
c Revision 1.4  1999/10/14  03:48:39  avi
c subroutine ck6_to_t4()
c
c Revision 1.3  1999/10/08  06:58:56  avi
c subroutines gettxy and gettxy1
c
c Revision 1.2  1999/10/08  01:10:11  avi
      subroutine cksub_rcs
      write(*,"('$Id: CKsub.f,v 1.5 1999/10/18 02:20:51 avi Exp $')")
      return
      end

      subroutine gettxy_testd
      real*4 txa(6),tya(6),txy(6),txy1(6)
      real*4 dtxydtxa(6,6),dtxydtya(6,6),dummy(36)
      real*4 del/0.001/

      do 9 l=1,3
      do k=1,3
        txa(k)=rand(0)-0.5
        tya(k)=rand(0)-0.5
        txa(3+k)=10.*(rand(0)-0.5)
        tya(3+k)=10.*(rand(0)-0.5)
      enddo
      call gettxy(txa,tya,txy,dtxydtxa,dtxydtya)
      write(*,"('Tax       ',3f10.6,3f10.4)")txa
      write(*,"('Tay       ',3f10.6,3f10.4)")tya
      write(*,"('Txy       ',3f10.6,3f10.4)")txy
      do 8 j=1,6
      write(*,"('dTxy/dtxa',i1,3f10.6,3f10.4)")j,(dtxydtxa(k,j),k=1,6)
      txa(j)=txa(j)+del
      call gettxy(txa,tya,txy1,dummy,dummy)
      txa(j)=txa(j)-del
    8 write(*,"(10x,         3f10.6,3f10.4)")((txy1(k)-txy(k))/del,k=1,6)
      do 7 j=1,6
      write(*,"('dTxy/dtya',i1,3f10.6,3f10.4)")j,(dtxydtya(k,j),k=1,6)
      tya(j)=tya(j)+del
      call gettxy(txa,tya,txy1,dummy,dummy)
      tya(j)=tya(j)-del
    7 write(*,"(10x,         3f10.6,3f10.4)")((txy1(k)-txy(k))/del,k=1,6)
    9 write(*,"()")

      call exit(0)
      end

      subroutine gettxy1_testd
      real*4 tax(6),tay(6),txy(6),txy1(6)
      real*4 dtxydtax(6,6),dtxydtay(6,6),dummy(36)
      real*4 del/0.001/

      do 9 l=1,3
      do k=1,3
        tax(k)=rand(0)-0.5
        tay(k)=rand(0)-0.5
        tax(3+k)=10.*(rand(0)-0.5)
        tay(3+k)=10.*(rand(0)-0.5)
      enddo
      call gettxy1(tax,tay,txy,dtxydtax,dtxydtay)
      write(*,"('Tax       ',3f10.6,3f10.4)")tax
      write(*,"('Tay       ',3f10.6,3f10.4)")tay
      write(*,"('Txy       ',3f10.6,3f10.4)")txy
      do 8 j=1,6
      write(*,"('dTxy/dtax',i1,3f10.6,3f10.4)")j,(dtxydtax(k,j),k=1,6)
      tax(j)=tax(j)+del
      call gettxy1(tax,tay,txy1,dummy,dummy)
      tax(j)=tax(j)-del
    8 write(*,"(10x,         3f10.6,3f10.4)")((txy1(k)-txy(k))/del,k=1,6)
      do 7 j=1,6
      write(*,"('dTxy/dtay',i1,3f10.6,3f10.4)")j,(dtxydtay(k,j),k=1,6)
      tay(j)=tay(j)+del
      call gettxy1(tax,tay,txy1,dummy,dummy)
      tay(j)=tay(j)-del
    7 write(*,"(10x,         3f10.6,3f10.4)")((txy1(k)-txy(k))/del,k=1,6)
    9 write(*,"()")

      call exit(0)
      end

      subroutine ckmul_testd
      real*4 r(3,3),e(0:3),dr(3,3,3)
      real*4 a(3,3),b(3,3),c(3,3,3)
      data del/0.001/

      do 9 l=1,10
      do k=1,3
        e(k)=rand(0)-0.5
      enddo
      call e2rd(e,r,dr)
      do j=1,3
        e(j)=e(j)+0.5*del
        call e2r(e,b)
        e(j)=e(j)-1.0*del
        call e2r(e,a)
        e(j)=e(j)+0.5*del
        call rotdif(del,a,b,c(1,1,j))
      enddo
      do 8 k=1,3
      do 8 i=1,3
    8 write(*,"(9f10.6)")(dr(i,j,k),c(i,j,k),dr(i,j,k)-c(i,j,k),j=1,3)
    9 write(*,"()")

      call exit(0)
      end

      subroutine rotdif(del,a,b,c)
      real*4 a(3,3),b(3,3),c(3,3)
      do 1 i=1,3
      do 1 j=1,3
    1 c(i,j)=(b(i,j)-a(i,j))/del
      return
      end

      subroutine ckmul_test1
      real*4 e(0:3),f(0:3),g(0:3),gt(0:3)
      real*4 dgde(0:3,3),dgdf(0:3,3),a(0:3,3),b(0:3,3)
      real*4 del/0.001/

      do 9 l=1,10
      do k=1,3
        e(k)=rand(0)-0.5
        f(k)=rand(0)-0.5
      enddo
      call ckmuld(e,f,g,dgde,dgdf)
      do j=1,3
        e(j)=e(j)+0.5*del
        call ckmul(e,f,gt)
        do i=0,3
          a(i,j)=gt(i)
        enddo
        e(j)=e(j)-1.0*del
        call ckmul(e,f,gt)
        do i=0,3
          a(i,j)=(a(i,j)-gt(i))/del
        enddo
        e(j)=e(j)+0.5*del
      enddo
      do j=1,3
        f(j)=f(j)+0.5*del
        call ckmul(e,f,gt)
        do i=0,3
          b(i,j)=gt(i)
        enddo
        f(j)=f(j)-1.0*del
        call ckmul(e,f,gt)
        do i=0,3
          b(i,j)=(b(i,j)-gt(i))/del
        enddo
        f(j)=f(j)+0.5*del
      enddo
      do 8 i=0,3
    8 write(*,"(9f10.6)")(dgde(i,j),a(i,j),dgde(i,j)-a(i,j),j=1,3)
      do 7 i=0,3
    7 write(*,"(9f10.6)")(dgdf(i,j),b(i,j),dgdf(i,j)-b(i,j),j=1,3)
    9 write(*,"()")

      call exit(0)
      end

      subroutine ckmul_test
      real*4 e(0:3),f(0:3),g(0:3)
      real*4 r(3,3),q(3,3),p(3,3)

      do 9 l=1,10
      do k=1,3
        e(k)=rand(0)-0.5
        f(k)=rand(0)-0.5
      enddo
      call ckmul(e,f,g)
      call e2r(e,p)
      call e2r(f,q)
      call matmul(p,q,r,3)
      call rotlst(r)
      call e2r(g,p)
      call rotlst(p)
    9 write(*,"()")

      call exit(0)
      end

      subroutine gettxy(txa,tya,txy,dtxydtxa,dtxydtya)
      real*4 txa(6),tya(6),txy(6)
      real*4 dtxydtxa(6,6),dtxydtya(6,6)
      real*4 exa(0:3),eay(0:3),exy(0:3),dexydexa(0:3,3),dexydeay(0:3,3)
      real*4 rxa(3,3),drxa(3,3,3)
      real*4 ray(3,3),dray(3,3,3)
      real*4 rxy(3,3),drxydexa(3,3,3),drxydeay(3,3,3)
      logical*4 ldebug/.false./

      if(ldebug)write(*,"(3f10.6,3f10.4,'  txa')")txa
      if(ldebug)write(*,"(3f10.6,3f10.4,'  tya')")tya
      do 1 k=1,6
      do 1 j=1,6
      dtxydtxa(k,j)=0.
    1 dtxydtya(k,j)=0.

      do k=1,3
        exa(k)= txa(k)
        eay(k)=-tya(k)
      enddo
      call ckmuld(exa,eay,exy,dexydexa,dexydeay)
      call e2rd(exa,rxa,drxa)
      call e2rd(eay,ray,dray)
      call matmul(rxa,ray,rxy,3)
      do j =1,3
        call matmul(drxa(1,1,j),ray,drxydexa(1,1,j),3)
        call matmul(rxa,dray(1,1,j),drxydeay(1,1,j),3)
      enddo
      do k=1,3
        txy(k)=exy(k)
        txy(3+k)=txa(3+k)-rxy(k,1)*tya(4)-rxy(k,2)*tya(5)-rxy(k,3)*tya(6)
        do j=1,3
          dtxydtxa(k,j)= dexydexa(k,j)
          dtxydtya(k,j)=-dexydeay(k,j)
          dtxydtxa(3+k,j)=-drxydexa(k,1,j)*tya(4)-drxydexa(k,2,j)*tya(5)-drxydexa(k,3,j)*tya(6)
          dtxydtya(3+k,j)= drxydeay(k,1,j)*tya(4)+drxydeay(k,2,j)*tya(5)+drxydeay(k,3,j)*tya(6)
          dtxydtya(3+k,3+j)=-rxy(k,j)
        enddo
        dtxydtxa(3+k,3+k)=1.
      enddo
      if(ldebug)write(*,"(3f10.6,3f10.4,'  txy')")txy

      return
      end

      subroutine gettxy1(tax,tay,txy,dtxydtax,dtxydtay)
      real*4 tax(6),tay(6),txy(6)
      real*4 dtxydtax(6,6),dtxydtay(6,6)
      real*4 exa(0:3),eay(0:3),exy(0:3),dexydexa(0:3,3),dexydeay(0:3,3)
      real*4 rxa(3,3),drxa(3,3,3)

      do 1 k=1,6
      do 1 j=1,6
      dtxydtax(k,j)=0.
    1 dtxydtay(k,j)=0.

      do k=1,3
        exa(k)=-tax(k)
        eay(k)= tay(k)
      enddo
      call ckmuld(exa,eay,exy,dexydexa,dexydeay)
      call e2rd(exa,rxa,drxa)
      do k=1,3
        txy(k)=exy(k)
        txy(3+k)=rxa(k,1)*(-tax(4)+tay(4))+rxa(k,2)*(-tax(5)+tay(5))+rxa(k,3)*(-tax(6)+tay(6))
        do j=1,3
          dtxydtax(k,j)=-dexydexa(k,j)
          dtxydtay(k,j)= dexydeay(k,j)
          dtxydtax(3+k,j)=-drxa(k,1,j)*(-tax(4)+tay(4))
     &                    -drxa(k,2,j)*(-tax(5)+tay(5))
     &                    -drxa(k,3,j)*(-tax(6)+tay(6))
          dtxydtax(3+k,3+j)=-rxa(k,j)
          dtxydtay(3+k,3+j)= rxa(k,j)
        enddo
      enddo

      return
      end

      subroutine ckmuld(e,f,g,dgde,dgdf)
      real*4 e(0:3),f(0:3),g(0:3)
      real*4 dgde(0:3,3),dgdf(0:3,3)
      integer*4 klc(3,3)/0,  3,  2,  3, 0,  1,  2,  1, 0/
      real*4 slc(3,3)/   0.,+1.,-1.,-1.,0.,+1.,+1.,-1.,0./

      call ckmul(e,f,g)
      do 1 j=1,3
      dgde(0,j)=-f(j)           -f(0)*e(j)/e(0)
      dgdf(0,j)=-e(j)           -e(0)*f(j)/f(0)
      do 1 i=1,3
      if(i.eq.j)then
        dgde(i,j)=f(0)          -f(i)*e(j)/e(0)
        dgdf(i,j)=e(0)          -e(i)*f(j)/f(0)
      else
        l=klc(i,j)
        dgde(i,j)= slc(i,j)*f(l)-f(i)*e(j)/e(0)
        dgdf(i,j)=-slc(i,j)*e(l)-e(i)*f(j)/f(0)
      endif
    1 continue
    
      return
      end

      subroutine ck6_to_param6(ed,param6)
      real*4 ed(6),param6(6)
      real*4 e(0:3),r(3,3)

      e(1)= ed(1)
      e(2)=-ed(2)
      e(3)=-ed(3)
      call e2r(e,r)		! in CKsub.o
      call rot2ang(r,param6(4))	! in librms.a
      param6(1)= ed(4)
      param6(2)=-ed(5)
      param6(3)=-ed(6)

      return
      end

      subroutine ck6_to_t4(ed,t4)
      reaL*4 t4(4,4),ed(6)
      real*4 r(3,3),e(0:3)

      call t4_init(t4)
      do 3 k=1,3
      e(k)=ed(k)
    3 t4(k,4)=ed(3+k)
      call e2r(e,r)
      do 1 i=1,3
      do 1 j=1,3
    1 t4(i,j)=r(i,j)

      return
      end

      subroutine t4_to_ck6(t4,ed)
      reaL*4 t4(4,4),ed(6)
      real*4 r(3,3),e(0:3)

      do 1 i=1,3
      do 1 j=1,3
    1 r(i,j)=t4(i,j)
      call r2e(r,e)
      do 2 k=1,3
      ed(k)=e(k)
    2 ed(3+k)=t4(k,4)

      return
      end

      subroutine ck6mul(a,b,c)
      real*4 a(6),b(6),c(6)
      real*4 ta(4,4),tb(4,4),tc(4,4)

      call ck6_to_t4(a,ta)
      call ck6_to_t4(b,tb)
      call matmul(ta,tb,tc,4)
      call t4_to_ck6(tc,c)

      return
      end

      subroutine ckmul(e,f,g)
      real*4 e(0:3),f(0:3),g(0:3)

      e(0)=sqrt(1.-e(1)**2-e(2)**2-e(3)**2)
      f(0)=sqrt(1.-f(1)**2-f(2)**2-f(3)**2)
      g(0)=f(0)*e(0)-f(1)*e(1)-f(2)*e(2)-f(3)*e(3)
      g(1)=f(0)*e(1)+f(1)*e(0)+f(2)*e(3)-f(3)*e(2)
      g(2)=f(0)*e(2)-f(1)*e(3)+f(2)*e(0)+f(3)*e(1)
      g(3)=f(0)*e(3)+f(1)*e(2)-f(2)*e(1)+f(3)*e(0)
      return
      end

      subroutine r2e(r,e)
      real*4 r(3,3),e(0:3)

      e(0)=sqrt(r(1,1)+r(2,2)+r(3,3)+1.)/2.
      e(1)=0.25*(r(2,3)-r(3,2))/e(0)
      e(2)=0.25*(r(3,1)-r(1,3))/e(0)
      e(3)=0.25*(r(1,2)-r(2,1))/e(0)
      return
      end

      subroutine e2r(e,r)
      real*4 r(3,3),e(0:3)

      e(0)=sqrt(1.-e(1)**2-e(2)**2-e(3)**2)
      r(1,1)=e(0)**2+e(1)**2-e(2)**2-e(3)**2
      r(2,2)=e(0)**2-e(1)**2+e(2)**2-e(3)**2
      r(3,3)=e(0)**2-e(1)**2-e(2)**2+e(3)**2
      r(1,2)=2.*(e(1)*e(2)+e(0)*e(3))
      r(2,1)=2.*(e(1)*e(2)-e(0)*e(3))
      r(1,3)=2.*(e(1)*e(3)-e(0)*e(2))
      r(3,1)=2.*(e(1)*e(3)+e(0)*e(2))
      r(2,3)=2.*(e(2)*e(3)+e(0)*e(1))
      r(3,2)=2.*(e(2)*e(3)-e(0)*e(1))
      return
      end

      subroutine e2rd(e,r,dr)
      real*4 r(3,3),e(0:3),dr(3,3,3)

      call e2r(e,r)
      dr(1,1,1)=0.
      dr(1,2,1)=-r(1,3)/e(0)
      dr(1,3,1)= r(1,2)/e(0)
      dr(2,1,1)= r(3,1)/e(0)
      dr(2,2,1)=-4.*e(1)
      dr(2,3,1)= (r(2,2)+r(3,3))/e(0)
      dr(3,1,1)=-r(2,1)/e(0)
      dr(3,2,1)=-dr(2,3,1)
      dr(3,3,1)= dr(2,2,1)

      dr(1,1,2)=-4.*e(2)
      dr(1,2,2)=-r(3,2)/e(0)
      dr(1,3,2)=-(r(1,1)+r(3,3))/e(0)
      dr(2,1,2)= r(2,3)/e(0)
      dr(2,2,2)=0.
      dr(2,3,2)= dr(3,1,1)
      dr(3,1,2)=-dr(1,3,2)
      dr(3,2,2)= dr(1,3,1)
      dr(3,3,2)= dr(1,1,2)

      dr(1,1,3)=-4.*e(3)
      dr(1,2,3)= (r(1,1)+r(2,2))/e(0)
      dr(1,3,3)= dr(2,1,2)
      dr(2,1,3)=-dr(1,2,3)
      dr(2,2,3)= dr(1,1,3)
      dr(2,3,3)= dr(1,2,1)
      dr(3,1,3)= dr(1,2,2)
      dr(3,2,3)= dr(2,1,1)
      dr(3,3,3)=0.

      return
      end

      subroutine rotlst(r)
      real*4 r(3,3)
      do 1 i=1,3
    1 write(*,"(3f10.6)")(r(i,j),j=1,3)
      return
      end
