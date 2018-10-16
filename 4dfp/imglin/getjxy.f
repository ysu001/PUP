      subroutine getjxy_testd
      real*4 txa(6),tya(6),txy(6),txy1(6),terr(6),terr1(6),that(6)
      real*4 j_xy_xa(6,6),j_xy_ya(6,6),dummy(36)
      real*4 del/0.001/

      do 9 l=1,3
      do k=1,3
        txa(k)=rand(0)-0.5
        tya(k)=rand(0)-0.5
        txa(3+k)=10.*(rand(0)-0.5)
        tya(3+k)=10.*(rand(0)-0.5)
      enddo
      write(*,"('Tax       ',3f10.6,3f10.4)")txa
      write(*,"('Tay       ',3f10.6,3f10.4)")tya
      call gettxy(txa,tya,that,dummy,dummy)
      do k=1,3
        that(k)=  that(k)  +0.5*(rand(0)-0.5)
        that(3+k)=that(3+k)+9.0*(rand(0)-0.5)
      enddo
      write(*,"('That      ',3f10.6,3f10.4)")that
      call getjxy(txa,tya,txy,that,terr,j_xy_xa,j_xy_ya)
      write(*,"('Txy       ',3f10.6,3f10.4)")txy
      write(*,"('Terr      ',3f10.6,3f10.4)")terr
      do 8 j=1,6
      write(*,"('J_xy_xa',i3,3f10.6,3f10.4)")j,(j_xy_xa(k,j),k=1,6)
      txa(j)=txa(j)+del
      call getjxy(txa,tya,txy1,that,terr1,dummy,dummy)
      txa(j)=txa(j)-del
    8 write(*,"(10x,         3f10.6,3f10.4)")((terr1(k)-terr(k))/del,k=1,6)
      do 7 j=1,6
      write(*,"('J_xy_ya',i3,3f10.6,3f10.4)")j,(j_xy_ya(k,j),k=1,6)
      tya(j)=tya(j)+del
      call getjxy(txa,tya,txy1,that,terr1,dummy,dummy)
      tya(j)=tya(j)-del
    7 write(*,"(10x,         3f10.6,3f10.4)")((terr1(k)-terr(k))/del,k=1,6)
    9 write(*,"()")

      call exit(0)
      end

      subroutine getjxy(txa,tya,txy,that,terr,j_xy_xa,j_xy_ya)
      real*4 txa(6),tya(6),txy(6),that(6),terr(6)
      real*4 j_xy_xa(6,6),j_xy_ya(6,6)
      real*4 exa(0:3),eay(0:3),exy(0:3),dexydexa(0:3,3),dexydeay(0:3,3)
      real*4 e4hat(0:3,0:3),ehat(0:3)
      real*4 rxa(3,3),drxa(3,3,3)
      real*4 ray(3,3),dray(3,3,3)
      real*4 rxy(3,3),drxydexa(3,3,3),drxydeay(3,3,3)
      logical*4 ldebug/.false./

      if(ldebug)write(*,"(3f10.6,3f10.4,'  txa')")txa
      if(ldebug)write(*,"(3f10.6,3f10.4,'  tya')")tya

      do 2 k=1,3
      exa(k)= txa(k)
      eay(k)=-tya(k)
    2 ehat(k)=that(k)
      ehat(0)=sqrt(1.0-ehat(1)**2-ehat(2)**2-ehat(3)**2)
      do 3 l=0,3
      e4hat(0,l)= ehat(l)
    3 e4hat(l,l)= ehat(0)
      e4hat(1,2)=-ehat(3)
      e4hat(1,3)= ehat(2)
      e4hat(2,3)=-ehat(1)
      do 4 l=0,2
      do 4 m=l+1,3
    4 e4hat(m,l)=-e4hat(l,m)
      if(.false.)then
        do l=0,3
          write(*,"(4f10.6)")(e4hat(l,m),m=0,3)
        enddo
        call exit(1)
      endif

      call ckmuld(exa,eay,exy,dexydexa,dexydeay)
      do 6 k=1,3
      terr(k)=0.
      do 6 m=0,3
    6 terr(k)=terr(k)+e4hat(k,m)*exy(m)

      do 1 k=1,6
      do 1 j=1,6
      j_xy_xa(k,j)=0.
    1 j_xy_ya(k,j)=0.

      call e2rd(exa,rxa,drxa)
      call e2rd(eay,ray,dray)
      call matmul(rxa,ray,rxy,3)
      do j=1,3
        call matmul(drxa(1,1,j),ray,drxydexa(1,1,j),3)
        call matmul(rxa,dray(1,1,j),drxydeay(1,1,j),3)
      enddo
      do k=1,3
        txy(k)=exy(k)
        txy (3+k)=txa(3+k)-rxy(k,1)*tya(4)-rxy(k,2)*tya(5)-rxy(k,3)*tya(6)
        terr(3+k)=txy(3+k)-that(3+k)
        do j=1,3
          do m=0,3
            j_xy_xa(k,j)=j_xy_xa(k,j)+e4hat(k,m)*dexydexa(m,j)
            j_xy_ya(k,j)=j_xy_ya(k,j)-e4hat(k,m)*dexydeay(m,j)
          enddo
        enddo
        do j=1,3
          j_xy_xa(3+k,j)=-drxydexa(k,1,j)*tya(4)-drxydexa(k,2,j)*tya(5)-drxydexa(k,3,j)*tya(6)
          j_xy_ya(3+k,j)= drxydeay(k,1,j)*tya(4)+drxydeay(k,2,j)*tya(5)+drxydeay(k,3,j)*tya(6)
          j_xy_ya(3+k,3+j)=-rxy(k,j)
        enddo
        j_xy_xa(3+k,3+k)=1.
      enddo
      if(ldebug)write(*,"(3f10.6,3f10.4,'  txy')")txy

      return
      end
