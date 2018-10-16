      subroutine fnegdef(a)
c     force negative definite 2nd order polynomial
      real*4 a(10)
      real*4 w(3,3),wt(3,3),g(3,3),b(3,3)
      logical*4 ltest

      g(1,1)=a(10)
      g(2,2)=a(7)
      g(3,3)=a(5)
      g(1,2)=a(9)/2.
      g(1,3)=a(8)/2.
      g(2,3)=a(6)/2.
      g(2,1)=g(1,2)
      g(3,1)=g(1,3)
      g(3,2)=g(2,3)

      ltest=.false.
      call eigen(g,w,3)
               write(*,"('fnegdef: eigenvalues',3f10.6)")(g(i,i),i=1,3)
      do 1 i=1,3
      if(g(i,i).gt.0.0)then
        ltest=.true.
        g(i,i)=0.0
      endif
    1 continue
      if(ltest)write(*,"('      -> eigenvalues',3f10.6)")(g(i,i),i=1,3)

      call transpos(w,wt,3)
      call matmul(w,g,b,3)
      call matmul(b,wt,g,3)

      a(10)=g(1,1)
      a(7)= g(2,2)
      a(5)= g(3,3)
      a(9)= g(1,2)*2.
      a(8)= g(1,3)*2.
      a(6)= g(2,3)*2.

      return
      end
