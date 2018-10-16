c$Header: /data/petsun4/data1/src_solaris/nonlin/RCS/fvfield.f,v 1.5 2004/10/14 23:45:56 rsachs Exp $
c$Log: fvfield.f,v $
c Revision 1.5  2004/10/14  23:45:56  rsachs
c Corrected an error.
c
c Revision 1.4  2003/01/13  20:47:38  rsachs
c Corrected the error, which had voxdim(i) in the numerator rather than the denominator.
c
c Revision 1.3  2002/12/19  23:25:46  rsachs
c Modified definition of the Norm (anorm) in subroutine bodyforce, by dividin
c anorm by the total nr. of pixels (nx*ny*nz).
c
c Revision 1.2  2002/12/19  03:01:55  avi
c use voxdim in computation of bodyforce
c

      subroutine bodyforce (img1,nx,ny,nz,img2,voxdim,b)
      real*4 img1(nx,ny,nz),img2(nx,ny,nz),voxdim(3),b(nx,ny,nz,3)
      real*4 aux(3),diff

      anorm=0.0
      do 2 iz=1,nz
      do 2 iy=1,ny
      do 2 ix=1,nx
    2 anorm=anorm+img1(ix,iy,iz)**2
      anorm=anorm/(nx*ny*nz)

      do 3 i=1,3
3     aux(i)=0.5/(anorm*voxdim(i))

      do 1 iz=2,nz-1
      do 1 iy=2,ny-1
      do 1 ix=2,nx-1
         diff=img1(ix,iy,iz)-img2(ix,iy,iz)
         b(ix,iy,iz,1)=aux(1)*(img1(ix+1,iy,iz)-img1(ix-1,iy,iz))*diff
         b(ix,iy,iz,2)=aux(2)*(img1(ix,iy+1,iz)-img1(ix,iy-1,iz))*diff
         b(ix,iy,iz,3)=aux(3)*(img1(ix,iy,iz+1)-img1(ix,iy,iz-1))*diff
1     continue
      return
      end

      subroutine fainv(mu,lambda,f,ainv,det)
      real*4 mu,lambda,f(3),ainv(3,3)

      rho2 = f(1)**2+f(2)**2+f(3)**2
      if(rho2.eq.0.0)then
        do i=1,3
        do j=1,3
          ainv(i,j)=0.0
        enddo
        enddo
        det=0.0
        return
      endif
      rho4 = rho2**2
      q = mu*rho4*(2.*mu+lambda)
      det = rho2*mu*q
      c = (mu+lambda)/q
      do 1 i = 1,3
      do 1 j = 1,3
    1 ainv(i,j) = -c*f(i)*f(j)
      do 2 i = 1,3
    2 ainv(i,i) = ainv(i,i)+1./(mu*rho2)
      return
      end

      subroutine faset(mu,lambda,f,a)
      real*4 mu,lambda,f(3),a(3,3)

      do 1 i = 1,3
      do 1 j = 1,3
    1 a(i,j) = (mu+lambda)*f(i)*f(j)
      rho2 = f(1)**2+f(2)**2+f(3)**2
      do 2 i = 1,3
    2 a(i,i) = a(i,i)+mu*rho2
      return
      end

      subroutine matlst(a,n)
      real*4 a(n,n)

      do 1 i = 1,n
    1 write (*,"(8f12.4)")(a(i,j),j=1,n)
      return
      end

      subroutine mat3inv(a,ainv,det)
      real*4 a(3,3),ainv(3,3)

      det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     1     +a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))
     2     +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      ainv(1,1)= (a(2,2)*a(3,3)-a(3,2)*a(2,3))/det
      ainv(1,2)=-(a(1,2)*a(3,3)-a(3,2)*a(1,3))/det
      ainv(1,3)= (a(1,2)*a(2,3)-a(2,2)*a(1,3))/det
      ainv(2,1)=-(a(2,1)*a(3,3)-a(3,1)*a(2,3))/det
      ainv(2,2)= (a(1,1)*a(3,3)-a(3,1)*a(1,3))/det
      ainv(2,3)=-(a(1,1)*a(2,3)-a(2,1)*a(1,3))/det
      ainv(3,1)= (a(2,1)*a(3,2)-a(2,2)*a(3,1))/det
      ainv(3,2)=-(a(1,1)*a(3,2)-a(3,1)*a(1,2))/det
      ainv(3,3)= (a(1,1)*a(2,2)-a(2,1)*a(1,2))/det
      return
      end

      subroutine mat3inv_test
      real*4 a(3,3),ainv(3,3),test(3,3)

      do 2 k=1,5
      do 1 i=1,3
      do 1 j=1,3
    1 a(i,j)=rand(0)
      call mat3inv(a,ainv,det)
      call matmul(ainv,a,test,3)
      call matlst(test,3)
      write(*,"()")
      call matmul(a,ainv,test,3)
      call matlst(test,3)
    2 write(*,"()")
      call exit(0)
      end

      subroutine mbodyforce (img1,msk1,nx,ny,nz,img2,msk2,voxdim,b)
      real*4 img1(nx,ny,nz),img2(nx,ny,nz),voxdim(3),b(nx,ny,nz,3)
      real*4 aux(3),diff
      integer*2 msk1(nx,ny,nz),msk2(nx,ny,nz)

      anorm=0.0
      do 2 iz=1,nz
      do 2 iy=1,ny
      do 2 ix=1,nx
    2 anorm=anorm+img1(ix,iy,iz)**2
      anorm=anorm/(nx*ny*nz)

      do 3 i=1,3
3     aux(i)=0.5/(anorm*voxdim(i))

      do 1 iz=2,nz-1
      do 1 iy=2,ny-1
      do 1 ix=2,nx-1
      if(msk1(ix,iy,iz).eq.0.or.msk2(ix,iy,iz).eq.0) then
         b(ix,iy,iz,1) = 0.0
         b(ix,iy,iz,2) = 0.0
         b(ix,iy,iz,3) = 0.0
      else
         diff=img1(ix,iy,iz)-img2(ix,iy,iz)
         b(ix,iy,iz,1)=aux(1)*(img1(ix+1,iy,iz)-img1(ix-1,iy,iz))*diff
         b(ix,iy,iz,2)=aux(2)*(img1(ix,iy+1,iz)-img1(ix,iy-1,iz))*diff
         b(ix,iy,iz,3)=aux(3)*(img1(ix,iy,iz+1)-img1(ix,iy,iz-1))*diff
      endif
1     continue
      return
      end
