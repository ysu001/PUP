
      subroutine mdu(m,nm,du,voxdim,ndx)
      integer*4 m,nm,ndx
      real*4 du,voxdim(3)

      if(m.ne.0) then
         nm=m-1
      else
         nm=0 
      endif
      du=voxdim(ndx)

      return
      end              ! mdu

      subroutine chkvel(v,b,err,vn,voxdim,nx,ny,nz,lambda,mu)

c1    sub to check the velocity vector (v). This is done by using the 
c1    v-vector, computing its various derivatives & plugging them in
c1    the Elliptic PDE. A perfect solution should yield 0 on every line.
c1    Input  Arguments:
c1               v,b - the velocity & body force vectors respectively.
c1    Output Arguments:
c1               the 3 norms of the resulting vector, in the vector vn.

      real*4 v(0:nx-1,0:ny-1,0:nz-1,3),b(0:nx-1,0:ny-1,0:nz-1,3)
      real*4 err(0:nx-1,0:ny-1,0:nz-1,3),voxdim(3),vn(3),lambda,mu,lpm
      real*4 dx,dy,dz,r4dxdy,r4dxdz,r4dydz,xx(3)

      lpm=lambda+mu

      do 1 iv=1,3
1     vn(iv)=0.0

      do 10 i=0,nx-1
         ip=mod(i+1,nx)
         call mdu(i,im,dx,voxdim,1)
         rdx2=1.0/(dx*dx)
      do 10 j=0,ny-1
         jp=mod(j+1,ny)
         call mdu(j,jm,dy,voxdim,2) 
         rdy2=1.0/(dy*dy)
      do 10 k=0,nz-1
         kp=mod(k+1,nz)
         call mdu(k,km,dz,voxdim,3) 
         rdz2=1.0/(dz*dz)
         call rs(i,j,k,voxdim,dx,dy,dz,r4dxdy,r4dxdz,r4dydz)
         if((i*j*k).gt.0) then
           D2Vxdx2 =(v(ip,j,k,1)-2.*v(i,j,k,1)+v(im,j,k,1))*rdx2
           D2Vxdy2 =(v(i,jp,k,1)-2.*v(i,j,k,1)+v(i,jm,k,1))*rdy2
           D2Vxdz2 =(v(i,j,kp,1)-2.*v(i,j,k,1)+v(i,j,km,1))*rdz2
           D2Vydx2 =(v(ip,j,k,2)-2.*v(i,j,k,2)+v(im,j,k,2))*rdx2
           D2Vydy2 =(v(i,jp,k,2)-2.*v(i,j,k,2)+v(i,jm,k,2))*rdy2
           D2Vydz2 =(v(i,j,kp,2)-2.*v(i,j,k,2)+v(i,j,km,2))*rdz2
           D2Vzdx2 =(v(ip,j,k,3)-2.*v(i,j,k,3)+v(im,j,k,3))*rdx2
           D2Vzdy2 =(v(i,jp,k,3)-2.*v(i,j,k,3)+v(i,jm,k,3))*rdy2
           D2Vzdz2 =(v(i,j,kp,3)-2.*v(i,j,k,3)+v(i,j,km,3))*rdz2
         else
           D2Vxdx2 =0.0
           D2Vxdy2 =0.0
           D2Vxdz2 =0.0
           D2Vydx2 =0.0
           D2Vydy2 =0.0
           D2Vydz2 =0.0
           D2Vzdx2 =0.0
           D2Vzdy2 =0.0
           D2Vzdz2 =0.0 
         endif
           D2Vydxdy=(v(ip,jp,k,2)-v(ip,jm,k,2)-v(im,jp,k,2)+
     1               v(im,jm,k,2))*r4dxdy     
           D2Vzdxdz=(v(ip,j,kp,3)-v(ip,j,km,3)-v(im,j,kp,3)+
     1               v(im,j,km,3))*r4dxdz
           xx(1)=mu*(D2Vxdx2+D2Vxdy2+D2Vxdz2)+lpm*(D2Vxdx2
     1                 +D2Vydxdy+D2Vzdxdz)+b(i,j,k,1)
           D2Vxdxdy=(v(ip,jp,k,1)-v(ip,jm,k,1)-v(im,jp,k,1)+
     1               v(im,jm,k,1))*r4dxdy   
           D2Vzdydz=(v(i,jp,kp,3)-v(i,jp,km,3)-v(i,jm,kp,3)+
     1               v(i,jm,km,3))*r4dydz   
           xx(2)=mu*(D2Vydx2+D2Vydy2+D2Vydz2)+lpm*(D2Vxdxdy
     1                  +D2Vydy2+D2Vzdydz)+b(i,j,k,2)
           D2Vxdxdz=(v(ip,j,kp,1)-v(ip,j,km,1)-v(im,j,kp,1)+
     1               v(im,j,km,1))*r4dxdz    
           D2Vydydz=(v(i,jp,kp,2)-v(i,jp,km,2)-v(i,jm,kp,2)+
     1               v(i,jm,km,2))*r4dydz       
           xx(3)=mu*(D2Vzdx2+D2Vzdy2+D2Vzdz2)+lpm*(D2Vxdxdz
     1               +D2Vydydz+D2Vzdz2)+b(i,j,k,3)
         do iv=1,3
            err(i,j,k,iv)=xx(iv)
            xx(iv)=abs(xx(iv))
            vn(1)=vn(1)+xx(iv)
            vn(2)=vn(2)+xx(iv)*xx(iv)
            if(xx(iv).gt.vn(3)) vn(3)=xx(iv)
         enddo
10    continue
      m=nx*ny*nz
      vn(1)=vn(1)/float(m)
      vn(2)=sqrt(vn(2))/float(m)
      write(*,"('chkvel: Norms of the Solution  ',3(e10.4,3x))"),vn(1),vn(2),vn(3)

      return
      end              ! sub chkvel 

C---------------------------------------

      subroutine norms(v,vn,n)
      real*4 v(n),vn(3)
      integer*4 n

      do 1 i=1,3
1     vn(i)=0.0

      do i=1,n
         x=abs(v(i))
         vn(1)=vn(1)+x
         vn(2)=vn(2)+x*x
         if(x.gt.vn(3)) vn(3)=x
      enddo
      vn(1)=vn(1)/float(n)
      vn(2)=sqrt(vn(2))/float(n)
c1      write(*,"('Nr of variables= ', i8)") n
      write(*,"('The norms=', 3e13.4)") vn(1), vn(2), vn(3)
      return
      end 
C---------------------------------------------

      subroutine rs1(i,j,k,voxdim,dx,dy,dz,r4dxdy,r4dxdz,r4dydz)
      integer*4 i,j,k
      real*4 voxdim(3),dx,dy,dz,r4dxdy,r4dxdz,r4dydz

      dx=voxdim(1)
      dy=voxdim(2)
      dz=voxdim(3)
      r4dxdy=1.0/(4.0*dx*dy) 
      r4dxdz=1.0/(4.0*dx*dz)
      r4dydz=1.0/(4.0*dy*dz)

      if(i.eq.0.and.j.eq.0) then           ! z-axis
         r4dxdz=1.0/(dx*dz)
         r4dydz=1.0/(dy*dz)
      endif
      if(i.eq.0.and.k.eq.0) then           ! y-axis
         r4dxdy=1.0/(dx*dy)
         r4dydz=1.0/(dy*dz)
      endif
      if(j.eq.0.and.k.eq.0) then           ! x-axis
         r4dxdy=1.0/(dx*dy)
         r4dxdz=1.0/(dx*dz)
      endif
  
      return
      end            ! rs1

      subroutine rs(i,j,k,voxdim,dx,dy,dz,r4dxdy,r4dxdz,r4dydz)
      integer*4 i,j,k
      real*4 voxdim(3),dx,dy,dz,r4dxdy,r4dxdz,r4dydz

      dx=voxdim(1)
      dy=voxdim(2)
      dz=voxdim(3)
      r4dxdy=1.0/(4.0*dx*dy) 
      r4dxdz=1.0/(4.0*dx*dz)
      r4dydz=1.0/(4.0*dy*dz)

      if(i.eq.0) then           ! the x=0 plane (yz plane)
         r4dxdy=1.0/(dx*dy)
         r4dxdz=1.0/(dx*dz)
      endif
      if(j.eq.0) then           ! the y=0 plane (xz plane)
         r4dxdy=1.0/(dx*dy)
         r4dydz=1.0/(dy*dz)
      endif
      if(k.eq.0) then           ! the z=0 plane (xy plane)
         r4dxdz=1.0/(dx*dz)
         r4dydz=1.0/(dy*dz)
      endif
  
      return
      end            ! rs

