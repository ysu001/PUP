c$Header: /data/petsun4/data1/src_solaris/nonlin/RCS/inc.f,v 1.3 2004/11/01 20:42:18 rsachs Exp $
c$Log: inc.f,v $
c Revision 1.3  2004/11/01  20:42:18  rsachs
c Corrected vn(1) in sub norms1. Now it's divided by float(n).
c
c Revision 1.2  2004/10/28  21:52:08  rsachs
c Initial Revision.
c
c1    2003.01.29 - inc.f - the Incrementor module.
c1    Subs: r2s, delu, ntarg, norms1, clear, copyimg.
      subroutine r2s(i,j,k,voxdim,dx,dy,dz,r2dx,r2dy,r2dz)
      integer*4 i,j,k
      real*4 voxdim(3),dx,dy,dz,r2dx,r2dy,r2dz

      dx=voxdim(1)
      dy=voxdim(2)
      dz=voxdim(3)
      r2dx=0.5/dx
      r2dy=0.5/dy
      r2dz=0.5/dz

      if(i.eq.0) r2dx=1.0/dx    ! x=0 plane (yz-plane)
      if(j.eq.0) r2dy=1.0/dy    ! y=0   "   (xz-  "  )
      if(k.eq.0) r2dz=1.0/dz    ! z=0   "   (xy-  "  )

      return
      end                       ! r2s

      subroutine delu(v,u,du,vn,nx,ny,nz,voxdim,dt)
      real*4 v(0:nx-1,0:ny-1,0:nz-1,3),u(0:nx-1,0:ny-1,0:nz-1,3)
      real*4 du(0:nx-1,0:ny-1,0:nz-1,3),voxdim(3),vn(3)

c1    The Incrementor - computes the incremental vector displacement 
c1    change(du), given the displacement vector (u), the velocity vector
c1    (v), dx,dy,dz (specified by voxdim) and the time increment dt.
c1    Computation is done bmo the hydrodynamic derivative.
c1    Input Arguments:
c1    u,v,nx,ny,nz,voxdim,dt.
c1    Output Arguments:
c1      u,du (u is updated by adding du)

      do 10 i=0,nx-1
         ip=mod(i+1,nx)
         call mdu(i,im,dx,voxdim,1)
      do 10 j=0,ny-1
         jp=mod(j+1,ny)
         call mdu(j,jm,dy,voxdim,2)
      do 10 k=0,nz-1
         kp=mod(k+1,nz)
         call mdu(k,km,dz,voxdim,3)
         call r2s(i,j,k,voxdim,dx,dy,dz,r2dx,r2dy,r2dz) 
            VxDUxdx=v(i,j,k,1)*(u(ip,j,k,1)-u(im,j,k,1))*r2dx
            VyDUxdy=v(i,j,k,2)*(u(i,jp,k,1)-u(i,jm,k,1))*r2dy
            VzDUxdz=v(i,j,k,3)*(u(i,j,kp,1)-u(i,j,km,1))*r2dz 
            du(i,j,k,1)=(v(i,j,k,1)-(VxDUxdx+VyDUxdy+VzDUxdz))*dt
            VxDUydx=v(i,j,k,1)*(u(ip,j,k,2)-u(im,j,k,2))*r2dx
            VyDUydy=v(i,j,k,2)*(u(i,jp,k,2)-u(i,jm,k,2))*r2dy
            VzDUydz=v(i,j,k,3)*(u(i,j,kp,2)-u(i,j,km,2))*r2dz 
            du(i,j,k,2)=(v(i,j,k,2)-(VxDUydx+VyDUydy+VzDUydz))*dt
            VxDUzdx=v(i,j,k,1)*(u(ip,j,k,3)-u(im,j,k,3))*r2dx
            VyDUzdy=v(i,j,k,2)*(u(i,jp,k,3)-u(i,jm,k,3))*r2dy
            VzDUzdz=v(i,j,k,3)*(u(i,j,kp,3)-u(i,j,km,3))*r2dz
            du(i,j,k,3)=(v(i,j,k,3)-(VxDUzdx+VyDUzdy+VzDUzdz))*dt 
10    continue

      do 20 iv=1,3         ! update the displacement vector.
      do 20 i=0,nx-1
      do 20 j=0,ny-1
      do 20 k=0,nz-1
         u(i,j,k,iv)=u(i,j,k,iv)+du(i,j,k,iv)    
20    continue
      call norms1(u,vn,3*nx*ny*nz,0)
      write(*,"('delu: The time step= ', f10.4)") dt
      write(*,"('delu: Norms of the Displacement Vector ', 3(e10.3,x))") vn(1), vn(2), vn(3)

      return
      end             ! sub delu

      subroutine ntarg(u,img1,img2,voxdim,nx,ny,nz)

c1    sub to create img2 from img1 using the displacement vector u &
c1    linear interpolation.

      real*4 u(0:nx-1,0:ny-1,0:nz-1,3), voxdim(3),x(3),val
      real*4 img1(0:nx-1,0:ny-1,0:nz-1)                 ! the old image
      real*4 img2(0:nx-1,0:ny-1,0:nz-1)                 !  "  new   " 

      do 10 i=0,nx-1
      do 10 j=0,ny-1
      do 10 k=0,nz-1
         x(1)=float(i)*voxdim(1)+u(i,j,k,1)
         x(2)=float(j)*voxdim(2)+u(i,j,k,2)
         x(3)=float(k)*voxdim(3)+u(i,j,k,3)
         call imgval0(img1,nx,ny,nz,voxdim,x,val)
         img2(i,j,k)=val
10    continue

      return
      end


      subroutine norms1(v,vn,n,lprint)
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
      if(lprint.gt.0) then 
         write(*,"('The norms=', 3f10.4)") vn(1), vn(2), vn(3)
      endif
      return
      end 

      subroutine clear(u,nx,ny,nz)
      real*4 u(0:nx-1,0:ny-1,0:nz-1,3)

      do 1 iv=1,3
      do 1 i=0,nx-1
      do 1 j=0,ny-1
      do 1 k=0,nz-1
1        u(i,j,k,iv)=0.0

      return
      end 

      subroutine copyimg(img0,img1,nx,ny,nz)
      integer*4 nx,ny,nz
      real*4 img0(0:nx-1,0:ny-1,0:nz-1),img1(0:nx-1,0:ny-1,0:nz-1)

      do 1 i=0,nx-1
      do 1 j=0,ny-1
      do 1 k=0,nz-1
1     img0(i,j,k)=img1(i,j,k)

      return
      end 

      subroutine dvrgnc(v,div,voxdim,nx,ny,nz)
c1    sub to compute the divergence (div) of a given vector (v).
c1    Input Arguments:
c1       v - the vector whose divergence is to be computed.
c1       voxdim - a 3-element vector whose elements serve as dx,dy,dz.
c1       nx,ny,nz - the nr. of voxels in the x,y,z directions.
c1    Output Argument:
c1       div - the divergence - an array having same dims as v.

      real*4 v(0:nx-1,0:ny-1,0:nz-1,3),voxdim(3)
      real*4 div(0:nx-1,0:ny-1,0:nz-1),dVxdx,dVydy,dVzdz

      do 10 i=0,nx-1
      do 10 j=0,ny-1
      do 10 k=0,nz-1

      if(i.eq.0)then
         dVxdx=(v(i+1,j,k,1)-v(i,j,k,1))/voxdim(1)
      elseif(i.eq.nx-1)then
         dVxdx=(v(i,j,k,1)-v(i-1,j,k,1))/voxdim(1)
      else
         dVxdx=(v(i+1,j,k,1)-v(i-1,j,k,1))/(2*voxdim(1))
      endif 
      
      if(j.eq.0)then
         dVydy=(v(i,j+1,k,2)-v(i,j,k,2))/voxdim(2)
      elseif(j.eq.ny-1)then
         dVydy=(v(i,j,k,2)-v(i,j-1,k,2))/voxdim(2)
      else
         dVydy=(v(i,j+1,k,2)-v(i,j-1,k,2))/(2*voxdim(2))
      endif 
      
      if(k.eq.0)then
         dVzdz=(v(i,j,k+1,3)-v(i,j,k,3))/voxdim(3)
      elseif(k.eq.nz-1)then
         dVzdz=(v(i,j,k,3)-v(i,j,k-1,3))/voxdim(3)
      else
         dVzdz=(v(i,j,k+1,3)-v(i,j,k-1,3))/(2*voxdim(3))
      endif

      div(i,j,k)=dVxdx+dVydy+dVzdz
10    continue

      return
      end 
