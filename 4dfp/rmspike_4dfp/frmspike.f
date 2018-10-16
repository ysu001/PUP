      subroutine frmspikev(stack,nx,ny,nz,nv,avg_vol,var_vol)
      real*4 stack(nx,ny,nz,nv),avg_vol(nx,ny,nz),var_vol(nx,ny,nz+1)

      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
        avg_vol(ix,iy,iz)=0.
        var_vol(ix,iy,iz)=0.
      enddo
      enddo
      enddo

      do iv=1,nv
        do iz=1,nz
        do iy=1,ny
        do ix=1,nx
          avg_vol(ix,iy,iz)=avg_vol(ix,iy,iz)+stack(ix,iy,iz,iv)
          var_vol(ix,iy,iz)=var_vol(ix,iy,iz)+stack(ix,iy,iz,iv)**2
        enddo
        enddo
        enddo
      enddo

      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
        avg_vol(ix,iy,iz)=avg_vol(ix,iy,iz)/float(nv)
        var_vol(ix,iy,iz)=var_vol(ix,iy,iz)/float(nv)-avg_vol(ix,iy,iz)**2
      enddo
      enddo
      enddo

      do iy=1,ny
      do ix=1,nx
        var_vol(ix,iy,nz+1)=0.
        do iz=1,nz
          var_vol(ix,iy,nz+1)=var_vol(ix,iy,nz+1)+var_vol(ix,iy,iz)
        enddo
        var_vol(ix,iy,nz+1)=var_vol(ix,iy,nz+1)/float(nz)
      enddo
      enddo

      return
      end

      subroutine frmspiker(stack,nx,ny,nz,nv,ix,iy)
      real*4 stack(nx,ny,nz,nv)

      do 1 iv=1,nv
      do 1 iz=1,nz
      t=0.25*(stack(ix-1,iy  ,iz,iv)
     &       +stack(ix  ,iy-1,iz,iv)
     &       +stack(ix+1,iy  ,iz,iv)
     &       +stack(ix  ,iy+1,iz,iv))
    1 stack(ix,iy,iz,iv)=t

      return
      end

      subroutine freorder(stack,ns,nz,nv,vol)
      real*4 stack(ns,nz,nv),vol(ns,nv)
c     reorder slices 16,8,15,7,14... => 1,2,3,4,5...

      do 1 iv=1,nv
      do 2 iz=1,16
      k=8*(1+mod(iz,2))-(iz-1)/2
      do 2 j=1,ns
    2 vol(j,iz)=stack(j,k,iv)
      do 3 iz=1,16
      do 3 j=1,ns
    3 stack(j,iz,iv)=vol(j,iz)
    1 continue

      return
      end

      
