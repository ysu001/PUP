      subroutine fsquash(stack,nvox,nvol,imag)
      real*4 stack(nvox,nvol),imag(nvox)

      do 2 i=1,nvox
      imag(i)=0.
      do 1 j=1,nvol
    1 imag(i)=imag(i)+stack(i,j)
    2 imag(i)=imag(i)/float(nvol)
      return
      end

      subroutine fts_var(stack,np,nv,mask,vol,var)
      real*4 stack(np,nv),vol(np),var(nv)
      integer*2 mask(np)

      do iv=1,nv
        var(iv)=0.
        n=0
        do ip=1,np
          if(mask(ip).gt.0)then
            var(iv)=var(iv)+(stack(ip,iv)-vol(ip))**2
            n=n+1
          endif
        enddo
        var(iv)=var(iv)/float(n)
      enddo

      return
      end

      subroutine frmo(stack,np,nv)
      real*4 stack(np,nv)

      do ip=1,np
        u0=0.
        do iv=1,nv
          u0=u0+stack(ip,iv)
        enddo
        u0=u0/float(nv)
        do iv=1,nv
          stack(ip,iv)=stack(ip,iv)-u0
        enddo
      enddo

      return
      end

      subroutine fslice_means(imag,np,nz,mask,n_z,v_z)
      real*4 imag(np,nz),v_z(nz)
      integer*2 mask(np,nz)
      integer*4 n_z(nz)

      do iz=1,nz
        v_z(iz)=0.
        n_z(iz)=0
        do ip=1,np
          if(mask(ip,iz).gt.0)then
            v=imag(ip,iz)
            v_z(iz)=v_z(iz)+v
            n_z(iz)=n_z(iz)+1
          endif
        enddo
        if(n_z(iz).gt.0)v_z(iz)=v_z(iz)/float(n_z(iz))
c       write(*,"(i10,f10.4,i10)")iz,v_z(iz),n_z(iz)
      enddo

      return
      end

      subroutine fnorm_volume(stack,np,nv,mask,vol,amean)
      real*4 stack(np,nv)	! np = pixels/volume
      real*4 vol(np)		! functional frame squashed volume (brain masked)
      integer*2 mask(np)

      amean=0.
      n=0
      do 1 ip=1,np
      if(mask(ip).gt.0)then
        amean=amean+vol(ip)
        n=n+1
      endif
    1 continue
      amean=amean/float(n)
      write(*,"('Normalizing functional frames to brain mean=',f10.4)")amean

      do 2 iv=1,nv
      sum=0.
      do ip=1,np
        if(mask(ip).gt.0)then
          sum=sum+stack(ip,iv)
        endif
      enddo
      f=amean*float(n)/sum
      do ip=1,np
        stack(ip,iv)=f*stack(ip,iv)
      enddo
    2 continue

      return
      end

      subroutine fnorm_slice(stack,np,nz,nv,mask,v_z)
      real*4 stack(np,nz,nv),v_z(nz)
      integer*2 mask(np,nz)

      do 1 iz=1,nz
      write(*,"('Normalizing slice',i3,' to mean=',f10.4)")iz,v_z(iz)
      do 2 iv=1,nv
      sum=0.
      n=0
      do ip=1,np
        if(mask(ip,iz).gt.0)then
          sum=sum+stack(ip,iz,iv)
          n=n+1
        endif
      enddo
      if(n.gt.0)then
        f=v_z(iz)*float(n)/sum
        do ip=1,np
          stack(ip,iz,iv)=f*stack(ip,iz,iv)
        enddo
      endif
    2 continue
    1 continue

      return
      end

      subroutine ffind_mode(hist,nbin,binval,amode)
      integer*4 hist(nbin)
      real*4 binval(nbin)
      parameter (nbinmax=1024)
      real*4 thist(nbinmax),shist(nbinmax)
      character*80 string
      logical*4 ldebug/.false./

      if(nbin.gt.nbinmax)then
        write(*,"('ffind_mode error: number of bins exceeds', i6)")nbinmax
        call exit (1)
      endif

c     copy hist to shist
      do 7 i=1,nbin
    7 shist(i)=float(hist(i))

c     smooth three times
      do 6 k=1,3
      do 4 i=2,nbin-1
    4 thist(i)=(shist(i-1)+2.0*shist(i)+shist(i+1))/4.0
      do 5 i=2,nbin-1
    5 shist(i)=thist(i)
    6 continue

c     find maximum bin
      hist_max=0.
      do 3 i=1,nbin
      if(shist(i).gt.hist_max)then
        hist_max=shist(i)
        imax=i
      endif
    3 continue

c     type image of hitogram
      write(*,"('voxel value histogram')")
      do 11 l=15,0,-1
      string=' '
      do 12 k=1,80
      i=nint(float(k-1)*float(nbin)/80.)+1
   12 if(shist(i).ge.float(l)*hist_max/15.)string(k:k)='x'
   11 write(*,"(a80)")string
      string=' '
      k=nint(float(imax-1)*80./float(nbin))+1
      string(k:k)='|'
      write(*,"(a80)")string

c     find analytic maximum by parabolic interpolation
      d=(binval(nbin)-binval(1))/float(nbin-1)
      amode=binval(imax+1)-d*(0.5+(shist(imax+1)-shist(imax))/(shist(imax+1)-2.0*shist(imax)+shist(imax-1)))

c     write(*,"('voxel_value_range count')")
c     do 8 i=imax-14,imax+14,2
c   8 write(*,"(f7.1,' to',f7.1,i7)")binval(i),binval(i+2),hist(i)+hist(i+1)

      if(ldebug)then
        do i=1,nbin
          write(*,"(f10.2,i10,f10.2)")binval(i),hist(i),shist(i)
        enddo
        write(*,"('mode value= ', f10.4)")amode
      endif

      return
      end


      
