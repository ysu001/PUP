      subroutine get_flipz(imag,np,nz,n_z,v_z,flipz)
      real*4 imag(np,nz),v_z(nz)
      integer*4 n_z(nz)
      logical*4 flipz

      e=0.
      o=0.
      ne=0
      no=0
      do iz=1,nz-1,2
        e=e+v_z(iz+1)*float(n_z(iz+1))
        ne=ne+n_z(iz+1)
        o=o+v_z(iz)*float(n_z(iz))
        no=no+n_z(iz)
      enddo
      e=e/float(ne)
      o=o/float(no)
      flipz=e.gt.o
      write(*,"('mean even/odd',f10.4,'/',f10.4,'  flipz=',l1)")e,o,flipz

      return
      end

      subroutine fdebandh(flipz,imag,np,nz,nv,v_z)
      logical*4 flipz
      real*4 imag(np,nz,nv),v_z(nz)

      iz0=2
      if(flipz)iz0=3
      do 9 iv=1,nv
      do iz=iz0,nz-1,2
        a=0.5*(v_z(iz-1)+v_z(iz+1))/v_z(iz)
        do ip=1,np
          imag(ip,iz,iv)=imag(ip,iz,iv)*a
        enddo
      enddo

      if(flipz)then
        a=v_z(2)/v_z(1)
        do ip=1,np
          imag(ip,1,iv)=imag(ip,1,iv)*a
        enddo
      endif

      if(flipz.eqv.mod(nz,2).eq.1)then
        a=v_z(nz-1)/v_z(nz)
        do ip=1,np
          imag(ip,nz,iv)=imag(ip,nz,iv)*a
        enddo
      endif
    9 continue

      return
      end

      subroutine fbande(flipz,imag,np,nz,mask,gamma)
      logical*4 flipz
      real*4 imag(np,nz)
      integer*4 mask(np,nz)

      if(flipz)then
        iz0=nz-1
        iz1=2
        jz=-2
      else
        iz0=2
        iz1=nz-1
        jz=2
      endif

      gamma=0.
      do 9 iter=1,7
      err=0.
      g=0.
      dg=0.
      npix=0
      kz=1
      do 8 iz=iz0,iz1,jz
      az=float(kz)
      e=exp(az*gamma)
      do ip=1,np
        if((mask(ip,iz-1).gt.0).and.(mask(ip,iz).gt.0).and.(mask(ip,iz+1).gt.0))then
          vk=imag(ip,iz)
          vb=0.5*(imag(ip,iz-1)+imag(ip,iz+1))
          err=err+(e*vk-vb)**2
          g=g+az*vk*e*(e*vk-vb)
          dg=dg+az*az*vk*e*(2.*e*vk-vb)
          npix=npix+1
        endif
      enddo
    8 kz=kz+2
      write(*, "('gamma =',f10.6)")gamma
      err=err/float(npix)
      g=2.*g/float(npix)
      dg=2.*dg/float(npix)
      write(*, "('rmserr=',f10.4,'  g=',e12.4,'  dg=',e12.4)")sqrt(err),g,dg
    9 gamma=gamma-g/dg

      return
      end

      subroutine fdebande(flipz,stack,np,nz,nv,gamma)
      logical*4 flipz
      real*4 stack(np,nz,nv)

      if(flipz)then
        iz0=nz-1
        iz1=2
        jz=-2
      else
        iz0=2
        iz1=nz-1
        jz=2
      endif

      write(*,"('Debanding: even slices multiplied by exp[',f10.6,'*islice]')")gamma
      do 9 iv=1,nv
      kz=1
      do 8 iz=iz0,iz1,jz
      az=float(kz)
      e=exp(az*gamma)
      do ip=1,np
        stack(ip,iz,iv)=e*stack(ip,iz,iv)
      enddo
    8 kz=kz+2
    9 continue

      return
      end

      subroutine fbandg(flipz,imag,np,nz,mask,beta)
      logical*4 flipz
      real*4 imag(np,nz)
      integer*4 mask(np,nz)

      if(flipz)then
        iz0=nz-1
        iz1=2
        jz=-2
      else
        iz0=2
        iz1=nz-1
        jz=2
      endif

      top=0.
      bot=0.
      kz=1
      do 8 iz=iz0,iz1,jz
      az=float(kz)
      do ip=1,np
        if((mask(ip,iz-1).gt.0).and.(mask(ip,iz).gt.0).and.(mask(ip,iz+1).gt.0))then
          vk=imag(ip,iz)
          vb=0.5*(imag(ip,iz-1)+imag(ip,iz+1))
          top=top+vk*(vb-vk)*az
          bot=bot+(az*vk)**2
        endif
      enddo
    8 kz=kz+2

      beta=top/bot
      write(*, "('beta =',f10.6)")beta
      return
      end

      subroutine fdebandg(flipz,stack,np,nz,nv,beta)
      logical*4 flipz
      real*4 stack(np,nz,nv)

      if(flipz)then
        iz0=nz-1
        iz1=2
        jz=-2
      else
        iz0=2
        iz1=nz-1
        jz=2
      endif

      write(*,"('Debanding: even slices multiplied by 1.0 + ',f10.6,'*islice')")beta
      do 9 iv=1,nv
      kz=1
      do 8 iz=iz0,iz1,jz
      az=float(kz)
      f=1.+az*beta
      do ip=1,np
        stack(ip,iz,iv)=f*stack(ip,iz,iv)
      enddo
    8 kz=kz+2
    9 continue

      return
      end
