C$Id: fdeband.f,v 1.2 2009/02/22 02:17:30 avi Exp $
C$Log: fdeband.f,v $
c Revision 1.2  2009/02/22  02:17:30  avi
c mask integer*2 -> integer*4
c
c Revision 1.1  1997/04/28  00:46:06  yang
c Initial revision
c
      subroutine fsquash(stack,nvox,nvol,imag)
      real*4 stack(nvox,nvol),imag(nvox)

      do 2 i=1,nvox
      imag(i)=0.
      do 1 j=1,nvol
    1 imag(i)=imag(i)+stack(i,j)
    2 imag(i)=imag(i)/float(nvol)
      return
      end

      subroutine ffluct(stack,np,nz,nv,imag,mask,n_z,v_z)
      real*4 stack(np,nz,nv),imag(np,nz),v_z(nz)
      integer*4 mask(np,nz)
      integer*4 n_z(nz)
      real*4 dv_z(256)

      do iv=1,nv
        do iz=1,nz
          sum=0.
          do ip=1,np
            if(mask(ip,iz).gt.0)sum=sum+stack(ip,iz,iv)
          enddo
          dv_z(iz)=sum/float(n_z(iz))-v_z(iz)
        enddo
        write(*,"(i3,1x,19f6.2)")iv,(dv_z(iz),iz=1,nz)
      enddo
      return
      end

      subroutine fband(imag,np,nz,mask,n_z,v_z)
      real*4 imag(np,nz),v_z(nz)
      integer*4 mask(np,nz)
      integer*4 n_z(nz)
      logical*4 lodd

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
        write(*,"(i10,f10.4,i10)")iz,v_z(iz),n_z(iz)
      enddo

      sumtot=0.
      vartot=0.
      sumodd=0.
      sumeve=0.
      ntot=0
      neve=0
      nodd=0
      
      do iz=2,nz-1
      lodd=mod(iz,2).eq.1
        do ip=1,np
          if((mask(ip,iz-1).gt.0).and.(mask(ip,iz).gt.0).and.(mask(ip,iz+1).gt.0))then
            v=imag(ip,iz)
            sumtot=sumtot+v
            vartot=vartot+v**2
            ntot=ntot+1
            if(lodd)then
              sumodd=sumodd+v
              nodd=nodd+1
            else
              sumeve=sumeve+v
              neve=neve+1
            endif
          endif
        enddo
      enddo

      write(*,"('Npix  even, odd, total',3i10)")neve,nodd,ntot
      write(*,"('Means even, odd, total',3f10.4)")sumeve/float(neve),sumodd/float(nodd),sumtot/float(ntot)
      f=sumeve**2/float(neve)+sumodd**2/float(nodd)-sumtot**2/float(ntot)
      write(*,"('Fraction of variance attributable to banding = ',e10.4)")f/vartot

      return
      end

      subroutine fbandr(imag,np,nz,mask,alpha)
      real*4 imag(np,nz)
      integer*4 mask(np,nz)

      top=0.
      bot=0.
      do iz=2,nz-1
        s=(-1.)**mod(iz,2)
        do ip=1,np
          if((mask(ip,iz-1).gt.0).and.(mask(ip,iz).gt.0).and.(mask(ip,iz+1).gt.0))then
            vk=imag(ip,iz)
            vb=0.5*(imag(ip,iz-1)+imag(ip,iz+1))
            top=top+s*vk*(vb-vk)
            bot=bot+vk**2
          endif
        enddo
      enddo

      alpha=0.5*top/bot
c     write (*, "('alpha =',f10.6)",alpha
      return
      end

      subroutine fdebandr(stack,np,nz,nv,alpha)
      real*4 stack(np,nz,nv)

      fe=1.+alpha
      fo=1.-alpha
      write(*,"('Debanding: slice multipliers even =',f10.6,' odd =',f10.6)")fe,fo
      do iv=1,nv
        do iz=1,nz
          if(mod(iz,2).eq.0)then
            f=fe
          else
            f=fo
          endif
          do ip=1,np
            stack(ip,iz,iv)=f*stack(ip,iz,iv)
          enddo
        enddo
      enddo

      return
      end
