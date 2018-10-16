c$Header: /data/petsun4/data1/src_solaris/nonlin/RCS/imgval0.f,v 1.1 2004/10/14 23:59:49 rsachs Exp $
c$Log: imgval0.f,v $
c Revision 1.1  2004/10/14  23:59:49  rsachs
c Initial revision
c      
      subroutine imgval0(imgt,nx,ny,nz,voxdim,x,val)
c     returns in val the interpolated value of imgt at locus x.
c     val = 0 if imgt is undefined at x.

      real*4 imgt(0:nx-1,0:ny-1,0:nz-1),voxdim(3),x(3),val
      logical*4 ldefined,lok
      real*4 wx(3),xf(3)
      integer*4 ix(3),mx(3)

      mx(1)=nx
      mx(2)=ny
      mx(3)=nz

      ldefined=.true.
      do 1 k=1,3
c1      center(k)=sign(center(k),mmppix(k))
c1      xf(k)=(center(k)+x(k))/mmppix(k)
      xf(k)=x(k)/voxdim(k)            ! the fractional index 
      ix(k)=nint(xf(k)-0.5)    
      wx(k)=xf(k)-float(ix(k))
      if((ix(k).eq.-1.and.wx(k).gt.0.999).or.mx(k).eq.1)then
        ix(k)=0       !  ix(k)=1
        wx(k)=0.
      endif
c1    if(mx(k).gt.1.and.ix(k).eq.mx(k).and.wx(k).lt.0.001)then
      if(mx(k).gt.1.and.ix(k).eq.(mx(k)-1).and.wx(k).lt.0.001)then
        ix(k)=mx(k)-2         !    ix(k)=mx(k)-1
        wx(k)=1.
      endif
c1    1 ldefined=ldefined.and.((mx(k).eq.1).or.(ix(k).gt.0).and.(ix(k).lt.mx(k)))

    1  ldefined=ldefined.and.((mx(k).eq.1).or.(ix(k).ge.0).and.(ix(k).lt.mx(k)-1))
      if(.not.ldefined)then
        val=0.
        return
      endif

      lok=(imgt(ix(1)+0,ix(2)+0,ix(3)+0).ge.0.0.or.imgt(ix(1)+0,ix(2)+0,ix(3)+0).le.0.0).and.
     &    (imgt(ix(1)+1,ix(2)+0,ix(3)+0).ge.0.0.or.imgt(ix(1)+1,ix(2)+0,ix(3)+0).le.0.0).and.
     &    (imgt(ix(1)+0,ix(2)+1,ix(3)+0).ge.0.0.or.imgt(ix(1)+0,ix(2)+1,ix(3)+0).le.0.0).and.
     &    (imgt(ix(1)+1,ix(2)+1,ix(3)+0).ge.0.0.or.imgt(ix(1)+1,ix(2)+1,ix(3)+0).le.0.0).and.
     &    (imgt(ix(1)+0,ix(2)+0,ix(3)+1).ge.0.0.or.imgt(ix(1)+0,ix(2)+0,ix(3)+1).le.0.0).and.
     &    (imgt(ix(1)+1,ix(2)+0,ix(3)+1).ge.0.0.or.imgt(ix(1)+1,ix(2)+0,ix(3)+1).le.0.0).and.
     &    (imgt(ix(1)+0,ix(2)+1,ix(3)+1).ge.0.0.or.imgt(ix(1)+0,ix(2)+1,ix(3)+1).le.0.0).and.
     &    (imgt(ix(1)+1,ix(2)+1,ix(3)+1).ge.0.0.or.imgt(ix(1)+1,ix(2)+1,ix(3)+1).le.0.0)
      if(.not.lok)then
        val=0.
        return
      endif

      val=
     & +(1.-wx(3))*((1.-wx(2))*(
     &               (1.-wx(1))*imgt(ix(1)+0,ix(2)+0,ix(3)+0)
     &                  +wx(1) *imgt(ix(1)+1,ix(2)+0,ix(3)+0))
     &               +wx(2)*(
     &               (1.-wx(1))*imgt(ix(1)+0,ix(2)+1,ix(3)+0)
     &                  +wx(1) *imgt(ix(1)+1,ix(2)+1,ix(3)+0)))
       if(mx(3).gt.1) val=val
     &      +wx(3)*((1.-wx(2))*(
     &               (1.-wx(1))*imgt(ix(1)+0,ix(2)+0,ix(3)+1)
     &                  +wx(1) *imgt(ix(1)+1,ix(2)+0,ix(3)+1))
     &               +wx(2)*(
     &               (1.-wx(1))*imgt(ix(1)+0,ix(2)+1,ix(3)+1)
     &                  +wx(1) *imgt(ix(1)+1,ix(2)+1,ix(3)+1)))
      return
      end

