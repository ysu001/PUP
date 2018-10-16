c$Header: /data/petsun4/data1/src_solaris/interp_4dfp/RCS/pbutt1dc.f,v 1.1 2012/04/16 03:05:45 avi Exp $
c$Log: pbutt1dc.f,v $
c Revision 1.1  2012/04/16  03:05:45  avi
c Initial revision
c
      subroutine pbutt1dc_rcs
      write (*,"('$Id: pbutt1dc.f,v 1.1 2012/04/16 03:05:45 avi Exp $')")
      return
      end

      subroutine pbutt1dcs(data,n0,n,delta,fhalf_lo,iorder_lo,fhalf_hi,iorder_hi,a,b)
c     version of pbutt1dc with scratch arrays a and b passed by pointer
      real*4 data(n)
      real*4 a(n/2+1),b(n/2+1)

      if(mod(n,2).ne.0)then
        write(*,"('pbutt1dcs: illegal odd input array length',i6)")n
        call exit(-1)
      endif
      if(iorder_lo.lt.0.or.iorder_hi.lt.0)then
        write(*,"('pbutt1dcs: negative Butterworth filter orders not allowed')")
        call exit(-1)
      endif

      i=1
      do 21 k=1,n,2
      a(i)=data(k)
      b(i)=data(k+1)
   21 i=i+1

      call FFT  (a,b,1,n/2,1,-1)
      call REALT(a,b,1,n/2,1,-1)
      do 31 i=1,n/2+1
      f=float(i-1)/(float(n0)*delta)
      if(iorder_lo.gt.0)then
        r_lo=(f/fhalf_lo)**(2*iorder_lo)
        factor_lo=sqrt(r_lo/(1.0+r_lo))
      else
        factor_lo=1.0
      endif
      if(iorder_hi.gt.0)then
        r_hi=(f/fhalf_hi)**(2*iorder_hi)
        factor_hi=sqrt(1.0/(1.0+r_hi))
      else
        factor_hi=1.0
      endif
      a(i)=factor_lo*factor_hi*a(i)
   31 b(i)=factor_lo*factor_hi*b(i)
      call REALT(a,b,1,n/2,1,+1)
      call FFT  (a,b,1,n/2,1,+1)

      i=1
      do 41 k=1,n,2
      data(k)  =a(i)
      data(k+1)=b(i)
   41 i=i+1

      return
      end

      subroutine pbutt1dca(data,n0,n,delta,fhalf_lo,iorder_lo,fhalf_hi,iorder_hi)
c     version of pbutt1dc that allocates buffers on each call
      real*4 data(n)
      real*4 a(1),b(1)
      pointer (pa,a),(pb,b)

      if(mod(n,2).ne.0)then
        write(*,"('pbutt1dca: illegal odd input array length',i6)")n
        call exit(-1)
      endif
      if(iorder_lo.lt.0.or.iorder_hi.lt.0)then
        write(*,"('pbutt1dca: negative Butterworth filter orders not allowed')")
        call exit(-1)
      endif

      pa=malloc(4*(n/2+1))
      pb=malloc(4*(n/2+1))
      if(pa.eq.0.or.pb.eq.0)then
        write(*,"('pbutt1dca: memory allocation error')")
        call exit(-1)
      endif

      i=1
      do 21 k=1,n,2
      a(i)=data(k)
      b(i)=data(k+1)
   21 i=i+1

      call FFT  (a,b,1,n/2,1,-1)
      call REALT(a,b,1,n/2,1,-1)
      do 31 i=1,n/2+1
      f=float(i-1)/(float(n0)*delta)
      if(iorder_lo.gt.0)then
        r_lo=(f/fhalf_lo)**(2*iorder_lo)
        factor_lo=sqrt(r_lo/(1.0+r_lo))
      else
        factor_lo=1.0
      endif
      if(iorder_hi.gt.0)then
        r_hi=(f/fhalf_hi)**(2*iorder_hi)
        factor_hi=sqrt(1.0/(1.0+r_hi))
      else
        factor_hi=1.0
      endif
      a(i)=factor_lo*factor_hi*a(i)
   31 b(i)=factor_lo*factor_hi*b(i)
      call REALT(a,b,1,n/2,1,+1)
      call FFT  (a,b,1,n/2,1,+1)

      i=1
      do 41 k=1,n,2
      data(k)  =a(i)
      data(k+1)=b(i)
   41 i=i+1

      call free(pa)
      call free(pb)
      return
      end
