cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009
c Washington University School od Medicine, Mallinckrodt Institute of Radiology.
c All rights reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of Washington University. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Header: /data/petsun4/data1/src_solaris/librms/RCS/dmatopr.f,v 1.1 2009/01/20 01:01:16 avi Exp $
c $Log: dmatopr.f,v $
c Revision 1.1  2009/01/20  01:01:16  avi
c Initial revision
c
      subroutine dmatmul(a,b,c,n)
      real*8 a(n,n),b(n,n),c(n,n)
      do 2 i=1,n
      do 2 j=1,n
      c(i,j)=0.
      do 2 k=1,n
    2 c(i,j)=c(i,j)+a(i,k)*b(k,j)
      return
      end
      subroutine dmatmulg(a,l,m,b,n,c)
      real*8 a(l,m),b(m,n),c(l,n)
      do 1 i=1,l
      do 1 j=1,n
      c(i,j)=0.
      do 1 k=1,m
    1 c(i,j)=c(i,j)+a(i,k)*b(k,j)
      return
      end
      subroutine dmatcop(a,b,n)
      real*8 a(n,n),b(n,n)
      do 1 i=1,n
      do 1 j=1,n
    1 b(i,j)=a(i,j)
      return
      end
      subroutine dtranspos(a,b,n)
      real*8 a(n,n),b(n,n)
      do 1 i=1,n
      do 1 j=1,n
    1 b(j,i)=a(i,j)
      return
      end
