************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine hotred (nx,n,a,d,e,z)
      implicit real*8  (a-h,o-z)
      dimension a(nx*(nx+1)/2),z(nx,nx),d(nx),e(nx)
      if(n.le.2) go to 200
      ij=0
      do 10 i=1,n
      do 11 j=1,i
      ij=ij+1
      z(i,j)=a(ij)
   11 continue
   10 continue
      do 100 ip=2,n
      i=n-ip+2
      l=i-2
      f=z(i,i-1)
      g=0.0d0
      if(l.eq.0) go to 25
      do 20 k=1,l
      g=g+z(i,k)*z(i,k)
   20 continue
   25 h=g+f*f
      if(g.gt.1d-12) go to 30
      e(i)=f
      h=0.0d0
      go to 90
   30 l=l+1
      g=sqrt(h)
      if(f.ge.0.0d0) g=-g
      e(i)=g
      h=h-f*g
      z(i,i-1)=f-g
      f=0.0d0
      do 70 j=1,l
      z(j,i)=z(i,j)/h
      g=0.0d0
      do 40 k=1,j
      g=g+z(j,k)*z(i,k)
   40 continue
      jn=j+1
      if(l.lt.jn) go to 60
      do 50 k=jn,l
      g=g+z(k,j)*z(i,k)
   50 continue
   60 e(j)=g/h
      f=f+g*z(j,i)
   70 continue
      hh=f/(h+h)
      do 80 j=1,l
      f=z(i,j)
      g=e(j)-hh*f
      e(j)=g
      do 81 k=1,j
      z(j,k)=z(j,k)-f* e(k)-g*z(i,k)
   81 continue
   80 continue
   90 d(i)=h
  100 continue
      d(1)=z(1,1)
      z(1,1)=1.0d0
      e(1)=0.0d0
      do 150 i=2,n
      l=i-1
      if(d(i).eq.0.0d0) go to 130
      do 120 j=1,l
      g=0.0d0
      do 110 k=1,l
      g=g+z(i,k)*z(k,j)
  110 continue
      do 121 k=1,l
      z(k,j)=z(k,j)-g*z(k,i)
  121 continue
  120 continue
  130 d(i)=z(i,i)
      z(i,i)=1.0d0
      do 151 j=1,l
      z(i,j)=0.0d0
      z(j,i)=0.0d0
  151 continue
  150 continue
      return
  200 go to (210,220),n
  210 d(1)=a(1)
      e(1)=0.0d0
      z(1,1)=1.0d0
      return
  220 d(1)=a(1)
      d(2)=a(3)
      e(1)=0.0d0
      e(2)=a(2)
      z(1,1)=1.0d0
      z(2,2)=1.0d0
      z(1,2)=0.0d0
      z(2,1)=0.0d0
      return
      end
c
      subroutine qlcm(nx,n,d,e,z)
      implicit real*8 (a-h,o-z)
      dimension z(nx,nx),d(nx),e(nx)
      do 10 i=2,n
      e(i-1)=e(i)
   10 continue
      e(n)=0.0d0
      b=0.0d0
      f=0.0d0
      do 150 l=1,n
      j=0
      h=1d-12*(abs(d(l))+abs(e(l)))
      if(b.ge.h) go to 20
      b=h
   20 do 30 m=l,n
      if(abs(e(m)).le.b) go to 40
   30 continue
   40 if(m.eq.l) go to 140
   50 if(j.eq.nx+1) go to 200
      j=j+1
      g=d(l)
      p=(d(l+1)-g)/(2*e(l))
      r=sqrt(p*p+1.0d0)
      if(p.lt.0.0d0) go to 60
      pp=p+r
      go to 70
   60 pp=p-r
   70 d(l)=e(l)/pp
      h=g-d(l)
      if(l.eq.n)  go to 85
      ll=l+1
      do 80 i=ll,n
      d(i)=d(i)-h
   80 continue
   85 f=f+h
      p=d(m)
      c=1.0d0
      s=0.0d0
      mm=m-1
      do 110 kk=l,mm
      i=mm+l-kk
      g=c*e(i)
      h=c*p
      if(abs(p).lt.abs(e(i))) go to 90
        c=e(i)/p
        r=sqrt(c*c+1.0d0)
      e(i+1)=s*p*r
      s=c/r
      c=1.0d0/r
      go to 100
   90 c=p/e(i)
      r=sqrt(c*c+1.0d0)
      e(i+1)=s*e(i)*r
      s=1.0d0/r
      c=c/r
  100 p=c*d(i)-s*g
      d(i+1)=h+s*(c*g+s*d(i))
      do 111 k=1,n
      h=z(k,i+1)
      z(k,i+1)=s*z(k,i)+c*h
      z(k,i)=c*z(k,i)-s*h
  111 continue
  110 continue
      e(l)=s*p
      d(l)=c*p
      if(abs(e(l)).gt.b) go to 50
  140 d(l)=d(l)+f
  150 continue
      do 180 i=1,n
      k=i
      p=d(i)
      l=i+1
      if(l.gt.n) go to 165
      do 160 j=l,n
      if(d(j).ge.p) go to 160
      k=j
      p=d(j)
  160 continue
  165 if(k.eq.i) go to 180
      d(k)=d(i)
      d(i)=p
      do 170 j=1,n
      p=z(j,i)
      z(j,i)=z(j,k)
      z(j,k)=p
  170 continue
  180 continue
      return
  200 write(6,250)
  250 format(1x///5x,'***the subroutine qlcm is fail,so this',
     *       'computation must stop***')
#ifdef MOLPRO
#else
      call abend()
#endif
#ifdef _XIANEST_
#endif
!     call abend
!      stop
      end
