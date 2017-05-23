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
      subroutine ssibfn(nmax,x,ssi)
c
c  scaled spherical i bessel functions
c
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1=1.0d0, a3s2=1.5d0, a2=2.0d0, a3=3.0d0,
     1  a20=20.0d0)
      dimension ssi(*)
c
      x2=x*x
      xmin=(abs(3*nmax-1))
      if(x.gt.xmin) go to 5
      n=nmax
      f2np1=(n+n+1)
      f2kp3=f2np1
      pkm1=a0
      pk=a1
      qkm1=a1
      qk=a1
      aprod=a1
    1 f2kp1=f2kp3
      f2kp3=f2kp3+a2
      ak=x2/(f2kp1*f2kp3)
      akpkm2=ak*pkm1
      pkm1=pk
      pk=pkm1+akpkm2
      akqkm2=ak*qkm1
      qkm1=qk
      qk=qkm1+akqkm2
      aprod=ak*aprod
      if(((pk*qkm1)+aprod).ne.(pk*qkm1)) go to 1
      ssi(n+1)=pk/qk
    2 if(n.eq.0) go to 3
      n=n-1
      f2np3=f2np1
      f2np1=f2np1-a2
      ssi(n+1)=(f2np1*f2np3)/((f2np1*f2np3)+x2*ssi(n+2))
      go to 2
    3 ssi(1)=ssi(1)/(a1+x*ssi(1))
      do 4 n=1,nmax
        ssi(n+1)=ssi(n+1)*ssi(n)
    4 continue
      return
c
    5 if(x.ge.a20) then
        ex=a0
      else
        ex=exp(-(x+x))
      endif
      ssi(1)=(a1-ex)/(x+x)
      if(nmax.eq.0) return
      ssi(2)=a3s2*(a1+ex+(ex-a1)/x)/x2
      if(nmax.eq.1) return
      f2np1=a3
      do 8 n=2,nmax
        f2nm1=f2np1
        f2np1=f2np1+a2
        ssi(n+1)=(ssi(n-1)-ssi(n))*(f2nm1*f2np1)/x2
    8 continue
      return
      end
