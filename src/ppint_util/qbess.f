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
      subroutine qbess( alpha,      apwr,       aterm1,     aterm2,
     &      binom,      bpref,      bpwr,       bterm1,     dfac,
     &      l,          lambu,      lmahi,      lmbhi,      ltot1,
     &      nu,         prd,        qsum,       rka,        rkb,
     &      ssi )
c
c  compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
c  using the bessel function formula for
c  lama=max(l,nu) to lmahi, lamb=max(l,nu) to lmbhi, n=lama+lamb-l-l
c
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1=1.0d0, a2=2.0d0, a4=4.0d0)
      dimension apwr(*), aterm1(*), aterm2(*), binom(*), bpref(*),
     &  bpwr(*), bterm1(*), dfac(*), qsum(ltot1,lambu,*), ssi(*)
c
c     nu=l+1-npi/2
c     # bessel function formula applies to all npi=2 cases, no npi=1
c     # cases, and some npi=0 cases.
      num1=nu-1
      lmlo=max(l,nu)
      lmaphi=lmahi-num1
      lmbphi=lmbhi-num1
      fcta=rka/(alpha+alpha)
      fct=rka*fcta
      fctb=rkb/(alpha+alpha)
      bpref(1)=fctb**num1
      do 24 lambp=2,lmbphi
        bpref(lambp)=fctb*bpref(lambp-1)
   24 continue
      apwr(1)=a1
      apwr(2)=fct
      do 28 lambp=3,lmbphi
        apwr(lambp)=fct*apwr(lambp-1)
   28 continue
      fct=rkb*fctb
      bpwr(1)=a1
      bpwr(2)=fct
      do 32 lamap=3,lmaphi
        bpwr(lamap)=fct*bpwr(lamap-1)
   32 continue
      lmihi=lmaphi+lmbphi+(nu-2)
      call ssibfn(lmihi-1,rka*fctb,ssi)
      do 36 lami=1,lmihi
        ssi(lami)=ssi(lami)/dfac(lami+lami-1)
   36 continue
      lmplo=lmlo-num1
      fctra=(alpha+alpha)**(nu-2)*fcta**(lmlo-1)*prd/sqrt(a4*alpha)*
     &  ((dfac(2*(2*lmplo+num1)-1)/dfac(2*(lmplo+num1)+1))*
     &  dfac(2*num1+1))/dfac(2*(lmplo+num1)+1)
      fctran=(2*(2*lmplo+num1)-1)
      fctrad=(2*(lmplo+num1)+1)
      do 64 lamap=lmplo,lmaphi
        fctru=a1
        fctrun=(nu+num1)
        fctrud=(2*(lamap+num1)+1)
        do 40 iu=1,lmbphi
          aterm1(iu)=fctru*apwr(iu)
          fctru=fctru*fctrun/fctrud
          fctrun=fctrun+a2
          fctrud=fctrud+a2
   40   continue
        do 44 it=1,lamap
          bterm1(it)=binom((lamap*(lamap-1))/2+it)*bpwr(it)
   44   continue
        fctrb=fctra
        fctrbn=(2*(lamap+lmplo+num1)-1)
        fctrbd=(2*(lmplo+num1)+1)
        do 60 lambp=lmplo,lmbphi
          n=((2*(nu-l)-1)+lamap)+lambp
          do 48 iu=1,lambp
            aterm2(iu)=binom((lambp*(lambp-1))/2+iu)*aterm1(iu)
   48     continue
          sum=a0
          fctrt=a1
          fctrtn=(nu+num1)
          fctrtd=(2*(lambp+num1)+1)
          do 56 it=1,lamap
            do 52 iu=1,lambp
              sum=sum+aterm2(iu)*(fctrt*bterm1(it))*
     &          ssi((it+(nu-2))+iu)
   52       continue
            fctrt=fctrt*fctrtn/fctrtd
            fctrtn=fctrtn+a2
            fctrtd=fctrtd+a2
   56     continue
          qsum(n,lambp+num1,lamap+num1)=qsum(n,lambp+num1,lamap+num1)+
     &      fctrb * bpref(lambp) * sum
          fctrb=fctrb*fctrbn/fctrbd
          fctrbn=fctrbn+a2
          fctrbd=fctrbd+a2
   60   continue
        fctra=fcta*fctra*fctran/fctrad
        fctran=fctran+a2
        fctrad=fctrad+a2
   64 continue
      return
      end
