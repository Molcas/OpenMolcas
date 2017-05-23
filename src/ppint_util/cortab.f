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
      subroutine cortab_molcas(binom,dfac,eps,flmtx,hpt,
     &  hwt,lmf,lml,lmx,lmy,
     &  lmz,lmax,lmn1u,lproju,mc,mr,ndfac,zlm)
c     # tables for core potential and spin-orbit integrals.
      implicit real*8 (a-h,o-z)
      parameter (a0=0.0d0, a1=1.0d0, a3=3.0d0)
      dimension binom(*),dfac(*),flmtx(3,*),hpt(*),hwt(*),lmf(*),
     &  lml(*),lmx(*),lmy(*),lmz(*),mc(3,*),mr(3,*),zlm(*)
c
c     # compute gauss-hermite points and weights for c, z integrals.
      igh =1
      nn = 5
      do 120 i = 1, 3
        call hermit_molcas(nn,hpt(igh),hwt(igh),eps)
        igh = igh + nn
        nn = 2*nn
  120 continue
c     # compute double factorials.
      dfac(1)=a1
      dfac(2)=a1
      fi=a0
      do 130 i=1,ndfac-2
        fi=fi+a1
        dfac(i+2)=fi * dfac(i)
  130 continue
c     # compute binomial coefficients.
      inew=1
      binom(1)=a1
      do 150 j=1,lmn1u-1
        inew=inew+1
        binom(inew)=a1
        do 140 i=1,j-1
          inew=inew+1
          binom(inew)=(dble(j-i+1)*binom(inew-1))/dble(i)
  140   continue
        inew=inew+1
        binom(inew)=a1
  150 continue
c     # compute tables by recursion for real spherical harmonics.  they
c     # are indexed by l, m and sigma.  the sequence number of the
c     # harmonic with quantum numbers l, m and sigma is given by
c     #            l**2+2*m+1-sigma
c     # lmf(index) and lml(index) hold the positions of the first and
c     # last terms of the harmonic in the lmx, lmy, lmz, and zlm arrays.
c     # the harmonics with angular momentum l are generated from those
c     # with angular momenta l-1 and l-2.
c     # for m = 0,1,2,...,l-1, the recursion relation
c     z*Z(l-1,m,s) = sqrt(((l-m)*(l+m))/((2*l-1)*(2*l+1)))*Z(l,m,s)+
c                  sqrt(((l+m-1)*(l-m-1))/((2*l-3)*(2*l-1)))*Z(l-2,m,s)
c     # is used.
c     # for m = l, the recursion relation
c     x*Z(l-1,l-1,s)+(-1)**(1-s)*y*Z(l-1,l-1,1-s) =
c                  sqrt((2*l))/((2*l+1)))*Z(l,l,s)
c     # is used.
c     # l=0
      lmf(1) = 1
      lml(1) = 1
      lmx(1) = 0
      lmy(1) = 0
      lmz(1) = 0
      zlm(1) = a1
c     # l=1
      lmf(2) = 2
      lml(2) = 2
      lmx(2) = 0
      lmy(2) = 0
      lmz(2) = 1
      zlm(2) = sqrt(a3)
      lmf(3) = 3
      lml(3) = 3
      lmx(3) = 0
      lmy(3) = 1
      lmz(3) = 0
      zlm(3) = zlm(2)
      lmf(4) = 4
      lml(4) = 4
      lmx(4) = 1
      lmy(4) = 0
      lmz(4) = 0
      zlm(4) = zlm(2)
      nterm=4
      do 270 lang=2,lmax
        do 240 mang=0,lang-1
          anum = dble((2*lang-1)*(2*lang+1))
          aden = dble((lang-mang)*(lang+mang))
          coef1 = sqrt(anum/aden)
          anum = dble((lang+mang-1)*(lang-mang-1)*(2*lang+1))
          aden = dble(2*lang-3)*aden
          coef2 = sqrt(anum/aden)
          nsigma=min(1,mang)
          do 230 isigma=nsigma,0,-1
            indexh=lang**2+2*mang+1-isigma
            lone=lang-1
            ltwo=lang-2
            ione=lone**2+2*mang+1-isigma
            itwo=ltwo**2+2*mang+1-isigma
            lmf(indexh)=lml(indexh-1)+1
            lml(indexh)=lml(indexh-1)
            nxy=(mang-isigma+2)/2
            iu=lmf(ione)+nxy-1
            do 200 i=lmf(ione),iu
              lml(indexh)=lml(indexh)+1
              j=lml(indexh)
              lmx(j)=lmx(i)
              lmy(j)=lmy(i)
              lmz(j)=lmz(i)+1
              zlm(j)=zlm(i)*coef1
              nterm=nterm+1
  200       end do
            if(ltwo.ge.mang) then
              il=iu+1
              do 210 i=il,lml(ione)
                lml(indexh)=lml(indexh)+1
                j=lml(indexh)
                k=lmf(itwo)+i-il
                lmx(j)=lmx(k)
                lmy(j)=lmy(k)
                lmz(j)=lmz(k)
                zlm(j)=zlm(i)*coef1-zlm(k)*coef2
                nterm=nterm+1
  210         enddo
              il=lml(itwo)-nxy+1
              if(mod(lang-mang,2).eq.0) then
                do 220 i=il,lml(itwo)
                  lml(indexh)=lml(indexh)+1
                  j=lml(indexh)
                  lmx(j)=lmx(i)
                  lmy(j)=lmy(i)
                  lmz(j)=lmz(i)
                  zlm(j)=-zlm(i)*coef2
                  nterm=nterm+1
  220           end do
              endif
            endif
  230     end do
  240   end do
        anum = dble(2*lang+1)
        aden = dble(2*lang)
        coef = sqrt(anum/aden)
        mang=lang
        isigma=1
        indexh=lang**2+2*mang+1-isigma
        lmf(indexh)=lml(indexh-1)+1
        lml(indexh)=lml(indexh-1)
c       # isig:  index of the harmonic (l-1),(m-1),sigma
c       # isigm: index of the harmonic (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        k=lmf(isigm)
        do 250 i=lmf(isig),lml(isig)
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(i)+1
          lmy(j)=lmy(i)
          lmz(j)=lmz(i)
          zlm(j)=(zlm(i)+zlm(k))*coef
          k=k+1
          nterm=nterm+1
  250   enddo
        if(mod(mang,2).eq.1) then
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(k)
          lmy(j)=lmy(k)+1
          lmz(j)=lmz(k)
          zlm(j)=zlm(k)*coef
          nterm=nterm+1
        endif
        isigma=0
        indexh=lang**2+2*mang+1-isigma
c       # isig:  index of the harmonic (l-1),(m-1),sigma
c       # isigm: index of the harmonuc (l-1),(m-1),(1-sigma)
        isig=(lang-1)**2+2*(mang-1)+1-isigma
        isigm=(lang-1)**2+2*(mang-1)+isigma
        lmf(indexh)=lml(indexh-1)+1
        lml(indexh)=lmf(indexh)
        j=lml(indexh)
        i=lmf(isig)
        lmx(j)=lmx(i)+1
        lmy(j)=lmy(i)
        lmz(j)=lmz(i)
        zlm(j)=zlm(i)*coef
        nterm=nterm+1
        k=lmf(isigm)
        do 260 i=lmf(isig)+1,lml(isig)
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          lmx(j)=lmx(i)+1
          lmy(j)=lmy(i)
          lmz(j)=lmz(i)
          zlm(j)=(zlm(i)-zlm(k))*coef
          k=k+1
          nterm=nterm+1
  260   enddo
        if(mod(mang,2).eq.0) then
          lml(indexh)=lml(indexh)+1
          j=lml(indexh)
          k=lml(isigm)
          lmx(j)=lmx(k)
          lmy(j)=lmy(k)+1
          lmz(j)=lmz(k)
          zlm(j)=-zlm(k)*coef
          nterm=nterm+1
        endif
  270 end do
      ixy = 0
      iz = 0
      do 300 lang=1,lproju
        do 290 mang=0,lang-1
          nsigma=min(1,mang)
          ndelta=max(0,1-mang)
          anum = dble((lang-mang)*(lang+mang+1))
          aden = dble(2*(2-ndelta))
          coef=sqrt(anum/aden)
          do 280 isigma=nsigma,0,-1
            isign=2*isigma-1
            ixy = ixy+1
            flmtx(1,ixy) = dble((isign))*coef
            flmtx(2,ixy) = coef
            if(mang.ne.0) then
              iz=iz+1
              flmtx(3,iz) = -dble(mang*isigma)
            endif
  280     enddo
  290   enddo
        iz=iz+1
        flmtx(3,iz) = -dble(lang)
  300 enddo
c     # column and row indices for angular momentum matrix elements.
      iadd = 1
      do 310 i=1,2*lproju-1
        mc(1,i) = i
        mc(2,i) = i
        mc(3,i) = i+1
        mr(1,i) = i+iadd
        mr(2,i) = i+2
        mr(3,i) = i+2
        iadd = 4-iadd
  310 continue
      return
      end
