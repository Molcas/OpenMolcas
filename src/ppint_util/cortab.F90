!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine cortab(binom,dfac,eps,lmf,lml,lmx,lmy,lmz,lmax,lmn1u,ndfac,zlm)
! Tables for core potential and spin-orbit integrals.

use ppint_arrays, only: hpt, hwt
use stdalloc, only: mma_allocate
use Constants, only: Zero, One, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lmax, lmn1u, ndfac
#define _LSIZE_ ((lmax*(lmax+2)*(lmax+4)/3)*(lmax+3)+(lmax+2)**2*(lmax+4))/16
real(kind=wp), intent(out) :: binom(lmn1u*(lmn1u+1)/2), dfac(ndfac), &
                              zlm(_LSIZE_)
real(kind=wp), intent(in) :: eps
integer(kind=iwp), intent(out) :: lmf((lmax+1)**2), lml((lmax+1)**2), &
                                  lmx(_LSIZE_), &
                                  lmy(_LSIZE_), &
                                  lmz(_LSIZE_)
integer(kind=iwp) :: i, igh, il, indexh, inew, ione, isig, isigm, isigma, itwo, iu, j, k, lang, lone, ltwo, mang, nn, nsigma, &
                     nterm, nxy
real(kind=wp) :: aden, anum, coef, coef1, coef2, fi

! Compute Gauss-Hermite points and weights for c, z integrals.
igh = 1
nn = 5
call mma_allocate(hpt,nn*(2**3-1),label='hpt')
call mma_allocate(hwt,nn*(2**3-1),label='hwt')
do i=1,3
  call hermit(nn,hpt(igh),hwt(igh),eps)
  igh = igh+nn
  nn = 2*nn
end do
! Compute double factorials.
dfac(1) = One
dfac(2) = One
fi = Zero
do i=1,ndfac-2
  fi = fi+One
  dfac(i+2) = fi*dfac(i)
end do
! Compute binomial coefficients.
inew = 1
binom(1) = One
do j=1,lmn1u-1
  inew = inew+1
  binom(inew) = One
  do i=1,j-1
    inew = inew+1
    binom(inew) = real(j-i+1,kind=wp)*binom(inew-1)/real(i,kind=wp)
  end do
  inew = inew+1
  binom(inew) = One
end do
! Compute tables by recursion for real spherical harmonics. They
! are indexed by l, m and sigma. The sequence number of the
! harmonic with quantum numbers l, m and sigma is given by
!     l**2+2*m+1-sigma
! lmf(index) and lml(index) hold the positions of the first and
! last terms of the harmonic in the lmx, lmy, lmz, and zlm arrays.
! The harmonics with angular momentum l are generated from those
! with angular momenta l-1 and l-2.
! For m = 0,1,2,...,l-1, the recursion relation
!     z*Z(l-1,m,s) = sqrt(((l-m)*(l+m))/((2*l-1)*(2*l+1)))*Z(l,m,s)+
!                    sqrt(((l+m-1)*(l-m-1))/((2*l-3)*(2*l-1)))*Z(l-2,m,s)
! is used.
! For m = l, the recursion relation
!     x*Z(l-1,l-1,s)+(-1)**(1-s)*y*Z(l-1,l-1,1-s) = sqrt((2*l))/((2*l+1)))*Z(l,l,s)
! is used.

! l=0
lmf(1) = 1
lml(1) = 1
lmx(1) = 0
lmy(1) = 0
lmz(1) = 0
zlm(1) = One
! l=1
lmf(2) = 2
lml(2) = 2
lmx(2) = 0
lmy(2) = 0
lmz(2) = 1
zlm(2) = sqrt(Three)
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
nterm = 4
do lang=2,lmax
  do mang=0,lang-1
    anum = real((2*lang-1)*(2*lang+1),kind=wp)
    aden = real((lang-mang)*(lang+mang),kind=wp)
    coef1 = sqrt(anum/aden)
    anum = real((lang+mang-1)*(lang-mang-1)*(2*lang+1),kind=wp)
    aden = real(2*lang-3,kind=wp)*aden
    coef2 = sqrt(anum/aden)
    nsigma = min(1,mang)
    do isigma=nsigma,0,-1
      indexh = lang**2+2*mang+1-isigma
      lone = lang-1
      ltwo = lang-2
      ione = lone**2+2*mang+1-isigma
      itwo = ltwo**2+2*mang+1-isigma
      lmf(indexh) = lml(indexh-1)+1
      lml(indexh) = lml(indexh-1)
      nxy = (mang-isigma+2)/2
      iu = lmf(ione)+nxy-1
      do i=lmf(ione),iu
        lml(indexh) = lml(indexh)+1
        j = lml(indexh)
        lmx(j) = lmx(i)
        lmy(j) = lmy(i)
        lmz(j) = lmz(i)+1
        zlm(j) = zlm(i)*coef1
        nterm = nterm+1
      end do
      if (ltwo >= mang) then
        il = iu+1
        do i=il,lml(ione)
          lml(indexh) = lml(indexh)+1
          j = lml(indexh)
          k = lmf(itwo)+i-il
          lmx(j) = lmx(k)
          lmy(j) = lmy(k)
          lmz(j) = lmz(k)
          zlm(j) = zlm(i)*coef1-zlm(k)*coef2
          nterm = nterm+1
        end do
        il = lml(itwo)-nxy+1
        if (mod(lang-mang,2) == 0) then
          do i=il,lml(itwo)
            lml(indexh) = lml(indexh)+1
            j = lml(indexh)
            lmx(j) = lmx(i)
            lmy(j) = lmy(i)
            lmz(j) = lmz(i)
            zlm(j) = -zlm(i)*coef2
            nterm = nterm+1
          end do
        end if
      end if
    end do
  end do
  anum = real(2*lang+1,kind=wp)
  aden = real(2*lang,kind=wp)
  coef = sqrt(anum/aden)
  mang = lang
  isigma = 1
  indexh = lang**2+2*mang+1-isigma
  lmf(indexh) = lml(indexh-1)+1
  lml(indexh) = lml(indexh-1)
  ! isig:  index of the harmonic (l-1),(m-1),sigma
  ! isigm: index of the harmonic (l-1),(m-1),(1-sigma)
  isig = (lang-1)**2+2*(mang-1)+1-isigma
  isigm = (lang-1)**2+2*(mang-1)+isigma
  k = lmf(isigm)
  do i=lmf(isig),lml(isig)
    lml(indexh) = lml(indexh)+1
    j = lml(indexh)
    lmx(j) = lmx(i)+1
    lmy(j) = lmy(i)
    lmz(j) = lmz(i)
    zlm(j) = (zlm(i)+zlm(k))*coef
    k = k+1
    nterm = nterm+1
  end do
  if (mod(mang,2) == 1) then
    lml(indexh) = lml(indexh)+1
    j = lml(indexh)
    lmx(j) = lmx(k)
    lmy(j) = lmy(k)+1
    lmz(j) = lmz(k)
    zlm(j) = zlm(k)*coef
    nterm = nterm+1
  end if
  isigma = 0
  indexh = lang**2+2*mang+1-isigma
  ! isig:  index of the harmonic (l-1),(m-1),sigma
  ! isigm: index of the harmonic (l-1),(m-1),(1-sigma)
  isig = (lang-1)**2+2*(mang-1)+1-isigma
  isigm = (lang-1)**2+2*(mang-1)+isigma
  lmf(indexh) = lml(indexh-1)+1
  lml(indexh) = lmf(indexh)
  j = lml(indexh)
  i = lmf(isig)
  lmx(j) = lmx(i)+1
  lmy(j) = lmy(i)
  lmz(j) = lmz(i)
  zlm(j) = zlm(i)*coef
  nterm = nterm+1
  k = lmf(isigm)
  do i=lmf(isig)+1,lml(isig)
    lml(indexh) = lml(indexh)+1
    j = lml(indexh)
    lmx(j) = lmx(i)+1
    lmy(j) = lmy(i)
    lmz(j) = lmz(i)
    zlm(j) = (zlm(i)-zlm(k))*coef
    k = k+1
    nterm = nterm+1
  end do
  if (mod(mang,2) == 0) then
    lml(indexh) = lml(indexh)+1
    j = lml(indexh)
    k = lml(isigm)
    lmx(j) = lmx(k)
    lmy(j) = lmy(k)+1
    lmz(j) = lmz(k)
    zlm(j) = -zlm(k)*coef
    nterm = nterm+1
  end if
end do
!ixy = 0
!iz = 0
!do lang=1,lproju
!  do mang=0,lang-1
!    nsigma = min(1,mang)
!    ndelta = max(0,1-mang)
!    anum = real((lang-mang)*(lang+mang+1),kind=wp)
!    aden = real(2*(2-ndelta),kind=wp)
!    coef = sqrt(anum/aden)
!    do isigma=nsigma,0,-1
!      i_sign = 2*isigma-1
!      ixy = ixy+1
!      flmtx(1,ixy) = real(i_sign,kind=wp)*coef
!      flmtx(2,ixy) = coef
!      if (mang /= 0) then
!        iz = iz+1
!        flmtx(3,iz) = -real(mang*isigma,kind=wp)
!      end if
!    end do
!  end do
!  iz = iz+1
!  flmtx(3,iz) = -real(lang,kind=wp)
!end do
! Column and row indices for angular momentum matrix elements.
!iadd = 1
!do i=1,2*lproju-1
!  mc(1,i) = i
!  mc(2,i) = i
!  mc(3,i) = i+1
!  mr(1,i) = i+iadd
!  mr(2,i) = i+2
!  mr(3,i) = i+2
!  iadd = 4-iadd
!end do

return

end subroutine cortab
