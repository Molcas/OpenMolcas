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

subroutine qbess(alpha,apwr,aterm1,aterm2,binom,bpref,bpwr,bterm1,dfac,l,lambu,lmahi,lmbhi,ltot1,nu,prd,qsum,rka,rkb,ssi)
! Compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
! using the Bessel function formula for
! lama=max(l,nu) to lmahi, lamb=max(l,nu) to lmbhi, n=lama+lamb-l-l

#include "intent.fh"

use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: l, lambu, lmahi, lmbhi, ltot1, nu
real(kind=wp), intent(in) :: alpha, binom(*), dfac(*), prd, rka, rkb
real(kind=wp), intent(_OUT_) :: apwr(*), aterm1(*), aterm2(*), bpref(*), bpwr(*), bterm1(*), ssi(*)
real(kind=wp), intent(inout) :: qsum(ltot1,lambu,*)
integer(kind=iwp) :: it, iu, lamap, lambp, lami, lmaphi, lmbphi, lmihi, lmlo, lmplo, n, num1
real(kind=wp) :: dsum, fct, fcta, fctb, fctra, fctrad, fctran, fctrb, fctrbd, fctrbn, fctrt, fctrtd, fctrtn, fctru, fctrud, fctrun

! nu=l+1-npi/2
! Bessel function formula applies to all npi=2 cases, no npi=1
! cases, and some npi=0 cases.
num1 = nu-1
lmlo = max(l,nu)
lmaphi = lmahi-num1
lmbphi = lmbhi-num1
fcta = rka/(alpha+alpha)
fct = rka*fcta
fctb = rkb/(alpha+alpha)
bpref(1) = fctb**num1
do lambp=2,lmbphi
  bpref(lambp) = fctb*bpref(lambp-1)
end do
apwr(1) = One
apwr(2) = fct
do lambp=3,lmbphi
  apwr(lambp) = fct*apwr(lambp-1)
end do
fct = rkb*fctb
bpwr(1) = One
bpwr(2) = fct
do lamap=3,lmaphi
  bpwr(lamap) = fct*bpwr(lamap-1)
end do
lmihi = lmaphi+lmbphi+(nu-2)
call ssibfn(lmihi-1,rka*fctb,ssi)
do lami=1,lmihi
  ssi(lami) = ssi(lami)/dfac(lami+lami-1)
end do
lmplo = lmlo-num1
fctra = (alpha+alpha)**(nu-2)*fcta**(lmlo-1)*prd/sqrt(Four*alpha)*((dfac(2*(2*lmplo+num1)-1)/dfac(2*(lmplo+num1)+1))* &
        dfac(2*num1+1))/dfac(2*(lmplo+num1)+1)
fctran = (2*(2*lmplo+num1)-1)
fctrad = (2*(lmplo+num1)+1)
do lamap=lmplo,lmaphi
  fctru = One
  fctrun = (nu+num1)
  fctrud = (2*(lamap+num1)+1)
  do iu=1,lmbphi
    aterm1(iu) = fctru*apwr(iu)
    fctru = fctru*fctrun/fctrud
    fctrun = fctrun+Two
    fctrud = fctrud+Two
  end do
  do it=1,lamap
    bterm1(it) = binom((lamap*(lamap-1))/2+it)*bpwr(it)
  end do
  fctrb = fctra
  fctrbn = (2*(lamap+lmplo+num1)-1)
  fctrbd = (2*(lmplo+num1)+1)
  do lambp=lmplo,lmbphi
    n = ((2*(nu-l)-1)+lamap)+lambp
    do iu=1,lambp
      aterm2(iu) = binom((lambp*(lambp-1))/2+iu)*aterm1(iu)
    end do
    dsum = Zero
    fctrt = One
    fctrtn = (nu+num1)
    fctrtd = (2*(lambp+num1)+1)
    do it=1,lamap
      do iu=1,lambp
        dsum = dsum+aterm2(iu)*(fctrt*bterm1(it))*ssi((it+(nu-2))+iu)
      end do
      fctrt = fctrt*fctrtn/fctrtd
      fctrtn = fctrtn+Two
      fctrtd = fctrtd+Two
    end do
    qsum(n,lambp+num1,lamap+num1) = qsum(n,lambp+num1,lamap+num1)+fctrb*bpref(lambp)*dsum
    fctrb = fctrb*fctrbn/fctrbd
    fctrbn = fctrbn+Two
    fctrbd = fctrbd+Two
  end do
  fctra = fcta*fctra*fctran/fctrad
  fctran = fctran+Two
  fctrad = fctrad+Two
end do

return

end subroutine qbess
