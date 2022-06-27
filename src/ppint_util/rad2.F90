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

subroutine rad2(ccr,kcrl,kcru,l,lambu,lmahi,lmalo,lmbhi,lmblo,ltot1,ncr,qsum,rka,rkb,zcr,lit,ljt,ca,cb,tai,taj,aa,aarr2,fctr2)
! compute type 2 radial integrals.
!
! 28-nov-90 new version replaced old version. -rmp
! 19-jan-97 put Bessel formula into a separate subroutine, qbess. -rmp

use ppint_arrays, only: binom, dfac, hpt, hwt
use Constants, only: Zero, One, Two, Four, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: kcrl, kcru, l, lambu, lmahi, lmalo, lmbhi, lmblo, ltot1, ncr(*), lit, ljt
real(kind=wp), intent(in) :: ccr(*), rka, rkb, zcr(*), ca, cb, tai, taj, aa, aarr2, fctr2
real(kind=wp), intent(out) :: qsum(ltot1,lambu,lmahi)
integer(kind=iwp) :: kcr, lama, lamb, ldifa1, ldifb, n, nhi, nlim, nlo, npi, nu
real(kind=wp) :: alpha, arc2, dum, f2lma3, f2lmb3, prd, qlim, rc, rk, t
real(kind=wp), parameter :: eps1 = 1.0e-15_wp, tol = 20.0_wp*log(Ten)
real(kind=wp), external :: qcomp

qsum(:,:,:) = Zero
! sum over potential terms
do kcr=kcrl,kcru
  npi = ncr(kcr)
  alpha = aa+zcr(kcr)
  rc = (rka+rkb)/(alpha+alpha)
  arc2 = alpha*rc*rc
  dum = aarr2+zcr(kcr)*arc2/aa
  if (dum > tol) cycle
  prd = fctr2*ccr(kcr)*exp(-dum)

  if ((rka == Zero) .and. (rkb == Zero)) then
    ! rka=0 and rkb=0
    rk = Zero
    t = Zero
    qsum(ltot1,1,1) = qsum(ltot1,1,1)+prd*qcomp(alpha,dfac,npi+ltot1-1,0,t,rk)
  else if (rka == Zero) then
    ! rka=0 and rkb>0
    rk = rkb
    t = arc2
    do lamb=l,lmbhi
      qsum(lamb-l+lit,lamb,1) = qsum(lamb-l+lit,lamb,1)+prd*qcomp(alpha,dfac,npi+lamb-l+lit-1,lamb-1,t,rk)
    end do
  else if (rkb == Zero) then
    ! rka>0 and rkb=0
    rk = rka
    t = arc2
    do lama=l,lmahi
      qsum(lama-l+ljt,1,lama) = qsum(lama-l+ljt,1,lama)+prd*qcomp(alpha,dfac,npi+lama-l+ljt-1,lama-1,t,rk)
    end do
  else if (npi == 2) then
    ! rka>0 and rkb>0; use bessel function formula.
    ! To be applicable for a set of integrals, must have nu<=l and
    ! nu==(integer), where nu=l+1-npi/2, so it is used here only
    ! for the npi=2 case. It can't be used at all for npi=(odd) and
    ! only for partial sets for npi=0
    nu = l
    call qbess(alpha,binom,dfac,l,lambu,lmahi,lmbhi,ltot1,nu,prd,qsum,rka,rkb)
  else if (arc2 >= 50.0_wp) then
    ! rka>0 and rkb>0; use pts and wts method
    ! estimate radial integrals and compare to threshold
    qlim = abs(prd)/(max(One,(rc+rc)*rka)*max(One,(rc+rc)*rkb))*sqrt(Four*(tai+tai)**lit*(taj+taj)**ljt*sqrt(tai*taj)/alpha)
    if (rc < ca) then
      nlim = npi
      qlim = qlim*ca**(lit-1)
    else
      nlim = npi+(lit-1)
    end if
    if (rc < cb) then
      qlim = qlim*cb**(ljt-1)
    else
      nlim = nlim+(ljt-1)
    end if
    if (qlim*rc**nlim >= eps1) then
      call ptwt(arc2,dfac,npi,l,lambu,ltot1,lmahi,lmbhi,alpha,rc,rka,rkb,prd,hpt,hwt,qsum)
    end if
  else
    ! rka>0 and rkb>0; use partially asymptotic method
    call qpasy(alpha,dfac,npi,l,lambu,lmahi,lmbhi,ltot1,rka,rkb,fctr2*ccr(kcr),dum+arc2,qsum)
  end if
end do

if ((rka == Zero) .and. (rkb /= Zero)) then
  ! rka=0 and rkb>0
  f2lmb3 = real(2*lmbhi-3,kind=wp)
  do lamb=lmbhi-2,lmblo,-1
    nlo = abs(lamb-l+1)+lit+1
    nhi = ljt-mod((ljt-1)-abs(lamb-l),2)+lit-1
    do n=nlo,nhi,2
      qsum(n,lamb,1) = qsum(n,lamb+2,1)+(f2lmb3/rkb)*qsum(n-1,lamb+1,1)
    end do
    f2lmb3 = f2lmb3-Two
  end do
else if ((rka /= Zero) .and. (rkb == Zero)) then
  ! rka>0 and rkb=0
  f2lma3 = real(2*lmahi-3,kind=wp)
  do lama=lmahi-2,lmalo,-1
    nlo = abs(lama-l+1)+ljt+1
    nhi = lit-mod((lit-1)-abs(lama-l),2)+ljt-1
    do n=nlo,nhi,2
      qsum(n,1,lama) = qsum(n,1,lama+2)+(f2lma3/rka)*qsum(n-1,1,lama+1)
    end do
    f2lma3 = f2lma3-Two
  end do
else if ((rka /= Zero) .and. (rkb /= Zero)) then
  ! rka>0 and rkb>0
  f2lma3 = real(lmahi+lmahi+1,kind=wp)
  do lama=lmahi,lmalo,-1
    ldifa1 = abs(l-lama)+1
    f2lmb3 = real(2*lmbhi+1,kind=wp)
    do lamb=lmbhi,lmblo,-1
      ldifb = abs(l-lamb)
      nlo = ldifa1+ldifb
      nhi = (ltot1-mod(lit-ldifa1,2))-mod((ljt-1)-ldifb,2)
      do n=nlo,nhi,2
        if (n-(lama+lamb) == (1-l-l)) cycle
        if ((lama > lmahi-2) .or. (n <= abs(l-lama-2)+ldifb)) then
          ! lamb recursion
          qsum(n,lamb,lama) = qsum(n,lamb+2,lama)+(f2lmb3/rkb)*qsum(n-1,lamb+1,lama)
        else
          ! lama recursion
          qsum(n,lamb,lama) = qsum(n,lamb,lama+2)+(f2lma3/rka)*qsum(n-1,lamb,lama+1)
        end if
      end do
      f2lmb3 = f2lmb3-Two
    end do
    f2lma3 = f2lma3-Two
  end do
end if

return

end subroutine rad2
