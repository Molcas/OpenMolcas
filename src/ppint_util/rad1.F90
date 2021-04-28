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

subroutine rad1(aa,aarr1,alpt,arp2,ccr,dfac,fctr2,kcrl,kcru,lamu,ltot1,ncr,qsum,rk,tol,zcr)
! compute type 1 radial integrals

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: kcrl, kcru, lamu, ltot1, ncr(*)
real(kind=wp), intent(in) :: aa, aarr1, alpt, arp2, ccr(*), dfac(*), fctr2, rk, tol, zcr(*)
real(kind=wp), intent(inout) :: qsum(ltot1,lamu)
integer(kind=iwp) :: kcr, lam, n, npi
real(kind=wp) :: alpha, dum, f2lam3, prd, t
real(kind=wp), external :: qcomp

do kcr=kcrl,kcru
  npi = ncr(kcr)
  alpha = aa+zcr(kcr)
  ! exponential factor from q functions included in dum
  dum = aarr1+zcr(kcr)*arp2/alpha
  if (dum > tol) cycle
  prd = fctr2*ccr(kcr)*exp(-dum)
  if (rk == Zero) then
    t = Zero
    do n=1,ltot1-mod(ltot1-1,2),2
      qsum(n,1) = qsum(n,1)+prd*qcomp(alpha,dfac,npi+n-1,0,t,rk)
    end do
  else
    t = alpt/alpha
    do lam=1,lamu
      qsum(lam,lam) = qsum(lam,lam)+prd*qcomp(alpha,dfac,npi+lam-1,lam-1,t,rk)
    end do
  end if
end do

if (rk /= Zero) then
  f2lam3 = (lamu+lamu-3)
  do lam=lamu-2,1,-1
    do n=lam+2,lamu-mod(lamu-lam,2),2
      qsum(n,lam) = qsum(n,lam+2)+(f2lam3/rk)*qsum(n-1,lam+1)
    end do
    f2lam3 = f2lam3-Two
  end do
end if

return

end subroutine rad1
