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

subroutine ptwt(arc2,dfac,npi,l,lambu,ltot1,lmahi,lmbhi,alpha,rc,rka,rkb,prd,hpt,hwt,qsum)
! Compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
! using the points and weights method,
! for lama=l to lmahi, lamb=l to lmbhi, n=lama+lamb-l-l

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: npi, l, lambu, ltot1, lmahi, lmbhi
real(kind=wp), intent(in) :: arc2, dfac(*), alpha, rc, rka, rkb, hpt(*), hwt(*)
real(kind=wp), intent(inout) :: prd, qsum(ltot1,lambu,lmahi)
integer(kind=iwp) :: i, idif, lama, lamb, n, npt
real(kind=wp) :: fctr, pt, sqalp
real(kind=wp), allocatable :: abess(:), bbess(:), ptpow(:), q2(:,:)

call mma_allocate(abess,lmahi,label='abess')
call mma_allocate(bbess,lmbhi,label='bbess')
call mma_allocate(ptpow,ltot1,label='ptpow')
call mma_allocate(q2,lambu,lmahi,label='q2')

q2(:,:) = Zero
if (arc2 > 50000.0_wp) then
  npt = 5
  idif = 0
else if (arc2 > 500.0_wp) then
  npt = 10
  idif = 5
else
  npt = 20
  idif = 15
end if
sqalp = sqrt(alpha)
prd = prd/sqalp
do i=1,npt
  pt = rc+hpt(i+idif)/sqalp
  call ssibfn(lmahi-1,rka*pt,abess)
  call ssibfn(lmbhi-1,rkb*pt,bbess)
  if ((npi+l+l-2) == 0) then
    ptpow(1) = prd
  else
    ptpow(1) = prd*pt**(npi+l+l-2)
  end if
  do n=2,ltot1
    ptpow(n) = (pt*pt)*ptpow(n-1)
  end do
  do lama=l,lmahi
    do lamb=l,lmbhi
      n = ((1-l-l)+lama)+lamb
      q2(lamb,lama) = q2(lamb,lama)+(hwt(i+idif)*abess(lama))*bbess(lamb)*ptpow(n)
    end do
  end do
end do
fctr = rkb**(l-1)
do lamb=l,lmbhi
  bbess(lamb) = fctr/dfac(lamb+lamb+1)
  fctr = rkb*fctr
end do
fctr = rka**(l-1)
do lama=l,lmahi
  do lamb=l,lmbhi
    n = ((1-l-l)+lama)+lamb
    qsum(n,lamb,lama) = qsum(n,lamb,lama)+(fctr/dfac(lama+lama+1))*bbess(lamb)*q2(lamb,lama)
  end do
  fctr = rka*fctr
end do

call mma_deallocate(abess)
call mma_deallocate(bbess)
call mma_deallocate(ptpow)
call mma_deallocate(q2)

return

end subroutine ptwt
