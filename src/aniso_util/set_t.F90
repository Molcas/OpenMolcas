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

subroutine set_T(nT,nTempMagn,TINPUT,TempMagn,Tmin,Tmax,chit_exp,Texp,T,XTexp)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nT, nTempMagn
logical(kind=iwp), intent(in) :: TINPUT
real(kind=wp), intent(in) :: TempMagn(nTempMagn), Tmin, Tmax, chit_exp(nT), Texp(nT)
real(kind=wp), intent(out) :: T(nT+nTempMagn), XTexp(nT+nTempMagn)
integer(kind=iwp) :: i
real(kind=wp) :: dltt

! set nT, T(i) and XTexp(i) arrays:
T(:) = Zero
XTexp(:) = Zero
!----------------------------------------------------------------------!
if (TINPUT) then
  ! case 1:  T(iT) computed from input values: Texp, and nTempMagn
  !          nTempMagn = 0
  if (nTempMagn > 0) T(1:nTempMagn) = TempMagn(:)
  T(nTempMagn+1:) = Texp(:)
  XTexp(nTempMagn+1:) = chit_exp(:)
else
  ! case 2:  T(iT) computed from input values: Tmin, Tmax, nT
  !          and nTempMagn
  dltt = (tmax-tmin)/real(nT-1,kind=wp)

  T(nTempMagn+1:) = Texp(:)
  XTexp(nTempMagn+1:) = chit_exp(:)

  if (nTempMagn > 0) T(1:nTempMagn) = TempMagn(:)

  T(nTempMagn+1) = 1.0e-4_wp
  do i=2,nT
    T(nTempMagn+i) = Tmin+dltt*real(i-1,kind=wp)
  end do
end if !tinput

!----------------------------------------------------------------------!

! check for T=0 values and replace them
! with a definite nonzero value:
do i=1,nT+nTempMagn
  if (abs(T(i)) <= tiny(T)) T(i) = 1.0e-4_wp
end do

return

end subroutine set_T
