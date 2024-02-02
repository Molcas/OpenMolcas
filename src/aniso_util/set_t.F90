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

implicit none
integer, parameter :: wp = kind(0.d0)
! input:
integer, intent(in) :: nT, nTempMagn
logical, intent(in) :: TINPUT
real(kind=8), intent(in) :: Tmin, Tmax, TempMagn(nTempMagn), Texp(nT), chit_exp(nT)
real(kind=8), intent(out) :: T(nT+nTempMagn), XTexp(nT+nTempMagn)
! local variables:
integer :: i
real(kind=8) :: dltt

! set nT, T(i) and XTexp(i) arrays:
T = 0.0_wp
XTexp = 0.0_wp
!----------------------------------------------------------------------!
if (TINPUT) then
  ! case 1:  T(iT) computed from input values: Texp, and nTempMagn
  !          nTempMagn = 0
  if (nTempMagn > 0) then
    do i=1,nTempMagn
      T(i) = TempMagn(i)
    end do
    do i=1,nT
      T(i+nTempMagn) = Texp(i)
      XTexp(i+nTempMagn) = chit_exp(i)
    end do
  else
    do i=1,nT
      T(i) = Texp(i)
      XTexp(i) = chit_exp(i)
    end do
  end if
else
  ! case 2:  T(iT) computed from input values: Tmin, Tmax, nT
  !          and nTempMagn
  dltt = 0.0_wp
  dltt = (tmax-tmin)/(dble(nT-1))

  if (nTempMagn > 0) then
    do i=1,nTempMagn
      T(i) = TempMagn(i)
    end do

    T(1+nTempMagn) = 0.0001_wp
    do i=2,nT
      T(i+nTempMagn) = Tmin+dltt*dble(i-1)
    end do
  else !compute_magnetization

    T(1) = 0.0001_wp
    do i=2,nT
      T(i) = Tmin+dltt*dble(i-1)
    end do
  end if
end if !tinput

!----------------------------------------------------------------------!

! check for T=0 values and replace them
! with a definite nonzero value:
do i=1,nT+nTempMagn
  if (abs(T(i)) <= tiny(0.0_wp)) T(i) = 0.0001_wp
end do

return

end subroutine set_T
