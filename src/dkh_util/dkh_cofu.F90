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

subroutine dkh_cofu(DKH_order,DKH_para,dkcof)
! Calculate coefficients of different parameterization of DKH unitary transformation

use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: DKH_order, DKH_para
real(kind=wp), intent(out) :: dkcof(max(4,DKH_order))
integer(kind=iwp) :: i, k, n
real(kind=wp) :: c, d

N = max(4,DKH_order)
if (DKH_para == 2) then
  ! Exponential
  dkcof(1) = One
  do i=2,N
    dkcof(i) = dkcof(i-1)/real(i,kind=wp)
  end do
else if (DKH_para == 3) then
  ! Square root
  dkcof(:) = Zero
  dkcof(1) = One
  dkcof(2) = Half
  do i=4,N,2
    dkcof(i) = -dkcof(i-2)*(i-3)/real(i,kind=wp)
  end do
else if (DKH_para == 5) then
  ! Cayley
  dkcof(1) = One
  do i=2,N
    dkcof(i) = dkcof(i-1)*Half
  end do
else if (DKH_para == 4) then
  ! McWeeny
  dkcof(1) = One
  dkcof(2) = Half
  dkcof(3) = Half
  do i=4,N,2
    dkcof(i) = dkcof(i-2)*(i-1)/real(i,kind=wp)
    if (i < N) then
      dkcof(i+1) = dkcof(i)
    end if
  end do
else if (DKH_para == 1) then
  ! opt ( M. Reiher )
  dkcof(1) = One
  dkcof(2) = Half
  dkcof(3) = (Two-sqrt(Two))*Quart
  dkcof(4) = dkcof(3)-0.125_wp
  do i=5,N,2
    ! W(i+3) in terms of a(k),k<i
    d = Zero
    do k=(i+3)/2,i-1
      c = dkcof(k)*dkcof(i+3-k)
      if (k /= (i+3)/2) c = c*Two
      if (mod(k,2) == 1) c = -c
      d = d-c
    end do
    dkcof(i) = d*sqrt(Two)
    if (i < n) then
      dkcof(i+1) = dkcof(i)
    end if
  end do
end if

return

end subroutine dkh_cofu
