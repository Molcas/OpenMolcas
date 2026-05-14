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

subroutine Set_Basis_Mode_Atomic(i,j)

use Basis_Info, only: dbsc
use BasisMode, only: Atomic, Auxiliary_Mode, Basis_Mode, kCnttp, lCnttp, Valence_Mode
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: i, j
integer(kind=iwp) :: k

if (dbsc(i)%Aux) then
  Basis_Mode = Auxiliary_Mode
else
  Basis_Mode = Valence_Mode
end if

do k=i+1,j
  if (dbsc(i)%Aux .neqv. dbsc(k)%Aux) then
    call WarningMessage(2,'dbsc(i)%Aux /= dbsc(k)%Aux')
    call Abend()
  end if
end do

Atomic = .true.
kCnttp = i
lCnttp = j

end subroutine Set_Basis_Mode_Atomic
