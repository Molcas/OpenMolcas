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

function R_Stab_A(R,S,nS)

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: R_Stab_A
integer(kind=iwp), intent(in) :: nS, R, S(nS)
integer(kind=iwp) :: iS

R_Stab_A = .false.
do iS=1,nS
  if (R == S(iS)) then
    R_Stab_A = .true.
    return
  end if
end do

return

end function R_Stab_A
