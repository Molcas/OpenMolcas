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

subroutine ireorg2(symp,typp,pup,rc)
! this routine def. summation limits for given symp and typp
! i.e. number of indices for this symmetry and typ
!
! symp  - irrep of p index (I)
! typp  - typ of p index in v2 (I)
! pup   - summation limit
! rc    - return (error) code (O)

use ccsort_global, only: noa, nob, NORB, nva, nvb
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: symp, typp
integer(kind=iwp), intent(out) :: pup, rc

rc = 0

select case (typp)
  case (1)
    pup = noa(symp)
  case (2)
    pup = nob(symp)
  case (3)
    pup = nva(symp)
  case (4)
    pup = nvb(symp)
  case (5)
    pup = norb(symp)
  case default
    rc = 1
    ! RC=1 : bad typp (Stup)
end select

return

end subroutine ireorg2
