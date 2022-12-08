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

subroutine ireorg3(symp,typp,typpv1,paddv1,rc)
! this routine def. constants to be added to index from v2
! to determine proper index in v1
! N.B. typp and typpv1 must be compatible (this is testet
! in this version)
!
! symp  - irrep of p index (I)
! typp  - typ of p index in v2 (I)
! typpv1- typ of corresponding p index in v1 (I)
! paddv1- constant to be added (O) pv1 = pv2+paddv1
! rc    - return (error) code (O)

use ccsort_global, only: noa, nob, nva, nvb
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: symp, typp, typpv1
integer(kind=iwp), intent(out) :: paddv1, rc

rc = 0

select case (typp)
  case (1,2)
    if ((typpv1 == 1) .or. (typpv1 == 2) .or. (typpv1 == 5)) then
      paddv1 = 0
    else
      rc = 1
      ! RC=1 : typp=1 or 2, incompatible typpv1 (Stup)
    end if
  case (3)
    if (typpv1 == 3) then
      paddv1 = 0
    else if (typpv1 == 4) then
      paddv1 = nvb(symp)-nva(symp)
    else if (typpv1 == 5) then
      paddv1 = noa(symp)
    else
      rc = 2
      ! RC=2 : typp=3, incompatible typpv1 (Stup)
    end if
  case (4)
    if (typpv1 == 4) then
      paddv1 = 0
    else if (typpv1 == 5) then
      paddv1 = nob(symp)
    else
      rc = 3
      ! RC=3 : typp=4, incompatible typpv1 (Stup)
    end if
  case (5)
    if (typpv1 == 5) then
      paddv1 = 0
    else
      ! RC=4 : typp=5, incompatible typpv1 (Stup)
    end if
  case default
    rc = 5
    ! RC=5 : improper typp (Stup)
end select

return

end subroutine ireorg3
