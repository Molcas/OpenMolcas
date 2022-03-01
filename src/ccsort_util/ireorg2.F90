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

#include "ccsort.fh"
integer symp, typp, pup, rc

if (typp == 1) then
  pup = noa(symp)
else if (typp == 2) then
  pup = nob(symp)
else if (typp == 3) then
  pup = nva(symp)
else if (typp == 4) then
  pup = nvb(symp)
else if (typp == 5) then
  pup = norb(symp)
else
  rc = 1
  ! RC=1 : bad typp (Stup)
  return
end if

return

end subroutine ireorg2
