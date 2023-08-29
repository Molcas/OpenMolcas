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

subroutine Cho_Inp_SetDecAlg(ForceParallel)

use Cholesky, only: Cho_DecAlg, Cho_Real_Par
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: ForceParallel

if (Cho_Real_Par .or. ForceParallel) then
  if (Cho_DecAlg == 1) then
    Cho_DecAlg = 4  ! parallel one-step
  else if (Cho_DecAlg == 2) then
    Cho_DecAlg = 5 ! parallel two-step
  else if (Cho_DecAlg == 3) then
    Cho_DecAlg = 6  ! parallel naive
  end if
end if

end subroutine Cho_Inp_SetDecAlg
