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

subroutine Cho_P_GetGSP(nGSP)
!
! Purpose: get global number of contributing shell pairs.

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nGSP
#include "cho_para_info.fh"

if (Cho_Real_Par) then
  call Cho_P_GetGSP_P(nGSP)
else
  call Cho_P_GetGSP_S(nGSP)
end if

end subroutine Cho_P_GetGSP
