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

subroutine Cho_P_GetGSP_P(nGSP)

use Cholesky, only: nnShl_G
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nGSP

nGSP = nnShl_G

end subroutine Cho_P_GetGSP_P
