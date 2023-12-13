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

subroutine rm_AuxShell(iCnttp)
!***********************************************************************
!                                                                      *
!     Remove an auxiliary basis set by making it empty.                *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, Shells
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iCnttp
integer(kind=iwp) :: iShll, k

!                                                                      *
!***********************************************************************
!                                                                      *
do k=0,dbsc(iCnttp)%nShells-1
  iShll = dbsc(iCnttp)%iVal+k

  Shells(iShll)%nExp = 0
  Shells(iShll)%nBasis = 0
  Shells(iShll)%nBasis_c = 0

end do
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine rm_AuxShell
