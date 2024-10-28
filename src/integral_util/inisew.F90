!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1999, Roland Lindh                                     *
!***********************************************************************

subroutine IniSew(DSCF,nDiff)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Chemical Physics,                 *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Basis_Info, only: Seward_Activated
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: DSCF
integer(kind=iwp), intent(inout) :: nDiff

if (Seward_Activated) then
  call ClsSew()
  call xRlsMem_Ints()
end if

call Seward_Init()

call GetInf(DSCF,nDiff)

return

end subroutine IniSew
