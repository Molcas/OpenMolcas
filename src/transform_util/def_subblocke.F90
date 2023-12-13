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
! Copyright (C) Giovanni Ghigo                                         *
!***********************************************************************

subroutine Def_SubBlockE(iSymA,iSymB)
!***********************************************************************
! Author  :  Giovanni Ghigo                                            *
!            Lund University, Sweden                                   *
!----------------------------------------------------------------------*
! Define the SubBlocks to calculate in the Exchange matrix.            *
!***********************************************************************

use Cho_Tra, only: DoTCVA, nAsh, nIsh, nSsh, SubBlocks
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iSymA, iSymB

SubBlocks(:,:) = .false.
if (DoTCVA) then
  if (nIsh(iSymA) > 0) then
    if (nIsh(iSymB) > 0) SubBlocks(1,1) = .true.
    if (nAsh(iSymB) > 0) SubBlocks(1,2) = .true.
    if (nSsh(iSymB) > 0) SubBlocks(1,3) = .true.
  end if
  if (nAsh(iSymA) > 0) then
    if (nIsh(iSymB) > 0) SubBlocks(2,1) = .true.
    if (nAsh(iSymB) > 0) SubBlocks(2,2) = .true.
    if (nSsh(iSymB) > 0) SubBlocks(2,3) = .true.
  end if
  if (nSsh(iSymA) > 0) then
    if (nIsh(iSymB) > 0) SubBlocks(3,1) = .true.
    if (nAsh(iSymB) > 0) SubBlocks(3,2) = .true.
  end if
end if
if ((nSsh(iSymA) > 0) .and. (nSsh(iSymB) > 0)) SubBlocks(3,3) = .true.

return

end subroutine Def_SubBlockE
