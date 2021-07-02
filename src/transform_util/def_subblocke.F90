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

use Cho_Tra

implicit real*8(a-h,o-z)
implicit integer(i-n)

do i=1,3
  do j=1,3
    SubBlocks(i,j) = .false.
  end do
end do
if (DoTCVA .and. (nIsh(iSymA) > 0)) then
  if (nIsh(iSymB) > 0) SubBlocks(1,1) = .true.
  if (nAsh(iSymB) > 0) SubBlocks(1,2) = .true.
  if (nSsh(iSymB) > 0) SubBlocks(1,3) = .true.
end if
if (DoTCVA .and. (nAsh(iSymA) > 0)) then
  if (nIsh(iSymB) > 0) SubBlocks(2,1) = .true.
  if (nAsh(iSymB) > 0) SubBlocks(2,2) = .true.
  if (nSsh(iSymB) > 0) SubBlocks(2,3) = .true.
end if
if (DoTCVA .and. (nSsh(iSymA) > 0)) then
  if (nIsh(iSymB) > 0) SubBlocks(3,1) = .true.
  if (nAsh(iSymB) > 0) SubBlocks(3,2) = .true.
end if
if ((nSsh(iSymA)*nSsh(iSymB)) > 0) SubBlocks(3,3) = .true.

return

end subroutine Def_SubBlockE
