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

subroutine Size_SO_block_g(iSD4,nSD,nSO,No_batch)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4)
integer(kind=iwp), intent(out) :: nSO
logical(kind=iwp), intent(out) :: No_batch
integer(kind=iwp), external :: MemSO2_P

nSO = MemSO2_P(iSD4(2,1),iSD4(2,2),iSD4(2,3),iSD4(2,4),iSD4(7,1),iSD4(7,2),iSD4(7,3),iSD4(7,4))
No_batch = nSO == 0

return

end subroutine Size_SO_block_g
