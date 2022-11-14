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

subroutine DecideOnLocalDF(LocalDF)

use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(out) :: LocalDF
integer(kind=iwp) :: DFMode
logical(kind=iwp) :: DF

call DecideOnDF(DF)
if (DF) then
  call Get_iScalar('DF Mode',DFMode)
  LocalDF = DFMode == 1
else
  LocalDF = .false.
end if

end subroutine DecideOnLocalDF
