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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************
!  handle2name
!
!> @brief
!>   Retrieve the file name from molcas I/O
!> @author Valera Veryazov
!>
!> @details
!> The routine can be called from ::aixrd or ::aixwr.
!> Return the file name, or '``Unknown``'.
!>
!> @param[in]  handle   file handle
!> @param[out] filename file name
!***********************************************************************

subroutine handle2name(handle,filename)

use Fast_IO, only: CtlBlk, FCtlBlk, MxFile, pHndle
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: handle
character(len=*), intent(out) :: filename
integer(kind=iwp) :: i

filename = 'Unknown'
do i=1,MxFile
  if (CtlBlk(pHndle,i) == handle) then
    filename = FCtlBlk(i)
    exit
  end if
end do

return

end subroutine handle2name
