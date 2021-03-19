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

integer function LineNr(Lu)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

implicit none
integer Lu, ios
real*8 foo

LineNr = 0
do while (.true.)
  read(Lu,*,iostat=ios) foo,foo,foo,foo
  if (ios > 0) then
    call Abend() !'problem somewhere'
  else if (ios < 0) then ! end of file is reached
    return
  else
    LineNr = LineNr+1
  end if
end do

return
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_real(foo)
#endif

end function LineNr
