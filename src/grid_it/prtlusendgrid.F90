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

subroutine PRTLUSENDGRID(LUVAL)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

character LINE*128

write(LINE,'(A)') ' </INPORB>'
call PRINTLINE(LUVAL,LINE,10,0)
write(LINE,'(A)') ' </GRID>'
call PRINTLINE(LUVAL,LINE,8,0)

return

end subroutine PRTLUSENDGRID
