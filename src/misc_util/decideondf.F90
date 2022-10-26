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

subroutine DecideOnDF(DoDF)

implicit none
logical DoDF
integer iOption

call Get_iScalar('System BitSwitch',iOption)
DoDF = iand(iOption,1024) == 1024

#if defined (_DEBUGPRINT_)
write(6,*) '>>> Exit from DecideOnDF:'
write(6,*) '    System Bit Switch = ',iOption
write(6,*) '    DoDF = ',DoDF
#endif

end subroutine DecideOnDF
