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
      SUBROUTINE PRTLUSENDGRID(LUVAL)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
      CHARACTER LINE*128
      WRITE(LINE,'(A)') ' </INPORB>'
      CALL PRINTLINE(LUVAL, LINE,10,0)
      WRITE(LINE,'(A)') ' </GRID>'
      CALL PRINTLINE(LUVAL, LINE,8,0)
      END
