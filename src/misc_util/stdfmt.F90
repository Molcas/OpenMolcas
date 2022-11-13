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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               1993, Per-Olof Widmark                                 *
!***********************************************************************

subroutine StdFmt(NameIn,NameUt)
!***********************************************************************
!                                                                      *
!     Convert the first word of the incoming string, NameIn, to        *
!     standard format (no leading blanks, all upper case letters).     *
!                                                                      *
!     calling arguments                                                *
!     NameIn  : Type character string, input                           *
!     NameUt  : Type character string, output                          *
!               First token of NameIn in upper case letters            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and P.O. Widmark                                  *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: NameIn
character(len=*), intent(out) :: NameUt
integer(kind=iwp) :: i, j, lIn, lLeft, lUt

!----------------------------------------------------------------------*
! start                                                                *
!----------------------------------------------------------------------*
lIn = len(NameIn)
lUt = len(NameUt)
!----------------------------------------------------------------------*
! set NameUt to a blank string                                         *
!----------------------------------------------------------------------*
NameUt = ' '
!do i=1,lUt
!  NameUt(i:i) = ' '
!end do
!----------------------------------------------------------------------*
! get the number of leading blanks                                     *
!----------------------------------------------------------------------*
!lLeft = 0
do i=1,lIn
  if (NameIn(i:i) /= ' ') exit
end do
lLeft = i

!do i=lIn,1,-1
!  if (NameIn(i:i) /= ' ') lLeft = i-1
!end do
!----------------------------------------------------------------------*
! copy only the first token of NameIn to NameUt                        *
!----------------------------------------------------------------------*
j = 0
do i=lLeft,lIn
  if ((NameIn(i:i) == ' ') .or. (j == lUt)) exit
  j = j+1
  NameUt(j:j) = NameIn(i:i)
end do
!----------------------------------------------------------------------*
! change lower case character to upper case                            *
!----------------------------------------------------------------------*
call UpCase(NameUt)

!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine StdFmt
