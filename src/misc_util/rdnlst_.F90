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
!***********************************************************************

subroutine RdNLst_(iUnit,NameIn,No_Input_OK)
!***********************************************************************
!                                                                      *
!     Locate the beginning of an input stream                          *
!     (similar to FORTRAN NAMELIST read known to some systems)         *
!                                                                      *
!     calling arguments:                                               *
!     iUnit  : Type integer, input                                     *
!              FORTRAN unit number                                     *
!     NameIn : Type character string, input                            *
!              Character string marking the beginning of the input     *
!     No_Input_OK: Logical                                             *
!                  On input determines if an input has to be found     *
!                  On exit determines if an input was found.           *
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

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iUnit
character(len=*), intent(in) :: NameIn
logical(kind=iwp), intent(inout) :: No_Input_OK
integer(kind=iwp) :: istatus, lStdNam
character(len=80) :: Line
character(len=8) :: StdNam

call ResetErrorLine()
!----------------------------------------------------------------------*
! push the entry name on the calling stack                             *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! convert the Name to internal standard format.                        *
!----------------------------------------------------------------------*
call StdFmt(NameIn,StdNam)
lStdNam = len_trim(StdNam)
!----------------------------------------------------------------------*
! read until an input Line is located which starts with                *
! the string, Name, not before the second column                       *
!----------------------------------------------------------------------*
do
  read(iUnit,'(A)',iostat=istatus) Line
  call Error(istatus)
  call UpCase(Line)
  Line = adjustl(Line)
  if ((Line(1:1) == '&') .and. (Line(2:lStdNam+1) == StdNam(1:lStdNam))) exit
end do

contains

!----------------------------------------------------------------------*
! error exit                                                           *
!----------------------------------------------------------------------*
subroutine Error(rc)

  integer(kind=iwp), intent(in) :: rc

  select case (rc)
    case (0)
      ! do nothing
    case (:-1)
      if (No_Input_OK) then
        No_Input_OK = .false.
      else
        write(u6,*) 'RdNLst: Input section not found in input file'
        write(u6,*) '        Looking for:',StdNam(1:lStdNam)
        call Quit_OnUserError()
      end if
    case default
      call Abend()
  end select

end subroutine Error

end subroutine RdNLst_
