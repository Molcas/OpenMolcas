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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!               2002, Roland Lindh                                     *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Copy all input to a local scratch file                           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!     Modified by:  R. Lindh to allow for use of old input file.       *
!                   Tokyo, Japan, June 2002.                           *
!                   J.W. Krogh to make auto.plx create the old input   *
!                   file. Lund, Sweden, October 2005.                  *
!***********************************************************************
Module Spool
Private

Logical, public :: Spool_On=.True.
Integer, public :: LuRd=5, LuWr=6

Public:: SpoolInp, Set_Spool, Disable_Spool, Close_LuSpool

contains
subroutine SpoolInp(LuSpool)

use UnixInfo, only: ProgName
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: LuSpool
integer(kind=iwp) :: iEnd
logical(kind=iwp) :: Exists
character(len=len(ProgName)) :: PName
character(len=128) :: FileName
integer(kind=iwp) :: IsFreeUnit

! Get the name of the module

PName = ProgName
call Upcase(PName)
PName = adjustl(PName)

iEnd = 1
do while (PName(iEnd:iEnd) /= ' ')
  iEnd = iEnd+1
end do
iEnd = min(iEnd-1,5)
FileName = PName(1:iend)//'INP'

! If Spool_on is true, then the input in StdIn is just used.
! Else we'll use the latest input. This is created by auto.plx

LuSpool = 17
if (Spool_On) then
  LuSpool = LuRd
else
  call f_inquire('LASTEN',Exists) ! customized Last_Energy input
  if (Exists) then
    LuSpool = IsFreeUnit(LuSpool)
    call Molcas_Open(LuSpool,'LASTEN')
  else
    call f_inquire(Filename,Exists)
    if (Exists) then
      LuSpool = IsFreeUnit(LuSpool)
      call Molcas_Open(LuSpool,Filename)
    end if
  end if
end if

end subroutine SpoolInp

subroutine Set_Spool()

implicit none

Spool_on = .True.

end subroutine Set_Spool

subroutine Disable_Spool()

implicit none

Spool_on = .false.

end subroutine Disable_Spool

subroutine Close_LuSpool(LuSpool)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LuSpool

if (.not. Spool_On) close(LuSpool)

end subroutine Close_LuSpool

End Module Spool
