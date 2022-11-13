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

subroutine ClsMCK(rc,option)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Close the one-electron integral file.                            *
!                                                                      *
!     input:                                                           *
!     option : Switch to set options                                   *
!              (not used at present)                                   *
!                                                                      *
!     output:                                                          *
!     rc     : Return code.                                            *
!              A value of 0 (zero) is returned upon successful         *
!              completion of the request. A nonzero value indicates    *
!              an error.                                               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.O. Widmark                                 *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use MckDat, only: AuxMck, pBas, pEnd, pFID, pNext, pOp, pSym, pSymOp, pTitle, pVersN, rcMck, sDmp, TocMck
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: option
integer(kind=iwp) :: LuMCK

!----------------------------------------------------------------------*
! Check file status                                                    *
!----------------------------------------------------------------------*
if (.not. AuxMck%Opn) then
  rc = rcMck%CL01
  call SysAbendMsg('ClsMCK','The MCK file has not been opened',' ')
end if
if (btest(Option,sDmp)) then
  write(u6,'(i6,z8)') pFID,TocMck(pFID)
  write(u6,'(i6,z8)') pVersN,TocMck(pVersN)
  write(u6,'(i6,z8)') pTitle,TocMck(pTitle)
  write(u6,'(i6,z8)') pOp,TocMck(pOp)
  write(u6,'(i6,z8)') pSym,TocMck(pSym)
  write(u6,'(i6,z8)') pSymOp,TocMck(pSymOp)
  write(u6,'(i6,z8)') pBas,TocMck(pBas)
  write(u6,'(i6,z8)') pNext,TocMck(pNext)
  write(u6,'(i6,z8)') pEnd,TocMck(pEnd)
end if
!----------------------------------------------------------------------*
! Reset error code,open flag and unit number. Close file.              *
!----------------------------------------------------------------------*
LuMCK = AuxMck%Lu
call DaClos(LuMCK)
AuxMck%Lu = 0
AuxMck%Opn = .false.
rc = rcMck%good

!----------------------------------------------------------------------*
! Terminate procedure                                                  *
!----------------------------------------------------------------------*
return

end subroutine ClsMCK
