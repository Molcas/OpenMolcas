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
!              completion of the request. A nonzero value indi-        *
!              cates an error.                                         *
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

implicit integer(A-Z)
#include "MckDat.fh"

!----------------------------------------------------------------------*
! Check file status                                                    *
!----------------------------------------------------------------------*
if (AuxMCK(pOpen) /= 1) then
  rc = rcCL01
  call SysAbendMsg('ClsMCK','The MCK file has not been opened',' ')
end if
if (iand(Option,1024) /= 0) then
  write(6,'(i6,z8)') pFID,TocOne(pFID)
  write(6,'(i6,z8)') pVersN,TocOne(pVersN)
  write(6,'(i6,z8)') pTitle,TocOne(pTitle)
  write(6,'(i6,z8)') pOp,TocOne(pOp)
  write(6,'(i6,z8)') pSym,TocOne(pSym)
  write(6,'(i6,z8)') pSymOp,TocOne(pSymOp)
  write(6,'(i6,z8)') pBas,TocOne(pBas)
  write(6,'(i6,z8)') pNext,TocOne(pNext)
  write(6,'(i6,z8)') pEnd,TocOne(pEnd)
end if
!----------------------------------------------------------------------*
! Reset error code,open flag and unit number. Close file.              *
!----------------------------------------------------------------------*
LuMCK = AuxMCK(pLu)
call DaClos(LuMCK)
AuxMCK(pLu) = 0
AuxMCK(pOpen) = 0
rc = rc0000

!----------------------------------------------------------------------*
! Terminate procedure                                                  *
!----------------------------------------------------------------------*
return

end subroutine ClsMCK
