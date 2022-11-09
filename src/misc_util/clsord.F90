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

subroutine ClsOrd(rc)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Close the two-electron integral file.                            *
!                                                                      *
!     input:                                                           *
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
!     M. P. Fuelscher                                                  *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use TwoDat, only: AuxTwo, iNoNum, lTocTwo, RAMD, rcTwo, TocTwo
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: iDisk, LuOrd

!----------------------------------------------------------------------*
! Start procedure:                                                     *
!----------------------------------------------------------------------*
rc = rcTwo%good
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
if (.not. AuxTwo%Opn) then
  rc = rcTwo%CL01
  call SysAbendMsg('ClsOrd','The ORDINT file has not been opened',' ')
end if
!----------------------------------------------------------------------*
! Reset error code,open flag and unit number. Close file.              *
!----------------------------------------------------------------------*
LuOrd = AuxTwo%Unt
iDisk = 0
call iDaFile(LuOrd,1,TocTwo,lTocTwo,iDisk)
call DaClos(LuOrd)
AuxTwo%Opn = .false.
AuxTwo%Unt = iNoNum
AuxTwo%DaDa = iNoNum

if (RAMD%act) RAMD%act = .false.

!----------------------------------------------------------------------*
! Terminate procedure                                                  *
!----------------------------------------------------------------------*
return

end subroutine ClsOrd
