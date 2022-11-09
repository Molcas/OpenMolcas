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

subroutine ClsOne(rc,Option)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Close the one-electron integral file.                            *
!                                                                      *
!     input:                                                           *
!     Option : Switch to set options                                   *
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
!     M. P. Fuelscher, University of Lund, Sweden, 1993                *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use OneDat, only: AuxOne, NaN, rcOne, sDmp, TocOne
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: option
integer(kind=iwp) :: LuOne

!----------------------------------------------------------------------*
! Start procedure                                                      *
!----------------------------------------------------------------------*
rc = rcOne%good
LuOne = AuxOne%Lu
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
if (.not. AuxOne%Opn) then
  rc = rcOne%CL01
  call SysAbendMsg('ClsOne','The ONEINT file has not been opened',' ')
end if
AuxOne%Opn = .false.
!----------------------------------------------------------------------*
! Dump the TOC upon request                                            *
!----------------------------------------------------------------------*
if (btest(Option,sDmp)) call DmpOne()
!----------------------------------------------------------------------*
! Reset error code,open flag and unit number. Close file.              *
!----------------------------------------------------------------------*
call DaClos(LuOne)
AuxOne%Lu = NaN
TocOne(:) = NaN

!----------------------------------------------------------------------*
!     Terminate procedure                                              *
!----------------------------------------------------------------------*
return

end subroutine ClsOne
