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
! Copyright (C) 1991, Per-Olof Widmark                                 *
!               1993, Markus P. Fuelscher                              *
!               1996, Luis Serrano-Andres                              *
!               2012, Victor P. Vysotskiy                              *
!***********************************************************************

subroutine FIOInit()
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     set initial values of all arguments in the common blocks         *
!     /FIO1/, /FIO2/ and /FIO3/                                        *
!     /PFIO1/,/PFIO2/,/PFIO3/,/PFIO4/                                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, IBM Sweden, 1991                                   *
!     M.P. Fuelscher, University of Lund, Sweden, 1993                 *
!     L. Serrano-Andres, University of Lund, Sweden, 1996              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
!     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
!                                                                      *
!***********************************************************************

use Fast_IO, only: Addr, FSCB, isOpen, LuName, MaxFileSize, MBL, MPUnit, Multi_File, MxFile, ProfData, Query, Trace
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i

isOpen(:) = 0
FSCB(:) = 0
Addr(:) = 0
Multi_File(:) = .false.
ProfData(:,:) = Zero
MPUnit(:,:) = 0
MBL(:) = 0
LuName(:) = 'FT__F001'
do i=1,MxFile
  write(LuName(i)(3:4),'(I2.2)') i
end do

MaxFileSize = 0
Trace = .false.
Query = .false.

return

end subroutine FIOInit
