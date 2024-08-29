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
!---------------------------------------------------------------------*
!                                                                     *
!     lDaRec       : minimal record (in 4Byte words) transfered by    *
!                    DAFILE                                           *
!     nSect        : integer multiple of lDaRec which constitutes the *
!                    standard record length                           *
!     lStRec       : standard record length in units of integers      *
!     DiskMx_Byte  : Max disk space -1 in units of Byte               *
!                                                                     *
!---------------------------------------------------------------------*

module IOBUF

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: lDaRec = 2**10, Mode_None = 111111111, Mode_Read = 987654321, Mode_Write = 198765432, &
                                nSect = 2**5, lStRec = nSect*lDaRec
integer(kind=iwp) :: iBuf, id, iPos, iStatIO, lBuf, LuTmp, nBuf
real(kind=wp) :: Disk, Disk_1, Disk_2, DiskMx_Byte
logical(kind=iwp) :: InCore, IODone, OnDisk
real(kind=wp), allocatable :: Buffer(:,:)

public :: Buffer, Disk, Disk_1, Disk_2, DiskMX_Byte, iBuf, ID, InCore, IODone, ipos, iStatIO, lBuf, lDaRec, lStRec, LuTmp, &
          Mode_None, Mode_Read, Mode_Write, nSect, nBuf, OnDisk

end module IOBUF
