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
!     lDaRec  : minimal record (in 4Byte words) transfered by         *
!               DAFILE                                                *
!     nSect   : integer multiple of lDaRec which consitutes the       *
!               standard record length                                *
!     lStRec  : standard record length in units of integers           *
!     nDiskMx : Max disk space in units of MByte                      *
!     iDiskMx : Max disk space -1 in units of Byte                    *
!                                                                     *
!---------------------------------------------------------------------*

module IOBUF

real*8 DiskMx_MByte, DiskMx_Byte
parameter(lDaRec=2**10) ! 1024
parameter(nSect=2**5) ! 32
parameter(lStRec=nSect*lDaRec)
parameter(Mode_Read=987654321)
parameter(Mode_Write=198765432)
parameter(Mode_None=111111111)
!parameter (DiskMx_MByte=2.0D00**11)
!parameter (DiskMx_Byte =2.0D00**31)
!real*8 Buf_IO(lStRec,2)
integer iPos, LuTmp, iStatIO, id, iBuf, ipBuf, nBuf, lBuf
real*8 Disk, Disk_1, Disk_2
logical IODone, InCore, OnDisk
real*8, allocatable :: Buffer(:,:)

end module IOBUF
