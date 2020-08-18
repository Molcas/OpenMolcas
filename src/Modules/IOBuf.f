************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
*---------------------------------------------------------------------*
*                                                                     *
*     lDaRec  : minimal record (in 4Byte words) transfered by         *
*               DAFILE                                                *
*     nSect   : integer multiple of lDaRec which consitutes the       *
*               standard record length                                *
*     lStRec  : standard record length in units of integers           *
*     nDiskMx : Max disk space in units of MByte                      *
*     iDiskMx : Max disk space -1 in units of Byte                    *
*                                                                     *
*---------------------------------------------------------------------*
      Module IOBUF
*
      Real*8 DiskMx_MByte, DiskMx_Byte
      Parameter ( lDaRec  = 2**10        ) ! 1024
      Parameter ( nSect   = 2**5         ) ! 32
      Parameter ( lStRec  = nSect*lDaRec )
      Parameter ( Mode_Read = 987654321  )
      Parameter ( Mode_Write= 198765432  )
      Parameter ( Mode_None = 111111111  )
C     Parameter ( DiskMx_MByte=2.0D00**11)
C     Parameter ( DiskMx_Byte =2.0D00**31)
*     Real*8 Buf_IO(lStRec,2)
      Integer iPos, LuTmp, iStatIO, id, iBuf, ipBuf, nBuf, lBuf
      Real*8 Disk, Disk_1, Disk_2
      Logical IODone,InCore,OnDisk
      Real*8, Allocatable:: Buffer(:,:)
*
      End Module IOBUF
