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
CSVC: modified to convert to the use of byte lengths/offests by the
C     underlying I/O routines (2016)
      Subroutine dDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

      Implicit None

#include "SysDef.fh"
#include "fio.fh"

      Integer Lu, iOpt, lBuf_, iDisk_
      Real*8 Buf(lBuf_)

      Integer lBuf, iDisk

      lBuf=lBuf_*RtoB
      iDisk=iDisk_*MBL(Lu)

      Call bDaFile(Lu,iOpt,Buf,lBuf,iDisk)


      iDisk_=(iDisk+MBL(Lu)-1)/MBL(Lu)

      Return
      End


      Subroutine cDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

      Implicit None

#include "SysDef.fh"
#include "fio.fh"

      Integer Lu, iOpt, lBuf_, iDisk_
      Character*1 Buf(lBuf_)

      Integer lBuf, iDisk

      lBuf=lBuf_
      iDisk=iDisk_*MBL(Lu)

      Call bDaFile(Lu,iOpt,Buf,lBuf,iDisk)

      iDisk_=(iDisk+MBL(Lu)-1)/MBL(Lu)

      Return
      End


      Subroutine iDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

      Implicit None

#include "SysDef.fh"
#include "fio.fh"

      Integer Lu, iOpt, lBuf_, iDisk_
      Integer Buf(lBuf_)

      Integer lBuf, iDisk

      lBuf=lBuf_*ItoB
      iDisk=iDisk_*MBL(Lu)

      Call bDaFile(Lu,iOpt,Buf,lBuf,iDisk)

      iDisk_=(iDisk+MBL(Lu)-1)/MBL(Lu)

      Return
      End
