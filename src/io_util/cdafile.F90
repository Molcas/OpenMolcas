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
!SVC: modified to convert to the use of byte lengths/offests by the
!     underlying I/O routines (2016)

subroutine cDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

implicit none

#include "SysDef.fh"
#include "fio.fh"
integer Lu, iOpt, lBuf_, iDisk_
character*1 Buf(lBuf_)
integer lBuf, iDisk

lBuf = lBuf_
iDisk = iDisk_*MBL(Lu)

call bDaFile(Lu,iOpt,Buf,lBuf,iDisk)

iDisk_ = (iDisk+MBL(Lu)-1)/MBL(Lu)

return

end subroutine cDaFile
