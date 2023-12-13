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
!SVC: modified to convert to the use of byte lengths/offsets by the
!     underlying I/O routines (2016)

subroutine cDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

use Fast_IO, only: MBL
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf_
character, intent(inout) :: Buf(lBuf_)
integer(kind=iwp), intent(inout) :: iDisk_
integer(kind=iwp) :: lBuf, iDisk

lBuf = lBuf_
iDisk = iDisk_*MBL(Lu)

call bDaFile(Lu,iOpt,Buf,lBuf,iDisk)

iDisk_ = (iDisk+MBL(Lu)-1)/MBL(Lu)

return

end subroutine cDaFile
