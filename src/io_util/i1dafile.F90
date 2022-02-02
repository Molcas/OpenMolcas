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
! like cdafile, but for single-byte integers

subroutine i1DaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Fast_IO, only: MBL
use Definitions, only: iwp, byte

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf_
integer(kind=byte), intent(inout) :: Buf(lBuf_)
integer(kind=iwp), intent(inout) :: iDisk_
integer(kind=iwp) :: lBuf, iDisk

call i1DaFile_Internal(Buf)

! This is to allow type punning without an explicit interface
contains

subroutine i1DaFile_Internal(Buf)

  integer(kind=byte), target, intent(inout) :: Buf(lBuf_)
  character, pointer :: cBuf(:)

  lBuf = lBuf_
  iDisk = iDisk_*MBL(Lu)

  call c_f_pointer(c_loc(Buf(1)),cBuf,[lBuf_])
  call bDaFile(Lu,iOpt,cBuf,lBuf,iDisk)
  nullify(cBuf)

  iDisk_ = (iDisk+MBL(Lu)-1)/MBL(Lu)

  return

end subroutine i1DaFile_Internal

end subroutine i1DaFile
