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

subroutine dDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Fast_IO, only: MBL
use Definitions, only: wp, iwp, RtoB

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf_
real(kind=wp), intent(inout) :: Buf(lBuf_)
integer(kind=iwp), intent(inout) :: iDisk_
integer(kind=iwp) :: lBuf, iDisk

call dDaFile_Internal(Buf)

! This is to allow type punning without an explicit interface
contains

subroutine dDaFile_Internal(Buf)
  real(kind=wp), target :: Buf(*)
  character, pointer :: cBuf(:)

  lBuf = lBuf_*RtoB
  iDisk = iDisk_*MBL(Lu)

  call c_f_pointer(c_loc(Buf(1)),cBuf,[lBuf])
  call bDaFile(Lu,iOpt,cBuf,lBuf,iDisk)
  nullify(cBuf)

  iDisk_ = (iDisk+MBL(Lu)-1)/MBL(Lu)

  return
end subroutine dDaFile_Internal

end subroutine dDaFile
