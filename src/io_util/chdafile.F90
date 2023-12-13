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

subroutine ChDaFile(Lu,iOpt,Buf,lBuf,iDisk)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf
integer(kind=iwp), intent(inout) :: iDisk
character, intent(inout) :: Buf(*)

call ChDaFile_Internal(Buf)

! This is to allow type punning without an explicit interface
contains

subroutine ChDaFile_Internal(Buf)
  character, target :: Buf(*)
  integer(kind=iwp), pointer :: iBuf(:)
  call c_f_pointer(c_loc(Buf(1)),iBuf,[1])
  call DaFile(Lu,iOpt,iBuf,lBuf,iDisk)
  nullify(iBuf)
end subroutine ChDaFile_Internal

end subroutine ChDaFile
