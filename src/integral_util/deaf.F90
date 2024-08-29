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

! This module simply provides an interface for the EAF subroutines with
! real buffers instead of integers, note that nBuf is still the length
! of the integer buffer.

module dEAF

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp

implicit none
private

public :: dEAFARead, dEAFAWrite, dEAFRead, dEAFWrite

contains

subroutine dEAFARead(Lu,Buf,nBuf,Disk,id)

  integer(kind=iwp), intent(in) :: Lu, nBuf
  real(kind=wp), target, intent(out) :: Buf(nBuf)
  real(kind=wp), intent(inout) :: Disk
  integer(kind=iwp), intent(out) :: id
  integer(kind=iwp), pointer :: iBuf(:)

  call c_f_pointer(c_loc(Buf),iBuf,[nBuf])
  call EAFARead(Lu,iBuf,nBuf,Disk,id)
  nullify(iBuf)

end subroutine dEAFARead

subroutine dEAFAWrite(Lu,Buf,nBuf,Disk,id)

  integer(kind=iwp), intent(in) :: Lu, nBuf
  real(kind=wp), target, intent(in) :: Buf(nBuf)
  real(kind=wp), intent(inout) :: Disk
  integer(kind=iwp), intent(out) :: id
  integer(kind=iwp), pointer :: iBuf(:)

  call c_f_pointer(c_loc(Buf),iBuf,[nBuf])
  call EAFAWrite(Lu,iBuf,nBuf,Disk,id)
  nullify(iBuf)

end subroutine dEAFAWrite

subroutine dEAFRead(Lu,Buf,nBuf,Disk)

  integer(kind=iwp), intent(in) :: Lu, nBuf
  real(kind=wp), target, intent(out) :: Buf(nBuf)
  real(kind=wp), intent(inout) :: Disk
  integer(kind=iwp), pointer :: iBuf(:)

  call c_f_pointer(c_loc(Buf),iBuf,[nBuf])
  call EAFRead(Lu,iBuf,nBuf,Disk)
  nullify(iBuf)

end subroutine dEAFRead

subroutine dEAFWrite(Lu,Buf,nBuf,Disk)

  integer(kind=iwp), intent(in) :: Lu, nBuf
  real(kind=wp), target, intent(in) :: Buf(nBuf)
  real(kind=wp), intent(inout) :: Disk
  integer(kind=iwp), pointer :: iBuf(:)

  call c_f_pointer(c_loc(Buf),iBuf,[nBuf])
  call EAFWrite(Lu,iBuf,nBuf,Disk)
  nullify(iBuf)

end subroutine dEAFWrite

end module dEAF
