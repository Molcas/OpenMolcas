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

subroutine datimx(TimeStamp)

implicit none
character(len=*), intent(out) :: TimeStamp
interface
  subroutine datimxc(TimeStamp) bind(C,name='datimxc_')
    use, intrinsic :: iso_c_binding, only: c_char
    character(kind=c_char) :: TimeStamp(*)
  end subroutine datimxc
end interface

call datimxc(TimeStamp)

end subroutine datimx
