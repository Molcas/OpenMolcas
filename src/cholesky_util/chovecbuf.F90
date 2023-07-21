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

! Information about the Cholesky vector buffer
module ChoVecBuf

implicit none
private

public :: CHVBUF, ip_CHVBUF_SYM, l_CHVBUF_SYM
public :: CHVBFI, ip_CHVBFI_SYM, l_CHVBFI_SYM
public :: nVec_in_Buf
real*8, allocatable, target :: CHVBUF(:)
real*8, allocatable, target :: CHVBFI(:)
integer ip_CHVBUF_SYM(8), l_CHVBUF_SYM(8)
integer ip_CHVBFI_SYM(8), l_CHVBFI_SYM(8)
integer nVec_in_Buf(8)

end module ChoVecBuf
