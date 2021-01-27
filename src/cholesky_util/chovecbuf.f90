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
!
! Information about the Cholesky vector buffer
!
Module ChoVecBuf

Type rPointers
  Real*8, Pointer:: A(:,:)
End Type rPointers

Type (rPointers):: CHVBUF(8)
Real*8, Allocatable, Target:: CHVBUF_T(:)

INTEGER ip_CHVBFI, l_CHVBFI
INTEGER ip_CHVBUF_SYM, l_CHVBUF_SYM
INTEGER ip_CHVBFI_SYM, l_CHVBFI_SYM

COMMON / CHVBUF / ip_CHVBFI, l_CHVBFI, &
                  ip_CHVBUF_SYM(8), l_CHVBUF_SYM(8), &
                  ip_CHVBFI_SYM(8), l_CHVBFI_SYM(8), &
                  NVEC_IN_BUF(8)
Save
End Module ChoVecBuf
