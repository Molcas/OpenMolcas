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

Real*8, Allocatable, Target:: CHVBUF(:)
Real*8, Allocatable, Target:: CHVBFI(:)

INTEGER ip_CHVBUF_SYM(8), l_CHVBUF_SYM(8)
INTEGER ip_CHVBFI_SYM(8), l_CHVBFI_SYM(8)
INTEGER nVec_in_Buf(8)
End Module ChoVecBuf
