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

subroutine Cho_P_SyncNumCho_Time(c,w)

implicit none
real*8 c, w
#include "cholesky.fh"

tMisc(1,5) = tMisc(1,5)+c
tMisc(2,5) = tMisc(2,5)+w

end subroutine Cho_P_SyncNumCho_Time
