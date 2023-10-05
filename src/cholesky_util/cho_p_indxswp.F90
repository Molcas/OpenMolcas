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

subroutine Cho_P_IndxSwp()
!
! Purpose: swap global and local reduced set index arrays.
!          Note that arrays InfRed and InfVec are swapped as well.
!
! NB: this procedure is inexpensive, as we are merely swapping
!     pointers, not actual data (except for the statically allocated
!     index arrays which amount to swapping 51 integers in total).

use Cholesky, only: iiBstR, iiBstR_G, iiBstRSh, iiBstRSh_G, IndRed, IndRed_G, IndRSh, IndRSh_G, InfRed, InfRed_G, InfVec, &
                    InfVec_G, mmBstRT, mmBstRT_G, nnBstR, nnBstR_G, nnBstRSh, nnBstRSh_G, nnBstRT, nnBstRT_G, nnShl, nnShl_G, &
                    pTemp, pTemp1, pTemp3
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iTmp, N

iTmp = nnShl_G
nnShl_G = nnShl
nnShl = iTmp

iTmp = mmBstRT_G
mmBstRT_G = mmBstRT
mmBstRT = iTmp

N = 8*3
call iSwap(N,iiBstR_G,1,iiBstR,1)
call iSwap(N,nnBstR_G,1,nnBstR,1)
call iSwap(3,nnBstRT_G,1,nnBstRT,1)

pTemp1 => InfRed_G
InfRed_G => InfRed
InfRed => pTemp1

pTemp3 => InfVec_G
InfVec_G => InfVec
InfVec => pTemp3

pTemp3 => iiBstRSh_G
iiBstRSh_G => iiBstRSh
iiBstRSh => pTemp3

pTemp3 => nnBstRSh_G
nnBstRSh_G => nnBstRSh
nnBstRSh => pTemp3

pTemp => IndRed_G
IndRed_G => IndRed
IndRed => pTemp

pTemp1 => IndRSh_G
IndRSh_G => IndRSh
IndRSh => pTemp1

end subroutine Cho_P_IndxSwp
