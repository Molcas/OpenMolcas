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

subroutine SOAO_g(iSD4,nSD,nSO,MemPrm,MemMax,iBsInc,jBsInc,kBsInc,lBsInc,iPrInc,jPrInc,kPrInc,lPrInc,ipMem1,ipMem2,Mem1,Mem2,iFnc, &
                  MemPSO)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4), nSO, MemPrm, MemMax, ipMem1
integer(kind=iwp), intent(out) :: iBsInc, jBsInc, kBsInc, lBsInc, iPrInc, jPrInc, kPrInc, lPrInc, ipMem2, Mem1, Mem2, iFnc(4), &
                                  MemPSO

call PSOAO1(nSO,MemPrm,MemMax,iFnc,iBsInc,jBsInc,kBsInc,lBsInc,iPrInc, &
            jPrInc,kPrInc,lPrInc,ipMem1,ipMem2,Mem1,Mem2,MemPSO,nSD,iSD4)

end subroutine SOAO_g
