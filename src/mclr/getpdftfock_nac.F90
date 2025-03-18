!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Paul B Calio                                     *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Based on cmsbk.f from Jie J. Bao &&                            *
! Additional work from  rhs_nac.f                                *
! ****************************************************************
      Subroutine GetPDFTFock_NAC(bk)
      use stdalloc, only : mma_allocate, mma_deallocate
      use MCLR_Data, only: nDens2, ipMat
      use input_mclr, only: nSym,nBas
      Implicit None
!*****Output
      Real*8,DIMENSION(nDens2)::bk
!*****Input
!*****Auxiliaries
      Real*8,DIMENSION(:),Allocatable::T,FT99,bktmp
      INTEGER IS,JS
      CALL mma_allocate(FT99,nDens2)
      CALL mma_allocate(bktmp,nDens2)
      CALL mma_allocate(T,nDens2)
      CALL Get_DArray('FxyMS           ',FT99 ,nDens2)
      CALL dcopy_(nDens2,FT99,1,T,1)

      DO IS=1,nSym
         jS=iEOR(iS-1,0)+1
         If (nBas(is)*nBas(jS).ne.0) then
           Call DGeSub(T(ipMat(iS,jS)),nBas(iS),'N',                    &
     &                 T(ipMat(jS,iS)),nBas(jS),'T',                    &
     &                 bktmp(ipMat(iS,jS)),nBas(iS),                    &
     &                 nBas(iS),nBas(jS))
         End If
      END DO
      CALL daxpy_(nDens2,-2.0d0,bktmp,1,bk,1)
      CALL mma_deallocate(T)
      CALL mma_deallocate(FT99)
      CALL mma_deallocate(bktmp)
      end subroutine GetPDFTFock_NAC
