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
! Based on cmsbk.f from Jie J. Bao                               *
! Additional work from  rhs_nac                                  *
! ****************************************************************

subroutine GetPDFTFock_NAC(bk)

use MCLR_Data, only: nDens2, ipMat
use input_mclr, only: nSym, nBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two

implicit none
! Output
real*8, dimension(nDens2) :: bk
! Auxiliaries
real*8, dimension(:), allocatable :: T, FT99, bktmp
integer IS, JS

call mma_allocate(FT99,nDens2)
call mma_allocate(bktmp,nDens2)
call mma_allocate(T,nDens2)
call Get_DArray('FxyMS',FT99,nDens2)
call dcopy_(nDens2,FT99,1,T,1)

do IS=1,nSym
  jS = ieor(iS-1,0)+1
  if (nBas(is)*nBas(jS) /= 0) &
    call DGeSub(T(ipMat(iS,jS)),nBas(iS),'N',T(ipMat(jS,iS)),nBas(jS),'T',bktmp(ipMat(iS,jS)),nBas(iS),nBas(iS),nBas(jS))
end do
call daxpy_(nDens2,-Two,bktmp,1,bk,1)
call mma_deallocate(T)
call mma_deallocate(FT99)
call mma_deallocate(bktmp)

end subroutine GetPDFTFock_NAC
