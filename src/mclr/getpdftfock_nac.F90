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

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMat, nDens
use input_mclr, only: nBas, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: bk(nDens)
integer(kind=iwp) :: IS, JS
real(kind=wp), allocatable :: T(:), bktmp(:)

call mma_allocate(bktmp,nDens)
call mma_allocate(T,nDens)
call Get_DArray('FxyMS',T,nDens)

do IS=1,nSym
  jS = Mul(iS,1)
  if (nBas(iS)*nBas(jS) /= 0) &
    call DGeSub(T(ipMat(iS,jS)),nBas(iS),'N',T(ipMat(jS,iS)),nBas(jS),'T',bktmp(ipMat(iS,jS)),nBas(iS),nBas(iS),nBas(jS))
end do
bk(:) = bk(:)-Two*bktmp(:)
call mma_deallocate(T)
call mma_deallocate(bktmp)

end subroutine GetPDFTFock_NAC
