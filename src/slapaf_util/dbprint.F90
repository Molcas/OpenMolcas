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

subroutine dBPrint(nQQ,nDim)

use Slapaf_Info, only: dBM, idBM, nqBM
use Slapaf_parameters, only: mq

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
real*8, allocatable :: dBQQ(:,:), K(:,:)

if (.not. allocated(dBM)) return
call mma_allocate(dBQQ,nDim,nDim,Label='dBQQ')
call mma_allocate(K,mq,nQQ,Label='K')
call Get_dArray('K',K,mq*nQQ)
do iQQ=1,nQQ
  dBQQ(:,:) = Zero
  idB = 1
  do iq=1,mq
    nElem = nqBM(iq)
    rK = K(iq,iQQ)
    do iElem=idB,idB+(nElem**2)-1
      dBqR = dBM(iElem)
      iDim = idBM(1+(iElem-1)*2)
      jDim = idBM(2+(iElem-1)*2)
      dBQQ(iDim,jDim) = dBQQ(iDim,jDim)+rK*dBqR
    end do
    idB = idB+nElem**2
  end do
  call RecPrt('dBQQ',' ',dBQQ,nDim,nDim)
end do
call mma_deallocate(dBQQ)
call mma_deallocate(K)

return

end subroutine dBPrint
