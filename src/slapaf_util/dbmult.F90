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

subroutine dBMult(dCdQ,QC,nQQ,nDim,nLambda)

use Slapaf_info, only: dBM, idBM, nqBM
use Slapaf_parameters, only: mq

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
real*8 dCdQ(nQQ,nLambda), QC(nDim**2,nLambda)
real*8, allocatable :: X(:,:), K(:,:)

QC(:,:) = 0.0d0

if (.not. allocated(dBM)) then
  !write(6,*) 'FAST out'
  return
end if

call mma_allocate(X,mq,nLambda,Label='X')
X(:,:) = 0.0d0
call mma_allocate(K,mq,nQQ,Label='K')
call Get_dArray('K',K,mq*nQQ)

call DGEMM_('N','N',mq,nLambda,nQQ,1.0d0,K,mq,dCdQ,nQQ,0.0d0,X,mq)
call mma_deallocate(K)

idB = 1
do iq=1,mq
  nElem = nqBM(iq)
  do iElem=idB,idB+(nElem**2)-1
    dBqR = dBM(iElem)
    iDim = idBM(1+(iElem-1)*2)
    jDim = idBM(2+(iElem-1)*2)
    ijDim = (jDim-1)*nDim+iDim
    do iLambda=1,nLambda
      QC(ijDim,iLambda) = QC(ijDim,iLambda)+X(iq,iLambda)*dBqR
    end do
  end do
  idB = idB+nElem**2
end do
call mma_deallocate(X)

return

end subroutine dBMult
