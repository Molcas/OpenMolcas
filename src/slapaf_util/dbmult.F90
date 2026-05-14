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

use Slapaf_Info, only: dBM, idBM, mq, nqBM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nQQ, nDim, nLambda
real(kind=wp), intent(in) :: dCdQ(nQQ,nLambda)
real(kind=wp), intent(out) :: QC(nDim**2,nLambda)
integer(kind=iwp) :: i_Dim, idB, iElem, ijDim, iq, jdim, nElem
real(kind=wp) :: dBqR
real(kind=wp), allocatable :: X(:,:), K(:,:)

QC(:,:) = Zero

if (.not. allocated(dBM)) then
  !write(u6,*) 'FAST out'
  return
end if

call mma_allocate(X,mq,nLambda,Label='X')
X(:,:) = Zero
call mma_allocate(K,mq,nQQ,Label='K')
call Get_dArray('K',K,mq*nQQ)

call DGEMM_('N','N',mq,nLambda,nQQ,One,K,mq,dCdQ,nQQ,Zero,X,mq)
call mma_deallocate(K)

idB = 1
do iq=1,mq
  nElem = nqBM(iq)
  do iElem=idB,idB+(nElem**2)-1
    dBqR = dBM(iElem)
    i_Dim = idBM(1+(iElem-1)*2)
    jDim = idBM(2+(iElem-1)*2)
    ijDim = (jDim-1)*nDim+i_Dim
    QC(ijDim,:) = QC(ijDim,:)+X(iq,:)*dBqR
  end do
  idB = idB+nElem**2
end do
call mma_deallocate(X)

return

end subroutine dBMult
