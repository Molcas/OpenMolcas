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

subroutine dBuu(uM12,nQQ,nDim,g,Hss,Inv)

use Slapaf_Info, only: dBM, idBM, mq, nqBM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nQQ, nDim
real(kind=wp), intent(in) :: uM12(nDim), g(nQQ)
real(kind=wp), intent(inout) :: Hss(nDim,nDim)
logical(kind=iwp), intent(in) :: Inv
integer(kind=iwp) :: i, i_Dim, idB, iElem, iq, iQQ, jDim, nElem
real(kind=wp) :: dBqR
real(kind=wp), allocatable :: K(:,:), Temp(:,:), Y(:)

if (.not. allocated(dBM)) then
  Hss(:,:) = Zero
  return
end if

call mma_allocate(Y,mq,Label='Y')
call mma_allocate(K,mq,nQQ,Label='K')
call Get_dArray('K',K,mq*nQQ)

! Compute Y(qR) = Sum_(Q) g(Q) rK(qR,Q)

Y(:) = Zero
do iQQ=1,nQQ
  Y(:) = Y(:)+g(iQQ)*K(:,iQQ)
end do
call mma_deallocate(K)

! Compute Temp = Sum_(qR) Y(qR) * dB(qR)

call mma_allocate(Temp,nDim,nDim,Label='Temp')
Temp(:,:) = Zero

idB = 1
do iq=1,mq
  nElem = nqBM(iq)
  do iElem=idB,idB+(nElem**2)-1
    dBqR = dBM(iElem)
    i_Dim = idBM(1+(iElem-1)*2)
    jDim = idBM(2+(iElem-1)*2)
    Temp(i_Dim,jDim) = Temp(i_Dim,jDim)+Y(iq)*dBqR
  end do
  idB = idB+nElem**2
end do
call mma_deallocate(Y)

if (Inv) Temp(:,:) = -Temp(:,:)

do i=1,nDim
  Hss(:,i) = Hss(:,i)+Temp(:,i)/sqrt(uM12(:)*uM12(i))
end do

call mma_deallocate(Temp)

return

end subroutine dBuu
