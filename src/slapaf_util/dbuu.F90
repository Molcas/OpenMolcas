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

use Slapaf_Info, only: dBM, idBM, nqBM
use Slapaf_parameters, only: mq

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
real*8 uM12(nDim), g(nQQ), Hss(nDim,nDim)
logical Inv
real*8, allocatable :: Y(:), K(:,:), Temp(:,:)

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
  call DaXpY_(mq,g(iQQ),K(:,iQQ),1,Y,1)
end do
call mma_deallocate(K)

! Compute Temp = Sum_(qR) Y(qR) * dB(qR)

call mma_allocate(Temp,nDim,nDim,Label='Temp')
Temp(:,:) = Zero

idB = 1
do iq=1,mq
  YqR = Y(iq)
  nElem = nqBM(iq)
  do iElem=idB,idB+(nElem**2)-1
    dBqR = dBM(iElem)
    iDim = idBM(1+(iElem-1)*2)
    jDim = idBM(2+(iElem-1)*2)
    Temp(iDim,jDim) = Temp(iDim,jDim)+YqR*dBqR
  end do
  idB = idB+nElem**2
end do
call mma_deallocate(Y)

if (Inv) call DScal_(nDim**2,-One,Temp,1)

do i=1,nDim
  do j=1,nDim
    xx = sqrt(uM12(i)*uM12(j))
    Hss(i,j) = Hss(i,j)+Temp(i,j)/xx
  end do
end do

call mma_deallocate(Temp)

return

end subroutine dBuu
