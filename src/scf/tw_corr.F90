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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

subroutine Tw_corr(irc,DeTW,CMOI,EOcc,EVir)

use InfSCF, only: nBas, nBT, nDel, nFro, nOcc, nSym
use ChoMP2, only: ChoAlg, DoDens
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iRC
real(kind=wp) :: DeTW, CMOI(*), EOcc(*), EVir(*)
integer(kind=iwp) :: i, nElk, nExt(8)
real(kind=wp) :: Grad(1), TW, TW0
real(kind=wp), allocatable :: DMAT(:,:), F_DFT(:)

DoDens = .false.
ChoAlg = 2

call mma_allocate(DMAT,nBT,2,Label='DMAT')

nElk = 0
do i=1,nSym
  nExt(i) = nBas(i)-nDel(i)-nOcc(i,1)-nFro(i)
  nElk = nElk+2*(nFro(i)+nOcc(i,1))
end do

call DM_FNO_RHF(irc,nSym,nBas,nFro,nOcc(1,1),nExt,nDel,CMOI,EOcc,EVir,DMAT(:,2),DMAT(:,1))
if (irc /= 0) then
  write(u6,*) 'DM_FNO_RHF returned ',irc
  call SysAbendMsg('DM_FNO_RHF','Non-zero return code from DM_FNO_RHF',' ')
end if

call mma_allocate(F_DFT,nBT,Label='F_DFT')

call Fold_tMat(nSym,nBas,DMAT(:,1),DMAT(:,1))
call dscal_(nBT,Half,DMAT(:,1),1)
call Fold_tMat(nSym,nBas,DMAT(:,2),DMAT(:,2))
call dscal_(nBT,Half,DMAT(:,2),1)
Grad = Zero

call wrap_DrvNQ('HUNTER',F_DFT,1,TW,DMAT(:,1),nBT,1,.false.,Grad,1,'SCF ')

call wrap_DrvNQ('HUNTER',F_DFT,1,TW0,DMAT(:,2),nBT,1,.false.,Grad,1,'SCF ')
DeTW = (TW-TW0)/real(nElk,kind=wp)

call mma_deallocate(F_DFT)
call mma_deallocate(DMAT)

return

end subroutine Tw_corr
