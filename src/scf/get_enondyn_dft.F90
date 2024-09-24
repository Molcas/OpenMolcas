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

subroutine Get_Enondyn_dft(nh1,Grad,nGrad,DFTFOCK)

use SCF_Arrays, only: CMO
use InfSCF, only: KSDFT, nBT, nSym, nBas, nOrb, nOcc
use DCSCF, only: Erest_xc
use Constants, only: Zero, One, Two
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer nh1, nGrad
real*8 Grad(nGrad)
character(len=4) DFTFOCK
integer i, iDji, iSym, j, ji, iOff, jOff
real*8, allocatable :: F_DFT(:,:), D_DS(:,:)

Erest_xc = Zero
call mma_allocate(D_DS,nBT,2,Label='D_DS ')
D_DS(:,:) = Zero

iOff = 1
jOff = 1
do iSym=1,nSym
  call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nOcc(iSym,1), &
                 One,CMO(iOff,1),nBas(iSym), &
                 CMO(iOff,1),nBas(iSym), &
                 Zero,D_DS(jOff:,1),nBas(iSym))
  call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nOcc(iSym,2), &
                 One,CMO(iOff,2),nBas(iSym), &
                 CMO(iOff,2),nBas(iSym), &
                 Zero,D_DS(jOff:,2),nBas(iSym))
  do j=1,nBas(iSym)
    do i=1,j-1
      ji = j*(j-1)/2+i
      iDji = iOff-1+ji
      D_DS(iDji,1) = Two*D_DS(iDji,1)
      D_DS(iDji,2) = Two*D_DS(iDji,2)
    end do
  end do
  iOff = iOff+nBas(iSym)*nOrb(iSym)
  jOff = jOff+nBas(iSym)*(nBas(iSym)+1)/2
end do

!----------------------------------------------------------------------*
call Get_Fmat_nondyn(D_DS(:,1),D_DS(:,2),nBT,.true.)
!----------------------------------------------------------------------*

!----------------------------------------------------------------------*
call mma_allocate(F_DFT,nBT,2,Label='F_DFT')
call Get_Exc_dft(nh1,Grad,nGrad,DFTFOCK,F_DFT,D_DS,nBT,2,KSDFT)
!----------------------------------------------------------------------*

call mma_deallocate(D_DS)
call mma_deallocate(F_DFT)

return

end subroutine Get_Enondyn_dft
