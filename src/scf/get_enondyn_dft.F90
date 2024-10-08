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

use Index_Functions, only: iTri, nTri_Elem
use InfSCF, only: CMO, Erest_xc, KSDFT, nBas, nBT, nOcc, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nh1, nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
character(len=4), intent(in) :: DFTFOCK
integer(kind=iwp) :: i, iDji, iOff, iSym, j, ji, jOff
real(kind=wp), allocatable :: D_DS(:,:), F_DFT(:,:)

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
      ji = iTri(j,i)
      iDji = iOff-1+ji
      D_DS(iDji,1) = Two*D_DS(iDji,1)
      D_DS(iDji,2) = Two*D_DS(iDji,2)
    end do
  end do
  iOff = iOff+nBas(iSym)*nOrb(iSym)
  jOff = jOff+nTri_Elem(nBas(iSym))
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
