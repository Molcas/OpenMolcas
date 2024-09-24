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

subroutine Get_DEcorr(nh1,Grad,nGrad,DFTFOCK)

use SCF_Arrays, only: CMO
use SpinAV, only: Do_SpinAV, DSC
use InfSCF, only: nBT, nSym, nOcc, nConstr, nBas, nOrb
use Constants, only: Zero, One, Two
use AddCorr, only: addc_KSDFT, DE_KSDFT_c
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer nh1, nGrad
real*8 Grad(nGrad), Ec_AB(2)
character(LEN=4) DFTFOCK
integer i, iAB, iDij, iDSC, iOff, iSym, iXoX, nXoX, j, jOff, ji, lOff, mAdCMOO, Nd
real*8, allocatable :: F_DFT(:,:), D_DS(:,:)

nD = 2

call mma_allocate(F_DFT,nBT,nD,Label='F_DFT')
call mma_allocate(D_DS,nBT,nD,Label='D_DS')

do iAB=1,2
  iOff = 1
  jOff = 1
  lOff = 0
  do iSym=1,nSym
    if (iAB == 1) then
      nXoX = nOcc(iSym,1)
      iXoX = 0
    else
      nXoX = nConstr(iSym)
      iXoX = nOcc(iSym,1)-nConstr(iSym)
    end if
    mAdCMOO = iOff+nBas(iSym)*iXoX
    call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nXoX, &
                   One,CMO(mAdCMOO,1),nBas(iSym), &
                   CMO(mAdCMOO,1),nBas(iSym), &
                   Zero,D_DS(jOff,1),nBas(iSym))
    if (iAB == 1) then
      nXoX = nOcc(iSym,2)
      iXoX = 0
    else
      nXoX = nConstr(iSym)
      iXoX = nOcc(iSym,2)-nConstr(iSym)
    end if
    mAdCMOO = iOff+nBas(iSym)*iXoX
    call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nXoX, &
                   One,CMO(mAdCMOO,2),nBas(iSym), &
                   CMO(mAdCMOO,2),nBas(iSym), &
                   Zero,D_DS(jOff,2),nBas(iSym))

    if (Do_SpinAV) then
      do j=1,nBas(iSym)
        do i=1,j
          iDSc = nBas(iSym)*(j-1)+i
          ji = j*(j-1)/2+i
          iDij = jOff-1+ji
          D_DS(iDij,1) = D_DS(iDij,1)-DSc(iDSc)
          D_DS(iDij,2) = D_DS(iDij,2)+DSc(iDSc)
        end do
      end do
      lOff = lOff+nBas(iSym)**2
    end if

    do j=1,nBas(iSym)
      do i=1,j-1
        ji = j*(j-1)/2+i
        iDij = jOff-1+ji
        D_DS(iDij,1) = Two*D_DS(iDij,1)
        D_DS(iDij,2) = Two*D_DS(iDij,2)
      end do
    end do
    iOff = iOff+nBas(iSym)*nOrb(iSym)
    jOff = jOff+nBas(iSym)*(nBas(iSym)+1)/2
  end do

  call Get_Ecorr_dft(nh1,Grad,nGrad,DFTFOCK,F_DFT,D_DS,nBT,nD,ADDC_KSDFT,Ec_AB(iAB))
end do
!----------------------------------------------------------------------*
DE_KSDFT_c = Ec_AB(1)-Ec_AB(2)
!----------------------------------------------------------------------*

call mma_deallocate(D_DS)
call mma_deallocate(F_DFT)

return

end subroutine Get_DEcorr
