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

use Index_Functions, only: iTri, nTri_Elem
use SpinAV, only: Do_SpinAV, DSC
use InfSCF, only: addc_KSDFT, CMO, DE_KSDFT_c, nBas, nBT, nConstr, nOcc, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nh1, nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
character(len=4), intent(in) :: DFTFOCK
integer(kind=iwp) :: i, iAB, iDij, iDSC, iOff, iSym, iXoX, j, ji, jOff, lOff, mAdCMOO, nXoX
real(kind=wp) :: Ec_AB(2)
real(kind=wp), allocatable :: D_DS(:,:), F_DFT(:,:)
integer(kind=iwp), parameter :: nD = 2

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
          ji = iTri(j,i)
          iDij = jOff-1+ji
          D_DS(iDij,1) = D_DS(iDij,1)-DSc(iDSc)
          D_DS(iDij,2) = D_DS(iDij,2)+DSc(iDSc)
        end do
      end do
      lOff = lOff+nBas(iSym)**2
    end if

    do j=1,nBas(iSym)
      do i=1,j-1
        ji = iTri(j,i)
        iDij = jOff-1+ji
        D_DS(iDij,1) = Two*D_DS(iDij,1)
        D_DS(iDij,2) = Two*D_DS(iDij,2)
      end do
    end do
    iOff = iOff+nBas(iSym)*nOrb(iSym)
    jOff = jOff+nTri_Elem(nBas(iSym))
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
