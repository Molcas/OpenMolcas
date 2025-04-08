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

subroutine Start_MCLR()
!***********************************************************************
!                                                                      *
!     Precompute whatever can be before starting the response section  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use transform_procedures, only: SetUp_CASPT2_Tra
use MCLR_Data, only: CMO_Inv, CMO
use MCLR_Data, only: nDens
use MCLR_Data, only: LuTri1, LuMotra, FnTri1, FnMotra, FnQDat, LuHlf2, LuHlf3, LuQDat, LuTri2
use input_mclr, only: StepType, TwoStep, NewCho, nSym, kPrint, nAsh, nBas, nDel, LuAChoVec, LuChoInt, LuIChoVec, nFro, nIsh, nOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, u6

implicit none
#include "warnings.h"
character(len=8) :: Label
character(len=5) Fname
real*8, allocatable :: STmat(:), Smat(:)
integer lTriDens, lSqrDens, nOrbBas, iSym, iSymLbl, iOpt, iComp, Index, iOff, i, j, iOff1, iOff2, iType, iSeed, iRC
integer, external :: IsFreeUnit
real*8 BufFrac

!----------------------------------------------------------------------*
!     start                                                            *
!----------------------------------------------------------------------*

!                                                                      *
!***********************************************************************
!                                                                      *
call setup_MCLR(1)

if ((StepType /= 'RUN2') .and. btest(kPrint,2)) write(u6,'(6X,A)') 'Transformation of integrals'
! For the mp2-gradient calculations we want the transformation
! routine to produce all integrals of the occupied and virtual
! orbitals so we tell it that the whole space is inactive and
! that the active is empty

! Use the driver from the CASPT2 code (only none-squared).

! re-direct the transformed integrals to the MOTRA file
! which is preserved at the end of the calculation.
! LuTri1 is deleted.
if (TwoStep) then
  LuTri1 = LuMOTRA
  FnTri1 = FnMOTRA
  call DaName_MF_wa(LuQDAT,FnQDAT)
end if

call DaName_MF_wa(LuTri1,FnTri1)

if (newCho) then

  ! Compute inverse CMO

  lTriDens = 0
  lSqrDens = 0
  nOrbBas = 0
  do iSym=1,nSym
    lTriDens = lTriDens+nTri_Elem(nBas(iSym))
    lSqrDens = lSqrDens+nBas(iSym)**2
    nOrbBas = nOrbBas+nOrb(iSym)*nBas(iSym)
  end do
  call mma_allocate(STmat,lTriDens,Label='STmat')
  call mma_allocate(Smat,lSqrDens,Label='Smat')

  iSymlbl = 1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  Label = 'Mltpl  0'
  iComp = 1
  call RdOne(irc,iOpt,Label,iComp,STmat,iSymlbl)

  index = 1
  iOff = 0
  do iSym=1,nSym
    do i=1,nBas(iSym)
      do j=1,i-1
        Smat(j+nBas(iSym)*(i-1)+iOff) = STmat(index)
        Smat(i+nBas(iSym)*(j-1)+iOff) = STmat(index)
        index = index+1
      end do
      Smat(i+nBas(iSym)*(i-1)+iOff) = STmat(index)
      index = index+1
    end do
    ioff = ioff+nBas(iSym)**2
  end do
  call mma_deallocate(STmat)

  call mma_allocate(CMO_inv,nOrbBas,Label='CMO_Inv')
  iOff1 = 1
  iOff2 = 1
  do iSym=1,nSym
    call dGemm_('T','N',nOrb(iSym),nBas(iSym),nBas(iSym),One,CMO(iOff2),nBas(iSym),Smat(iOff1),nBas(iSym),Zero,CMO_Inv(iOff2), &
                nOrb(iSym))

    iOff1 = iOff1+nBas(iSym)**2
    iOff2 = iOff2+nOrb(iSym)*nBas(iSym)
  end do

  call mma_deallocate(Smat)
end if

call SetUp_CASPT2_Tra(nSym,nBas,nOrb,nIsh,nAsh,nFro,nDel,CMO,nDens,LuTri1,LuTri2,LuHlf2,LuHlf3)
iType = 3  ! Means that TraCtl is called by MCLR

if ((.not. newCho) .and. (StepType /= 'RUN2')) call TraCtl_Drv(iType,.true.,1)

! fetch some data from existing file LuTri1
! (from a previous MCLR run)
! and make it available to the module intgrl
! (LuTRI1=LuMOTRA)
if (TwoStep .and. (StepType == 'RUN2')) call put_temp_data_on_intgrl(LuMOTRA,nSym,nOrb)

! Init Cholesky informations
if (NewCho) then
  BufFrac = 0.3_wp
  call Cho_X_Init(irc,BufFrac)
  iSeed = 10
  do i=1,nsym
    LuAChoVec(i) = IsFreeUnit(iSeed)
    iseed = LuAChoVec(i)+1
    write(Fname,'(A4,I1)') 'CHTA',i
    call DANAME_MF_WA(LuAChoVec(i),Fname)
  end do
  do i=1,nsym
    LuIChoVec(i) = IsFreeUnit(iSeed)
    iSeed = LuIChoVec(i)+1
    write(Fname,'(A4,I1)') 'CHTI',i
    call DANAME_MF_WA(LuIChoVec(i),Fname)
  end do
  LuChoInt(1) = IsFreeUnit(iSeed)
  write(Fname,'(A4)') 'CHIN'
  call DANAME_MF_WA(LuChoInt(1),Fname)
  LuChoInt(2) = IsFreeUnit(iSeed)
  write(Fname,'(A4)') 'CHTU'
  call DANAME_MF_WA(LuChoInt(2),Fname)
end if

call DaClos(LuTri2)
call DaClos(LuHlf2)
call DaClos(LuHlf3)

call FckMat()
call StPert()

! With Cholesky there is no other choice than computing some
! integrals used for the preconditioner

if (NewCho) call cho_prec_mclr(CMO,nIsh,nASh,LuAChoVec,LuChoInt)

!----------------------------------------------------------------------*
!     exit                                                             *
!----------------------------------------------------------------------*

end subroutine Start_MCLR
