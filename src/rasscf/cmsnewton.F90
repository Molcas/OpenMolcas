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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 11, 2022, created this file.               *
!*****************************************************************

subroutine CMSNewton(R,GDorbit,GDstate,Dgorbit,Dgstate,nGD)

use Index_Functions, only: nTri_Elem
use CMS, only: CMSNotConverged, CMSThres, LargestQaaGrad, NCMSScale, NeedMoreStep, nPosHess
use rasscf_global, only: CMSThreshold, iCMSIterMax, iCMSIterMin, lRoots, NAC
use PrintLevel, only: USUAL
use output_ras, only: IPRLOC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGD
real(kind=wp), intent(inout) :: R(lRoots**2), GDstate(nGD), Dgstate(nGD)
real(kind=wp), intent(out) :: GDorbit(nGD), Dgorbit(nGD)
integer(kind=iwp) :: iPrLev, iStep, lRoots2, NAC2, nDDg, nScr, nSPair, nSPair2
real(kind=wp) :: Qnew, Qold
logical(kind=iwp) :: Saved
real(kind=wp), allocatable :: DDg(:), deltaR(:), DgCopy(:), EigVal(:), GDCopy(:), Grad(:), GScr(:), Hess(:), RCopy(:), ScrDiag(:), &
                              X(:), XScr(:)

IPRLEV = IPRLOC(6)

! preparation
lRoots2 = lRoots**2
NAC2 = NAC**2
nDDg = lRoots2**2
nSPair = nTri_Elem(lRoots-1)
nSPair2 = nSPair**2
CMSThres = CMSThreshold
call mma_allocate(DDg,nDDg)
call mma_allocate(X,nSPair)
call mma_allocate(XScr,nSPair)
call mma_allocate(GScr,nSPair)
call mma_allocate(Grad,nSPair)
call mma_allocate(Hess,nSPair2)
call mma_allocate(EigVal,nSPair)
call mma_allocate(DeltaR,lRoots2)
call mma_allocate(GDCopy,nGD)
call mma_allocate(DgCopy,nGD)
call mma_allocate(RCopy,lRoots2)
! Step 0
iStep = 0
Qold = Zero
! Note that the following six lines appear as a group
call RotGD(GDstate,R,nGD,lRoots,NAC2)
call RotGD(Dgstate,R,nGD,lRoots,NAC2)
call Trnsps(lRoots2,NAC2,Dgstate,Dgorbit)
call Trnsps(lRoots2,NAC2,GDstate,GDorbit)
call DGEMM_('T','N',lRoots2,lRoots2,NAC2,One,Dgorbit,NAC2,GDorbit,NAC2,Zero,DDg,lRoots2)
call CalcQaa(Qnew,DDg,lRoots,nDDg)
nPosHess = 0
LargestQaaGrad = Zero
Qold = Qnew
call PrintCMSIter(iStep,Qnew,Qold,R,lRoots)
call CalcGradCMS(Grad,DDg,lRoots,nSPair)
call CalcHessCMS(Hess,DDg,lRoots,nSPair)
call GetDiagScr(nScr,Hess,EigVal,nSPair)
call mma_allocate(ScrDiag,nScr)

! Starting iteration
do while (CMSNotConverged)
  iStep = iStep+1
  Qold = Qnew
  if (iStep > iCMSIterMax) then
    write(6,'(4X,A)') 'NOT CONVERGED AFTER MAX NUMBER OF CYCLES'
    exit
  end if

  RCopy(:) = R(:)
  GDCopy(:) = GDState(:)
  DgCopy(:) = DgState(:)

  call CalcNewX(X,Hess,Grad,nSPair,XScr,GScr,EigVal,ScrDiag,nScr)
  call UpDateRotMat(R,DeltaR,X,lRoots,nSPair)

  call RotGD(GDstate,DeltaR,nGD,lRoots,NAC2)
  call RotGD(Dgstate,DeltaR,nGD,lRoots,NAC2)
  call Trnsps(lRoots2,NAC2,Dgstate,Dgorbit)
  call Trnsps(lRoots2,NAC2,GDstate,GDorbit)
  call DGEMM_('T','N',lRoots2,lRoots2,NAC2,One,Dgorbit,NAC2,GDorbit,NAC2,Zero,DDg,lRoots2)
  call CalcQaa(Qnew,DDg,lRoots,nDDg)

  NCMSScale = 0
  Saved = .true.
  if ((Qold-Qnew) > CMSThreshold) then
    ! When Onew is less than Qold, scale the rotation matrix
    if (iStep > ICMSIterMin) call CMSScaleX(X,R,DeltaR,Qnew,Qold,RCopy,GDCopy,DgCopy,GDstate,GDOrbit,Dgstate,DgOrbit,DDg,nSPair, &
                                            lRoots2,nGD,NAC2,nDDg,Saved)
  end if
  if (IPRLEV >= USUAL) call PrintCMSIter(iStep,Qnew,Qold,R,lRoots)
  call PrintMat('ROT_VEC','CMS-PDFT temp',R,lroots,lroots,7,13,'T')

  if (.not. Saved) then
    CMSNotConverged = .true.
    !exit
  end if
  ! sanity check
  if (abs(Qnew-Qold) < CMSThreshold) then
    CMSNotConverged = .false.
    if (NeedMoreStep) CMSNotConverged = .true.
    if (iStep < iCMSIterMin) CMSNotConverged = .true.
    if (NCMSScale > 0) CMSNotConverged = .true.
  end if
  if (CMSNotConverged) then
    call CalcGradCMS(Grad,DDg,lRoots,nSPair)
    call CalcHessCMS(Hess,DDg,lRoots,nSPair)
  else
    if (IPRLEV >= USUAL) write(6,'(4X,A)') 'CONVERGENCE REACHED'
  end if
end do

call mma_deallocate(DDg)
call mma_deallocate(X)
call mma_deallocate(XScr)
call mma_deallocate(GScr)
call mma_deallocate(Grad)
call mma_deallocate(Hess)
call mma_deallocate(EigVal)
call mma_deallocate(DeltaR)
call mma_deallocate(ScrDiag)
call mma_deallocate(GDCopy)
call mma_deallocate(DgCopy)
call mma_deallocate(RCopy)

end subroutine CMSNewton
