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
! Copyright (C) 1994, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RotMOs(Delta,nDelta)
!***********************************************************************
!                                                                      *
!     purpose: rotates MOs according to last displacement vector       *
!              delta after QNR step and DIIS extrapolation.            *
!              only called during second order update (QNR) opt.       *
!                                                                      *
!     input:                                                           *
!       Delta   : displacement vectors used to construct unitary       *
!                 rotation matrix via expKap                           *
!       CMO_ref : orthonormal vectors from reference iteration of      *
!                 length nCMO                                          *
!                                                                      *
!     output:                                                          *
!       CMO     : orthonormal vectors, rotated by U=exp(delta)         *
!                                                                      *
!                                                                      *
!***********************************************************************

use InfSCF, only: CMO, CMO_ref, kOV, nBas, nD, nFro, nOcc, nOFS, nOrb, nSym, TimFld
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDelta
real(kind=wp), intent(in) :: Delta(nDelta)
integer(kind=iwp) :: iCMOpt, iD, iEnd, iSt, iSyBlpt, iSym, nOccmF(8), nOF, nOfNBA, nSize, nVrt
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nCMO
#endif
real(kind=wp) :: Cpu1, CPU2, Tim1, Tim2, Tim3, WhatEver
real(kind=wp), allocatable :: RoM(:), Scratch(:)

call Timing(Cpu1,Tim1,Tim2,Tim3)

call mma_allocate(RoM,nOFS,Label='RoM')
nSize = 0
do iSym=1,nSym
  nOF = nOrb(iSym)-nFro(iSym)
  nOFnBa = nOF*nBas(iSym)
  nSize = max(nSize,nOFnBa)
end do
call mma_allocate(Scratch,nSize,Label='Scratch')

iEnd = 0
do iD=1,nD
  if (kOV(iD) < 1) cycle
  iSt = iEnd+1
  iEnd = iEnd+kOV(iD)
  ! compute rotation matrix via expkap
  nOccmF(:) = nOcc(:,iD)-nFro(:)
  call ExpKap(Delta(iSt:iEnd),kOV(id),RoM,nOccmF)
  iSyBlpt = 1
  iCMOpt = 1

  ! loop over all symmetry blocks

  do iSym=1,nSym

    nOF = nOrb(iSym)-nFro(iSym)
    nVrt = nOrb(iSym)-nOcc(iSym,iD)
    nOFnBa = nOF*nBas(iSym)
    iCMOpt = iCMOpt+nBas(iSym)*nFro(iSym)

    ! skip, if no orbitals within this irrep
    if ((nVrt > 0) .and. (nOccmF(iSym) > 0)) then
      Scratch(1:nOFnBA) = CMO_ref(iCMOpt:iCMOpt+nOFnBa-1,iD)
      call DGEMM_('N','N',nBas(iSym),nOF,nOF, &
                  One,Scratch,nBas(iSym), &
                  RoM(iSyBlpt),nOF, &
                  Zero,CMO(iCMOpt,iD),nBas(iSym))
#     ifdef _DEBUGPRINT_
      call NrmClc(Scratch,nBas(iSym)*nOrb(iSym),'RotMOs','Old CMOs')
      call NrmClc(CMO(iCMOpt,iD),nBas(iSym)*nOrb(iSym),'RotMOs','New CMOs')
      call RecPrt('RoM',' ',RoM(iSyBlpt),nOF,nOF)
      !call RecPrt('RotMOs: Old CMOs',' ',Scratch,nBas(iSym),nOrb(iSym))
      !call RecPrt('RotMOs: New CMOs',' ',CMO(iCMOpt,iD),nBas(iSym),nOrb(iSym))
#     endif
      iSyBlpt = iSyBlpt+nOF*nOF
    end if
    iCMOpt = iCMOpt+nOF*nBas(iSym)
  end do

  ! Check orthogonality

  call ChkOrt(iD,Whatever)

end do ! iD

call mma_deallocate(Scratch)
call mma_deallocate(RoM)

#ifdef _DEBUGPRINT_
nCMO = size(CMO,1)
call NrmClc(Delta,nDelta,'RotMos','Delta')
call NrmClc(CMO,nCMO*nD,'RotMos','CMO')
call RecPrt('RotMOs: Delta',' ',Delta,1,nDelta)
call RecPrt('RotMOs: CMO',' ',CMO,nCMO,nD)
#endif
call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(9) = TimFld(9)+(Cpu2-Cpu1)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine RotMOs
