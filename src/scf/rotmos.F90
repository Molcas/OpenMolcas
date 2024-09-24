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
      SubRoutine RotMOs(Delta,nDelta)
!***********************************************************************
!                                                                      *
!     purpose: rotates MOs according to last displacement vector       *
!              delta after QNR step and DIIS extrapolation.            *
!              only called during second order update (QNR) opt.       *
!                                                                      *
!     input:                                                           *
!       Delta   : displacement vectors used to construct unitary       *
!                 rotation matrix via expKap                           *
!       CMO     : orthonormal vectors from previous iteration of       *
!                 length nCMO                                          *
!                                                                      *
!     output:                                                          *
!       CMO     : orthonormal vectors, rotated by U=exp(delta)         *
!                                                                      *
!                                                                      *
!***********************************************************************
      use InfSCF, only: nSym, kOV, nBas, nFro, nOcc, NoFS, nOrb, TimFld, nD
      use stdalloc, only: mma_allocate, mma_deallocate
      use SCF_Arrays, only: CMO
      use Constants, only: Zero, One
      Implicit None
!
      Integer nDelta
      Real*8  Delta(nDelta)
!
!---- Define local variables
      Integer iSym,iSyBlpt,nOF,nVrt,nOccmF,iCMOpt, nSize, nOfNBA, iEnd, iD, iSt
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Integer nCMO
#endif
      Real*8, Dimension(:), Allocatable:: RoM, Scratch
      Real*8 Cpu1,CPU2,Tim1,Tim2,Tim3, WhatEver
!
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
!
      Call mma_allocate(RoM,nOFS,Label='RoM')
      nSize=0
      Do iSym=1,nSym
         nOF=nOrb(iSym)-nFro(iSym)
         nOFnBa=nOF*nBas(iSym)
         nSize=Max(nSize,nOFnBa)
      End Do
      Call mma_allocate(Scratch,nSize,Label='Scratch')
!
      iEnd = 0
      Do iD = 1, nD
         iSt = iEnd + 1
         iEnd = iEnd + kOV(iD)
!        compute rotation matrix via expkap
         Call ExpKap(Delta(iSt:iEnd),kOV(id),RoM,nOcc(1,iD))
         iSyBlpt=1
         iCMOpt=1
!
!        loop over all symmetry blocks
!
         Do iSym=1,nSym
!
            nOF=nOrb(iSym)-nFro(iSym)
            nVrt=nOrb(iSym)-nOcc(iSym,iD)
            nOccmF=nOcc(iSym,iD)-nFro(iSym)
            nOFnBa=nOF*nBas(iSym)
            iCMOpt=iCMOpt+nBas(iSym)*nFro(iSym)
!
            If ((nVrt.gt.0).AND.(nOccmF.gt.0)) Then
!              skip, if no orbitals within this irrep
               call dcopy_(nOFnBa,CMO(iCMOpt,iD),1,Scratch,1)
               Call DGEMM_('N','N',nBas(iSym),nOF,nOF,       &
                           One,Scratch,nBas(iSym),           &
                                 RoM(iSyBlpt),nOF,           &
                           Zero,CMO(iCMOpt,iD),nBas(iSym))
#ifdef _DEBUGPRINT_
               Call NrmClc(Scratch,nBas(iSym)*nOrb(iSym),'RotMOs','Old CMOs')
               Call NrmClc(CMo(iCMOpt,iD),nBas(iSym)*nOrb(iSym),'RotMOs','New CMOs')
               Call RecPrt('RoM',' ',RoM(iSyBlpt),nOF,nOF)
!              Call RecPrt('RotMOs: Old CMOs',' ',Scratch,nBas(iSym),nOrb(iSym))
!              Call RecPrt('RotMOs: New CMOs',' ',CMO(iCMOpt,iD),nBas(iSym),nOrb(iSym))
#endif
               iSyBlpt=iSyBlpt+nOF*nOF
            End If
            iCMOpt=iCMOpt+nOF*nBas(iSym)
         End Do
!
!----    Check orthogonality
!
         Call ChkOrt(iD,Whatever)
!
      End Do ! iD
!
      Call mma_deallocate(Scratch)
      Call mma_deallocate(RoM)
!
#ifdef _DEBUGPRINT_
      nCMO = Size(CMO,1)
      Call NrmClc(Delta,nDelta,'RotMos','Delta')
      Call NrmClc(CMO,nCMO*nD,'RotMos','CMO')
      Call RecPrt('RotMOs: Delta',' ',Delta,1,nDelta)
      Call RecPrt('RotMOs: CMO',' ',CMO,nCMO,nD)
#endif
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 9) = TimFld( 9) + (Cpu2 - Cpu1)
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Return
      End SubRoutine RotMOs
