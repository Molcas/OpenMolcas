************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Martin Schuetz                                   *
*               2017, Roland Lindh                                     *
************************************************************************
      SubRoutine RotMOs(Delta,nDelta,CMO,nCMO,nD,Ovrlp,mBT)
************************************************************************
*                                                                      *
*     purpose: rotates MOs according to last displacement vector       *
*              delta after QNR step and DIIS extrapolation.            *
*              only called during second order update (QNR) opt.       *
*              delta is taken as the last entry on LLDelt              *
*                                                                      *
*     input:                                                           *
*       Delta   : displacement vectors used to construct unitary       *
*                 rotation matrix via expKap                           *
*       CMO     : orthonormal vectors from previous iteration of       *
*                 length nCMO                                          *
*                                                                      *
*     output:                                                          *
*       CMO     : orthonormal vectors, rotated by U=exp(delta)         *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*     calls to: ExpKap                                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1994                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: VVUHF                                                   *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
#include "file.fh"
#include "llists.fh"
*
      Integer nDelta,nCMO
      Real*8 CMO(nCMO,nD),Delta(nDelta,nD),Ovrlp(mBT)
      Real*8 Cpu1,Tim1,Tim2,Tim3
*
*---- Define local variables
      Integer iSym,iSyBlpt,nOF,nVrt,nOccmF,iCMOpt
      Real*8, Dimension(:), Allocatable:: RoM, Scratch
*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*define _DEBUGPRINT_
*
      Call mma_allocate(RoM,nOFS,Label='RoM')
      nSize=0
      Do iSym=1,nSym
         nOF=nOrb(iSym)-nFro(iSym)
         nOFnBa=nOF*nBas(iSym)
         nSize=Max(nSize,nOFnBa)
      End Do
      Call mma_allocate(Scratch,nSize,Label='Scratch')
*
      Do iD = 1, nD
*        compute rotation matrix via expkap
         Call ExpKap(Delta(1,iD),RoM,nOcc(1,iD))
         iSyBlpt=1
         iCMOpt=1
*
*        loop over all symmetry blocks
*
         Do iSym=1,nSym
*
            nOF=nOrb(iSym)-nFro(iSym)
            nVrt=nOrb(iSym)-nOcc(iSym,iD)
            nOccmF=nOcc(iSym,iD)-nFro(iSym)
            nOFnBa=nOF*nBas(iSym)
            iCMOpt=iCMOpt+nBas(iSym)*nFro(iSym)
*
            If ((nVrt.gt.0).AND.(nOccmF.gt.0)) Then
*              skip, if no orbitals within this irrep
               call dcopy_(nOFnBa,CMO(iCMOpt,iD),1,Scratch,1)
               Call DGEMM_('N','N',
     &                     nBas(iSym),nOF,nOF,
     &                     1.0d0,Scratch,nBas(iSym),
     &                           RoM(iSyBlpt),nOF,
     &                     0.0d0,CMO(iCMOpt,iD),nBas(iSym))
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
               Call NrmClc(Scratch,nBas(iSym)*nOrb(iSym),'RotMOs',
     &                     'Old CMOs')
               Call NrmClc(CMo(iCMOpt,iD),nBas(iSym)*nOrb(iSym),
     &                     'RotMOs','New CMOs')
               Call RecPrt('RoM',' ',RoM(iSyBlpt),nOF,nOF)
*              Call RecPrt('RotMOs: Old CMOs',' ',Scratch,
*    &                     nBas(iSym),nOrb(iSym))
*              Call RecPrt('RotMOs: New CMOs',' ',CMO(iCMOpt,iD),
*    &                     nBas(iSym),nOrb(iSym))
#endif
               iSyBlpt=iSyBlpt+nOF*nOF
            End If
            iCMOpt=iCMOpt+nOF*nBas(iSym)
         End Do
*
*----    Check orthogonality
*
         Call ChkOrt(CMO(1,iD),nBO,Ovrlp,mBT,Whatever)
*
      End Do ! iD
*
      Call mma_deallocate(Scratch)
      Call mma_deallocate(RoM)
*
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call NrmClc(Delta,nDelta*nD,'RotMos','Delta')
      Call NrmClc(CMO,nCMO*nD,'RotMos','CMO')
      Call RecPrt('RotMOs: Delta',' ',Delta,nDelta,nD)
      Call RecPrt('RotMOs: CMO',' ',CMO,nCMO,nD)
#endif
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 9) = TimFld( 9) + (Cpu2 - Cpu1)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
