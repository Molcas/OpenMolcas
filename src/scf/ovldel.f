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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      SubRoutine OvlDel(Ovlp,nOvlp,TrMat,nTrMat)
************************************************************************
*                                                                      *
*     purpose: Remove near linear dependencies from basis set          *
*                                                                      *
*     input:                                                           *
*       Ovlp    : overlap in AO basis of length nOvlp                  *
*       TrMat   : unit matrix or matrix transforming from AO's to      *
*                 cartesian functions (dependently on input) of        *
*                 length nTrMat                                        *
*                                                                      *
*     output:                                                          *
*       TrMat   : input matrix modified such, that near linear de-     *
*                 pendencies (if exist) are removed                    *
*                                                                      *
*     called from: TrGen                                               *
*                                                                      *
*     calls to: PickUp                                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
      Real*8 Ovlp(nOvlp),TrMat(nTrMat)
*
      Real*8, Dimension(:), Allocatable:: OvlT, OvlH, OvlS, EVec, EVal,
     &                                    NewB, Scratch


*---- Allocate memory for transformed overlap matrix
      Call mma_allocate(OvlT,MaxOrF*(MaxOrF+1)/2,Label='OvlT')
*
*---- Allocate memory for half-transformed overlap matrix
      Call mma_allocate(OvlH,MaxBOF,Label='OvlH')
*
*---- Allocate memory for squared overlap matrix
      Call mma_allocate(OvlS,MaxBas**2,Label='OvlS')
*
*---- Allocate memory for eigenvectors of overlap
      Call mma_allocate(EVec,MaxOrF**2,Label='EVec')
*
*---- Allocate memory for eigenvalues of overlap
      Call mma_allocate(EVal,MaxOrF,Label='EVal')
*
*---- Allocate memory for 'basis' that diagonalizes overlap
      Call mma_allocate(NewB,MaxBOF,Label='NewB')
*
      ij   = 1
      iOld = 1
      iNew = 1
      Do iSym = 1, nSym
*
         iiBT = nBas(iSym)*(nBas(iSym) + 1)/2
         nOF  = nOrb(iSym) - nFro(iSym)
*
*------- Copy frozen vectors to the right position
         If(nFro(iSym)*nBas(iSym).gt.0) Then
*           Write(6,'(a,i2)') 'Copying symmetry block',iSym
*           Write(6,'(i8,a,i8)') iOld,' ->',iNew
*           Write(6,'(a,i8)') 'nTrMat =',nTrMat
            call dcopy_(nFro(iSym)*nBas(iSym),
     &                 TrMat(iOld),1,TrMat(iNew),1)
*        Else
*           Write(6,'(a,i2)') 'No copying of symmetry block',iSym
         End If
*...+....1....+....2....+....3....+....4....+....5....+....6....+....7.>..+....8
*
         iOld = iOld + nFro(iSym)*nBas(iSym)
         iNew = iNew + nFro(iSym)*nBas(iSym)
*
         If (nOF.gt.0) Then
*
*---------- Square overlap and transform to basis given by TrMat
            Call Square(Ovlp(ij),OvlS,1,nBas(iSym),nBas(iSym))
            Call DGEMM_('N','N',
     &                  nBas(iSym),nOF,nBas(iSym),
     &                  1.0d0,OvlS,nBas(iSym),
     &                        TrMat(iOld),nBas(iSym),
     &                  0.0d0,OvlH,nBas(iSym))
            Call MxMt(TrMat(iOld),nBas(iSym),1,
     &                OvlH,1,nBas(iSym),
     &                OvlT,
     &                nOF,nBas(iSym))
*
*---------- Diagonalize overlap and form eigenvalues vector
            Call mma_allocate(Scratch,nOF**2,Label='Scrtach')
            Dummy=0.0D0
            iDum=0
            Call Diag_Driver('V','A','L',nOF,OvlT,
     &                       Scratch,nOF,Dummy,Dummy,iDum,iDum,
     &                       EVal,EVec,nOF,1,0,'J',
     &                       nFound,iErr)
            Call mma_deallocate(Scratch)
C??         Call dCopy_(nOF*(nOF+1)/2,Zero,0,OvlT,1)
C??         iDiag=0
C??         Do i = 1, nOF
C??            OvlT(i+iDiag) = EVal(i)
C??            iDiag = iDiag + i
C??         End Do
*
*---------- Transform to basis that diagonalizes overlap
            Call DGEMM_('N','N',
     &                  nBas(iSym),nOF,nOF,
     &                  1.0d0,TrMat(iOld),nBas(iSym),
     &                        EVec,nOF,
     &                  0.0d0,NewB,nBas(iSym))
*
*---------- Remove linear dependencies
            nOrbi = nFro(iSym)
            ind   = 1
            Do iOrb = 1, nOF
               If (EVal(iOrb).gt.DelThr) Then
                  If(EVal(iOrb).lt.1.0d-5) MiniDn=.False.
                  call dcopy_(nBas(iSym),NewB(ind),1,TrMat(iNew),1)
                  iNew  = iNew  + nBas(iSym)
                  nOrbi = nOrbi + 1
               End If
               ind = ind + nBas(iSym)
            End Do
*
            iOld       = iOld + nOF*nBas(iSym)
            nDel(iSym) = nOrb(iSym) - nOrbi
            nOrb(iSym) = nOrbi
*
         End If
*
         ij = ij + iiBT
*
      End Do
      Call Put_iArray('nDel',nDel,nSym)
*
*---- Deallocate memory
      Call mma_deallocate(NewB)
      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
      Call mma_deallocate(OvlS)
      Call mma_deallocate(OvlH)
      Call mma_deallocate(OvlT)

      End subroutine OvlDel
