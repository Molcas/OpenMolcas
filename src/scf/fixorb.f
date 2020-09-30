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
      SubRoutine FixOrb(Ovrlp,CMO,TrMat,nCMO)
************************************************************************
*                                                                      *
*     purpose: Diagonalize core hamiltonian to get starting orbitals.  *
*                                                                      *
*     input:                                                           *
*       OneHam  : one-electron hamiltonian of length nOH               *
*       TrMat   : matrix transforming to orthonotmal basis of          *
*                 length nCMO                                          *
*                                                                      *
*     output:                                                          *
*       CMO     : starting vectors of length nCMO                      *
*                                                                      *
*     called from: Start1                                              *
*                                                                      *
*     calls to: ModFck, PickUp, Sort                                   *
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
      Real*8 Ovrlp(nCMO),CMO(nCMO),TrMat(nCMO)
*
      Real*8, Dimension(:), Allocatable:: S, TT, TTS, CMO0
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate memory
*
      Call mma_allocate(S,MaxBas**2,Label='S')
      Call mma_allocate(TT,MaxBas**2,Label='TT')
      Call mma_allocate(TTS,MaxBas**2,Label='TTS')
      Call mma_allocate(CMO0,MaxBas**2,Label='CMO0')
*                                                                      *
************************************************************************
*                                                                      *
*     Observe that TrGen has been called and that the dimensions might
*     have been change. Note also that at this point work in the full
*     basis.
*
*---- Loop over symmetry blocks
*
      iS     = 1
      iCMO   = 1
      iTrM   = 1
      iOvrlp = 1
      Do iSym = 1, nSym
#ifdef _DEBUG_
         Call RecPrt('FixOrb: CMO(in)',' ',CMO(iCMO),
     &                       nBas(iSym),nBas(iSym))
         Call RecPrt('FixOrb: TrMat',' ',TrMat(iCMO),
     &                       nBas(iSym),nOrb(iSym))
#endif
*
         nOF  = nOrb(iSym) - nFro(iSym)
         nBF  = nBas(iSym) - nFro(iSym)
*
*------- Skip the frozen orbitals
*
         iCMO   = iCMO + nBas(iSym)*nFro(iSym)
         iTrM   = iTrM + nBas(iSym)*nFro(iSym)
*
         If (nBF.gt.0) Then
*
*---------- Project to the reduced basis
*
*---------- Create the project operator
            Call DGEMM_('N','T',
     &                  nBas(iSym),nBas(iSym),nOF,
     &                  1.0d0,TrMat(iTrM),nBas(iSym),
     &                        TrMat(iTrM),nBas(iSym),
     &                  0.0d0,TT,nBas(iSym))
            Call Square(Ovrlp(iS),S,1,nBas(iSym),nBas(iSym))
            Call DGEMM_('N','N',
     &                  nBas(iSym),nBas(iSym),nBas(iSym),
     &                  1.0d0,TT,nBas(iSym),
     &                        S,nBas(iSym),
     &                  0.0d0,TTS,nBas(iSym))
#ifdef _DEBUG_
            Call RecPrt('FixOrb: TT',' ',TT,nBas(iSym),nBas(iSym))
            Call RecPrt('FixOrb: TTS',' ',TTS,nBas(iSym),nBas(iSym))
#endif
*
*---------- Project
            Call DGEMM_('N','N',
     &                  nBas(iSym),nOF,nBas(iSym),
     &                  1.0d0,TTS,nBas(iSym),
     &                        CMO(iCMO),nBas(iSym),
     &                  0.0d0,CMO0,nBas(iSym))
            call dcopy_(nBas(iSym)*nOF,CMO0,1,CMO(iCMO),1)
#ifdef _DEBUG_
            Call RecPrt('FixOrb: CMO(out)',' ',CMO(iCMO),
     &                             nBas(iSym),nBas(iSym))
#endif
*
         End If
*
*------- Update pointers
         iCMO = iCMO + nBF*nBas(iSym)
         iTrM = iTrM + nOF*nBas(iSym)
         iS   = iS   + nBas(iSym)*(nBas(iSym)+1)/2
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Deallocate memory
*
      Call mma_deallocate(CMO0)
      Call mma_deallocate(TTS)
      Call mma_deallocate(TT)
      Call mma_deallocate(S)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
#endif
      Return
      End
