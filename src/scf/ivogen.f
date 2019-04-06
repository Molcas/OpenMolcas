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
      SubRoutine IvoGen(OneHam,nOne,CMO,nCMO,EOrb,nEOrb,mynOcc)
************************************************************************
*                                                                      *
*     purpose: Generate improved virtual orbitals by diagonalization   *
*              of the one-electron hamiltonian in the subspace spanned *
*              by the virtual orbitals                                 *
*                                                                      *
*     input:                                                           *
*       OneHam  : one-electron hamiltonian of length nOne              *
*       CMO     : molecular orbital coefficients of length nCMO        *
*                                                                      *
*     output:                                                          *
*       CMO     : molecular orbital coefficients with virtual orbitals *
*                 modified                                             *
*       EOrb    : orbital energies (set to zero for virtual orbitals)  *
*                                                                      *
*     called from: Final                                               *
*                                                                      *
*     calls to: PickUp, Sort                                           *
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
      Real*8 OneHam(nOne),CMO(nCMO),EOrb(nEOrb)
      Integer mynOcc(*)
*
      Real*8, Dimension(:), Allocatable:: FckS, FckH, FckT,
     &                                    Scratch
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*---- Allocate memory for squared modified Fock matrix
      Call mma_allocate(FckS,MaxBas**2,Label='FckS')
*
*---- Allocate memory for half transformed Fock matrix
      Call mma_allocate(FckH,MaxBOO,Label='FckH')
*
*---- Allocate memory for transformed Fock matrix
      Call mma_allocate(FckT,MaxOrO*(MaxOrO + 1)/2,Label='FckT')
*
      ij   = 1
      iCMO = 1
      iEOr = 1
      Do iSym = 1, nSym
         nOrbi = nOrb(iSym) - mynOcc(iSym)
*
*------- If nOrbi.eq.0 - no virtual orbitals; iCMO and iEOr must be
*        updated anyway (occupied orbitals may exist)
         iCMO = iCMO + nBas(iSym)*mynOcc(iSym)
         iEOr = iEOr + mynOcc(iSym)
*
         If (nOrbi.gt.0) Then
*
*---------- Transform OneHam to space spanned by virtual orbitals
            Call Square(OneHam(ij),FckS,1,nBas(iSym),nBas(iSym))
            Call DGEMM_('N','N',
     &                  nBas(iSym),nOrbi,nBas(iSym),
     &                  1.0d0,FckS,nBas(iSym),
     &                        CMO(iCMO),nBas(iSym),
     &                  0.0d0,FckH,nBas(iSym))
            Call MxMt(CMO(iCMO),   nBas(iSym),1,
     &                FckH,1,nBas(iSym),
     &                FckT,
     &                nOrbi,nBas(iSym))
*
*---------- Diagonalize OneHam within virtual space and form orbital energies
            Call mma_allocate(Scratch,nOrbi**2,Label='Scratch')
            Dummy=0.0D0
            iDum=0
            Call Diag_Driver('V','A','L',nOrbi,FckT,
     &                       Scratch,nOrbi,Dummy,Dummy,iDum,iDum,
     &                       EOrb(iEOr),CMO(iCMO),nBas(iSym),0,-1,'J',
     &                       nFound,iErr)
            Call mma_deallocate(Scratch)
*
*---------- Orbital energies are now meaningless; set them to zero
            call dcopy_(nOrbi,[Zero],0,EOrb(iEOr),1)
*
         End If
*
*------- Update pointers
         iCMO = iCMO + nOrbi*nBas(iSym)
         iEOr = iEOr + nOrbi
         ij   = ij   + nBas(iSym)*(nBas(iSym) + 1)/2
*
      End Do
*
*---- Deallocate memory
      Call mma_deallocate(FckS)
      Call mma_deallocate(FckH)
      Call mma_deallocate(FckT)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
