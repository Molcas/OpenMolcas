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
      SubRoutine DCore(OneHam,nOH,CMO,TrMat,nCMO,EOr,nEOr,mynOcc,Ovrlp)
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
      Real*8 OneHam(nOH),CMO(nCMO),TrMat(nCMO),EOr(nEOr), Ovrlp(nOH)
      Integer mynOcc(*)
*
      Real*8, Dimension(:), Allocatable:: OMod, OHSq, OHHl, OHTr, EiVe,
     &                                    Scratch
      Integer, Dimension(:), Allocatable:: Fermi
      Data iseed/13/
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*---- Allocate memory for modified one-electron hamoltonian
      Call mma_allocate(OMod,nBT,Label='OMod')
*
*---- Allocate memory for squared one-el. Hamiltonian
      Call mma_allocate(OHSq,MaxBas**2,Label='OHSq')
*
*---- Allocate memory for half-transformed one-el. Hamiltonian
      Call mma_allocate(OHHl,MaxBOF,Label='OHHl')
*
*---- Allocate memory for transformed one-electron Hamiltonian
      Call mma_allocate(OHTr,MaxOrF*(MaxOrF+1)/2,Label='OHTr')
*
*---- Allocate memory for eigenvectors
      Call mma_allocate(EiVe,MaxOrF**2,Label='EiVe')
*
*---- Allocate memory for info on electronic and muonic basis sets
      Call mma_allocate(Fermi,nEOr,Label='Fermi')
      Call Get_iArray('Fermion IDs',Fermi,nEOr)
*
*
*
*---- Modify one-electron hamiltonian
      call dcopy_(nBT,OneHam,1,OMod,1)
      If (nnFr.gt.0)
     &   Call ModFck(OMod,Ovrlp,nBT,TrMat,nBO,mynOcc)
*
*---- Diagonalize core in non-frozen molecular basis
      ij   = 1
      iCMO = 1
      iEOr = 1
      Do iSym = 1, nSym
*
         iiBT = nBas(iSym)*(nBas(iSym) + 1)/2
         nOF  = nOrb(iSym) - nFro(iSym)
*
*------- Copy frozen vectors to CMO array
         if(nFro(iSym)*nBas(iSym).gt.0) then
          call dcopy_(nFro(iSym)*nBas(iSym),TrMat(iCMO),1,CMO(iCMO),1)
         endif
*
         iCMO = iCMO + nBas(iSym)*nFro(iSym)
         iEOr = iEOr + nFro(iSym)
*
         If (nOF.gt.0) Then
*
*---------- Square one-el. Ham. and transform to orthonormal basis.
*           Call Square(OMod(ij),nBas(iSym),OHSq,nBas(iSym))
            Call Square(OMod(ij),OHSq,1,nBas(iSym),nBas(iSym))
            Call DGEMM_('N','N',
     &                  nBas(iSym),nOF,nBas(iSym),
     &                  1.0d0,OHSq,nBas(iSym),
     &                        TrMat(iCMO),nBas(iSym),
     &                  0.0d0,OHHl,nBas(iSym))
            Call MxMt(TrMat(iCMO), nBas(iSym),1,OHHl,1,nBas(iSym),
     &                OHTr,nOF,nBas(iSym))
*
*---------- Put a unit matrix into the eigenvector matrix
            call dcopy_(nOF*nOF,Zero,0,EiVe,      1)
            call dcopy_(nOF,    One, 0,EiVe,nOF + 1)
*
*---------- Add small random number to the one-electron Hamiltonian
*           Not done here anymore, done in routine scram called
*           by sorb!
*
            If (Scrmbl .and. .false.) Then
c               iseed = 13
               ind   = 1
               Do i = 1, nOF
                  Do j = 1, i - 1
                     OHTr(ind) =
     &               OHTr(ind) + 0.050d+00*Random_Molcas(iseed)
                     ind = ind + 1
                  End Do
                  ind = ind + 1
               End Do
            End If
*
*---------- Diagonalize and form orbital energies
            Call mma_allocate(Scratch,nOF**2,Label='Scratch')
            Dummy=0.0D0
            iDum=0
            Call Diag_Driver('V','A','L',nOF,OHTr,
     &                       Scratch,nOF,Dummy,Dummy,iDum,iDum,
     &                       EOr(iEOr),EiVe,nOF,1,-1,'J',nFound,
     &                       iErr)
            Call mma_deallocate(Scratch)
*
*---------- Transform to AO basis
            Call DGEMM_('N','N',
     &                  nBas(iSym),nOF,nOF,
     &                  1.0d0,TrMat(iCMO),nBas(iSym),
     &                        EiVe,nOF,
     &                  0.0d0,CMO(iCMO),nBas(iSym))
*
         End If
*
*------- Update pointers
         iCMO = iCMO + nOF*nBas(iSym)
         iEOr = iEOr + nOF
         ij   = ij   + iiBT
*
      End Do
*
*---- Deallocate memory
      Call mma_deallocate(Fermi)
      Call mma_deallocate(EiVe)
      Call mma_deallocate(OHTr)
      Call mma_deallocate(OHHl)
      Call mma_deallocate(OHSq)
      Call mma_deallocate(OMod)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
