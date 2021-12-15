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
*               2019, Liviu Ungur                                      *
************************************************************************
      SubRoutine IvoGen_rasscf(nSym,nBas,nFro,nIsh,nAsh,nCMO,nEOrb,
     &                         CMO,EOrb)
************************************************************************
*                                                                      *
*     purpose: Generate improved virtual orbitals by diagonalization   *
*              of the one-electron hamiltonian in the subspace spanned *
*              by the virtual orbitals                                 *
*                                                                      *
*     input:                                                           *
*       CMO     : molecular orbital coefficients of length nCMO=NTOT2  *
*       EOrb    : on input contain orbital energies as obtained in the *
*                 NATORB_RASSCF, energies  of the virtual orbitals     *
*                 are non-zero, length nEOrb=NTOT                      *
*                                                                      *
*     output:                                                          *
*       CMO     : molecular orbital coefficients with virtual orbitals *
*                 modified                                             *
*       EOrb    : orbital energies (set to zero for virtual orbitals)  *
*                                                                      *
*     called from: OutCtl, only when IVO is used in the input          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Original MOLCAS/scf/ivogen.f written by:                         *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*                                                                      *
*     written by:                                                      *
*     L. Ungur  using the MOLCAS/scf/ivogen.f as start                 *
*     National University of Singapore,  Feb 2019                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit None
#include "real.fh"
#include "stdalloc.fh"
#include "warnings.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='IVOGEN_RASSCF   ')
*
      Integer  nCMO, nEOrb, nSym
      Integer  nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym)
      Real*8   CMO(nCMO),EOrb(nEOrb)
*
      Real*8, Dimension(:), Allocatable:: FckS, FckH, FckT, OneHam,
     &                                    Scratch
*
      ! clocal variables:
      Integer  iSym, ij, nOrbi, iCMO, iEOr, iDum, nFound, iErr
      Integer  nOcc(nSym)
      Integer  MaxBas, MaxBOO, MaxOrO, nBT
      Integer  iRc, iOpt, iComp, iSyLbl
      Real*8   Dummy
#ifdef _DEBUGPRINT_
      Integer  iOff, iBas
#endif

      nBT=0
      MaxBas=0
      MaxOrO=0
      MaxBOO=0
      nOcc(1:nSym)=0
      Do iSym=1, nSym
         nOcc(iSym) =  nFro(iSym) + nIsh(iSym) + nAsh(iSym)
         nBT    = nBT  + nBas(iSym)*(nBas(iSym) + 1)/2
         MaxBas = Max(MaxBas,nBas(iSym))
         MaxOrO = Max(MaxOrO,nBas(iSym) - nOcc(iSym))
         MaxBOO = Max(MaxBOO,nBas(iSym)*(nBas(iSym)-nOcc(iSym)))
      End Do
#ifdef _DEBUGPRINT_
      Do iSym=1, nSym
         Call RecPrt('IvoGen: CMO(in)',' ',CMO(iCMO),nBas(iSym),
     &                                               nBas(iSym) )
      End Do
#endif
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*---- Allocate memory for the core Hamiltonian
      Call mma_allocate(OneHam, nBT ,Label='OneHam')
      Call dcopy_(nBT,[Zero],0,OneHam,1)
*---- Load bare nuclei Hamiltonian

      iRc    = -1
      iOpt   =  6
      iComp  =  1
      iSyLbl =  1
      Call RdOne(iRc,iOpt,'OneHam  ',iComp,OneHam,iSyLbl)
      If ( iRc.ne.0 ) Then
        Write(LF,*)' RASSCF tried to construct compact virtual orbitals'
        Write(LF,*)' by diagonalization of core Hamiltonian, but ran   '
        Write(LF,*)' into a severe error: Failed to read the           '
        Write(LF,*)' Hamiltonian from the ONEINT file. Something may be'
        Write(LF,*)' wrong with the file.'
        Call Quit(_RC_IO_ERROR_READ_)
      End If
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*) ' OneHam in AO basis in RASSCF'
      Write(LF,*) ' ---------------------'
      Write(LF,*)
      iOff=0
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        Call TriPrt('IvoGen: OneHam:','(5G17.11)',OneHam(1+iOff),iBas)
        iOff = iOff + (iBas*iBas+iBas)/2
      End Do
#endif
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
         nOrbi = nBas(iSym) - nOcc(iSym)
*
*------- If nOrbi.eq.0 - no virtual orbitals; iCMO and iEOr must be
*        updated anyway (occupied orbitals may exist)
         iCMO = iCMO + nBas(iSym)*nOcc(iSym)
         iEOr = iEOr + nOcc(iSym)
*
         If (nOrbi.gt.0) Then
*
*---------- Transform OneHam to space spanned by virtual orbitals
            Call Square(OneHam(ij),FckS,1,nBas(iSym),nBas(iSym))
            ! multiply FckH = OneHam x CMO
            Call DGEMM_('N','N',
     &                  nBas(iSym),nOrbi,nBas(iSym),
     &                  1.0d0,FckS,nBas(iSym),
     &                        CMO(iCMO),nBas(iSym),
     &                  0.0d0,FckH,nBas(iSym))
            ! multiply FckT =  CMO x FckH
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
      Call mma_deallocate(OneHam)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
