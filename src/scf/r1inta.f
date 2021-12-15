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
      Subroutine R1IntA
************************************************************************
*                                                                      *
*     purpose: Read one-electron hamiltonian and overlap matrix        *
*                                                                      *
*     called from: ReadIn                                              *
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
      use SCF_Arrays
      Implicit Real*8 (a-h,o-z)
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
#include "real.fh"
#ifdef _FDE_
#include "embpotdata.fh"
#endif
*
*---- Define local variables
      Character*8 Label
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _FDE_
      ! Embedding
      iDummyEmb=0
      Call Get_iScalar('embpot', iDummyEmb)
      if (iDummyEmb.eq.1) embPot=.true.
      if (embPot) then
         call mma_allocate(Emb,nBT,Label='Emb')
         ipEmb=ip_of_Work(Emb(1))
         Call EmbPotRdRun
      end if
#endif
*---- Allocate memory for one-electron integrals
      Call mma_allocate(OneHam,nBT,Label='OneHam')
      Call mma_allocate(Ovrlp,nBT+4,Label='Ovrlp')
      Call FZero(OneHam,nBT)
      Call FZero(Ovrlp,nBT+4)
*
*---- Read core Hamiltonian
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      Label='OneHam  '
      Call RdOne(iRc,iOpt,Label,iComp,OneHam,iSyLbl)
      If (iRc.ne.0) Then
         Write (6,*) 'R1Inta: Error readin ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
#ifdef _FDE_
      ! Embedding
      if (embPot) then
       if (embPotInBasis) then
        ! If the potential is given in basis set representation it
        ! has not been calculated with a OneEl call and is just read
        ! from file here.
        iunit = isFreeUnit(1)
        call molcas_open(iunit, embPotPath)
        do iEmb=1, nBT
         read(iunit,*) Emb(iEmb)
        end do
        close(iunit)
       else
        ! Read one-electron integrals due to embedding potential
        iRc=-1
        Label='embpot  '
        Call RdOne(iRc,iOpt,Label,iComp,Emb,iSyLbl)
        If (iRc.ne.0) Then
           Write (6,*) 'R1Inta: Error readin ONEINT'
           Write (6,'(A,A)') 'Label=',Label
           Call Abend()
        End If
       end if
      end if
#endif
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      ist1Hm=1
      Write (6,*)
      Write (6,*) ' One electron Hamiltonian at start'
      Write (6,*) ' ---------------------------------'
      Do iSym=1,nSym
         Write (6,*) ' symmetry block',iSym
         Call TriPrt(' ',' ',OneHam(ist1Hm),nBas(iSym))
         ist1Hm=ist1Hm+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
#endif
*
*---- Read overlap integrals and total effective nuclear charge
      iRc=-1
      iOpt=2
      iComp=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,iComp,Ovrlp,iSyLbl)
      If (iRc.ne.0) Then
         Write (6,*) 'R1Inta: Error readin ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
      Tot_Nuc_Charge=Ovrlp(nBT+4)
#ifdef _DEBUGPRINT_
      istOvl=1
      Write (6,*)
      Write (6,*) ' Overlap matrix at start'
      Write (6,*) ' -----------------------'
      Do iSym=1,nSym
         Write (6,*) ' symmetry block',iSym
         Call TriPrt(' ',' ',Ovrlp(istOvl),nBas(iSym))
         istOvl=istOvl+nBas(iSym)*(nBAs(iSym)+1)/2
      End Do
#endif
*
*---- Perform Lowdin orthonormalization
*
      nD = iUHF+1
      Call mma_allocate(Lowdin,nBB,nD,Label='Lowdin')
      Call FZero(Lowdin,nBB*nD)
      ioff=1
      Do iSym=1,nSym
         Call dCopy_(nBas(iSym),[One],0,Lowdin(ioff,1),nBas(iSym)+1)
         ioff=ioff+nBas(iSym)*nBas(iSym)
      End Do
      Call xxLowdin(Ovrlp,Lowdin(1,1),nBas,nSym)
      If (nD.eq.2) Call DCopy_(nBB,Lowdin(1,1),1,Lowdin(1,2),1)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
