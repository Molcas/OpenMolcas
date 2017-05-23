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
* Copyright (C) 2006, Per-Olof Widmark                                 *
************************************************************************
      Subroutine NatoUHF(DensA,DensB,FockA,FockB,nBT,CMO,nBB,Ovl,
     &                   Nato,Eta,Eps,nnB,
     &                   nSym,nBas,nOrb)
************************************************************************
*                                                                      *
* This routine computes the natural orbitals for a UHF wavefunction.   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
* Written: September 2006                                              *
*                                                                      *
************************************************************************
      Implicit None
#include "stdalloc.fh"
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Integer nBT, nBB, nnB
      Real*8  DensA(nBT)
      Real*8  DensB(nBT)
      Real*8  FockA(nBT)
      Real*8  FockB(nBT)
      Real*8  CMO(nBB)
      Real*8  Ovl(nBT)
      Real*8  Nato(nBB)
      Real*8  Eta(nnB)
      Real*8  Eps(nnB)
      Integer nSym
      Integer nBas(nSym)
      Integer nOrb(nSym)
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer MaxTri
      Integer MaxSqr
      Integer nTri
      Integer nSqr
      Integer nCMO
      Integer nEta
      Integer iOffTri
      Integer iOffSqr
      Integer iOffCMO
      Integer iOffEta
      Integer iSym
      Integer iOrb
      Integer indx
      Integer i
      Real*8  tmp
      Real*8, Dimension(:), Allocatable:: Dens, Fock, SMat, Aux1, Aux2,
     &                                    Aux3
*----------------------------------------------------------------------*
* Setup                                                                *
*----------------------------------------------------------------------*
      MaxTri=0
      MaxSqr=0
      nTri=0
      nSqr=0
      nCMO=0
      nEta=0
      Do iSym=1,nSym
         MaxTri=Max(MaxTri,nBas(iSym)*(nBas(iSym)+1)/2)
         MaxSqr=Max(MaxSqr,nBas(iSym)*nBas(iSym))
         nTri=nTri+nBas(iSym)*(nBas(iSym)+1)/2
         nSqr=nSqr+nBas(iSym)*nBas(iSym)
         nCMO=nCMO+nBas(iSym)*nOrb(iSym)
         nEta=nEta+nOrb(iSym)
      End Do
*----------------------------------------------------------------------*
* Allocate arrays                                                      *
*----------------------------------------------------------------------*
      Call mma_allocate(Dens,nTri,Label='Dens')
      Call mma_allocate(Fock,nTri,Label='Fock')
      Call mma_allocate(SMat,MaxSqr,Label='SMat')
      Call mma_allocate(Aux1,MaxSqr,Label='Aux1')
      Call mma_allocate(Aux2,MaxSqr,Label='Aux2')
      Call mma_allocate(Aux3,MaxSqr,Label='Aux3')
*----------------------------------------------------------------------*
* Add up the densities                                                 *
*----------------------------------------------------------------------*
      Do i=1,nTri
        Dens(i)=DensA(i)+DensB(i)
      End Do
*----------------------------------------------------------------------*
* Copy orbitals                                                        *
*----------------------------------------------------------------------*
      Do i=1,nCMO
         Nato(i)=CMO(i)
      End Do
*----------------------------------------------------------------------*
* Average Fock matrix.                                                 *
*----------------------------------------------------------------------*
      Do i=1,nTri
         Fock(i)=0.5d0*(FockA(i)+FockB(i))
      End Do
*     iOffTri=0
*     Do iSym=1,nSym
*        Call TriPrt('natouhf: Fock A','(20f10.4)',
*    &      FockA(1+iOffTri),nBas(iSym))
*        Call TriPrt('natouhf: Fock B','(20f10.4)',
*    &      FockB(1+iOffTri),nBas(iSym))
*        Call TriPrt('natouhf: Average Fock','(20f10.4)',
*    &      FockB(1+iOffTri),nBas(iSym))
*        iOffTri=iOffTri+nBas(iSym)*(nBas(iSym)+1)/2
*     End Do
*----------------------------------------------------------------------*
* Form density in MO basis: C(t)SDSC                                   *
*----------------------------------------------------------------------*
      iOffSqr=0
      iOffTri=0
      iOffCMO=0
      iOffEta=0
      Do iSym=1,nSym
         If(nBas(iSym).le.0) Goto 200
*
* Compute C(t)S (=aux1), but square S first
*
         Call Square(Ovl(iOfftri+1),Smat,1,nBas(iSym),nBas(iSym))
*        Call RecPrt('natouhf: Smat','(12f12.6)',Smat,
*    &               nBas(iSym),nBas(iSym))
         Call DGEMM_('T','N',
     &               nOrb(iSym),nBas(iSym),nBas(iSym),
     &               1.0d0,CMO(iOffCMO+1),nBas(iSym),
     &                     Smat,nBas(iSym),
     &               0.0d0,Aux1,nOrb(iSym))
*        Call RecPrt('natouhf: C(t)S','(12f12.6)',Aux1,
*    &               nOrb(iSym),nBas(iSym))
*
* Compute C(t)SD (=aux3), but first unfold D (=aux2)
*
         Call Dsq(Dens(1+iOffTri),Aux2,1,nBas(iSym),nBas(iSym))
*        Call RecPrt('natouhf: Density','(12f12.6)',Aux2,
*    &               nBas(iSym),nBas(iSym))
         Call DGEMM_('N','N',
     &               nOrb(iSym),nBas(iSym),nBas(iSym),
     &               1.0d0,Aux1,nOrb(iSym),
     &                     Aux2,nBas(iSym),
     &               0.0d0,Aux3,nOrb(iSym))
*        Call RecPrt('natouhf: C(t)SD','(12f12.6)',Aux3,
*    &               nOrb(iSym),nBas(iSym))
*
* Compute C(t)SDS (=aux1)
*
         Call DGEMM_('N','N',
     &               nOrb(iSym),nBas(iSym),nBas(iSym),
     &               1.0d0,Aux3,nOrb(iSym),
     &                     Smat,nBas(iSym),
     &               0.0d0,Aux1,nOrb(iSym))
*        Call RecPrt('natouhf: C(t)SDS','(12f12.6)',Aux1,
*    &               nOrb(iSym),nBas(iSym))
*
* Compute C(t)SDSC (=aux2)
*
         Call MxMt(Aux1,   1,nOrb(iSym),
     &             CMO(iOffCMO+1), 1,nBas(iSym),
     &             Aux2,
     &             nOrb(iSym),nBas(iSym))
*        Call TriPrt('natouhf: C(t)SDSC','(12f12.6)',
*    &           Aux2,nOrb(iSym))
*
* Shift diagonal slightly
*
         tmp=1.0d-6
         Do iOrb=1,nOrb(iSym)
            indx=iOrb*(iOrb+1)/2
            Aux2(indx)=Aux2(indx)+tmp
            tmp=0.5d0*tmp
         End Do
*        Call TriPrt('natouhf: C(t)SDSC + shift','(12f12.6)',
*    &           Aux2,nOrb(iSym))
*
* Diagonalize
*
         Call NIdiag(Aux2,Nato(iOffCMO+1),
     &               nOrb(iSym),nBas(iSym),0)
         Call Pickup(Aux2,Eta(iOffEta+1),nOrb(iSym))
         Do i=1,nOrb(iSym)
            Eta(iOffEta+i)=-Eta(iOffEta+i)
         End Do
         Call Sort(Eta(iOffEta+1),Nato(iOffCMO+1),
     &             nOrb(iSym),nBas(iSym))
         Do i=1,nOrb(iSym)
            Eta(iOffEta+i)=-Eta(iOffEta+i)
         End Do
200      Continue
         iOffTri=iOffTri+nBas(iSym)*(nBas(iSym)+1)/2
         iOffSqr=iOffSqr+nBas(iSym)*nBas(iSym)
         iOffEta=iOffEta+nOrb(iSym)
         iOffCMO=iOffCMO+nBas(iSym)*nOrb(iSym)
      End Do
*----------------------------------------------------------------------*
* Compute diagonal of average Fock matrix                              *
*----------------------------------------------------------------------*
      Call MkEorb_(Fock,nBT,Nato,nBB,Eps,nnB,nSym,nBas,nOrb)
*     Call PriMO('natouhf: UHF nato',.true.,.true.,-1.0d0,1.0d6,
*    &           nSym,nBas,nOrb,Name,
*    &           Eps,Eta,Nato,3)
*----------------------------------------------------------------------*
* Deallocate arrays                                                    *
*----------------------------------------------------------------------*
      Call mma_deallocate(Aux3)
      Call mma_deallocate(Aux2)
      Call mma_deallocate(Aux1)
      Call mma_deallocate(SMat)
      Call mma_deallocate(Fock)
      Call mma_deallocate(Dens)
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Return
      End
