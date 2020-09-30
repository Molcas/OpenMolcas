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
      SubRoutine MODens(Dens,Ovlp,nDO,NumD,CMO,nCMO,nD)
************************************************************************
*                                                                      *
*     purpose: Compute density matrix in molecular orbital basis       *
*                                                                      *
*     input:                                                           *
*       Dens    : density matrix differences in AO basis (nDO,NumD)    *
*       Ovlp    : overlap in AO basis of length of length nDO          *
*       CMO     : molecular orbitals coefficients of length nCMO       *
*                                                                      *
*     called from: WfCtl                                               *
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
      Real*8 Dens(nDO,nD,NumD),Ovlp(nDO),CMO(nCMO,nD)
      Real*8, Dimension(:), Allocatable:: DnsS, OvlS, DMoO, Aux1, Aux2
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*define _DEBUG_
#ifdef _DEBUG_
#endif
*
*---- Allocate memory for squared density matrix
      Call mma_allocate(DnsS,MaxBas**2,Label='DnsS')
*
*---- Allocate memory for squared overlap matrix
      Call mma_allocate(OvlS,MaxBas**2,Label='OvlS')
*
*---- Allocate memory for density matrix in MO basis
      Call mma_allocate(DMoO,MaxOrb*(MaxOrb+1)/2,Label='DMoO')
*
*---- Allocate memory for auxiliary matrix
      Call mma_allocate(Aux1,MaxBxO,Label='Aux1')
*
*---- Allocate memory for auxiliary matrix
      Call mma_allocate(Aux2,MaxBxO,Label='Aux2')
*
*
      DMOMax = Zero
      Do jD = 1, nD
*
#ifdef _DEBUG_
         Call NrmClc(Dens(1,jD,nDens),nBT,'MoDens','D in AO   ')
         Call NrmClc(Ovlp         ,nBT,'MoDens','Overlap   ')
         Call NrmClc(CMO(1,jD)    ,nBO,'MoDens','CMOs      ')
         Write (6,*) 'nOcc=',(nOcc(i,jD),i=1,nSym)
C        Write (6,'(F16.8)') DXot(MaxBxO,CMO(1,jD),1,CMO(1,jD),1)
#endif
         it   = 1
         id   = 1
         iOvl = 1
         Do iSym = 1, nSym
*
            iiBB = nBas(iSym)*nBas(iSym)
            iiBO = nBas(iSym)*nOrb(iSym)
            iiBT = nBas(iSym)*(nBas(iSym) + 1)/2
*
            If (nOcc(iSym,jD).gt.0.or.(Teee.and.nBas(iSym).gt.0)) Then
               Call DSq(Dens(id,jD,nDens),DnsS,1,nBas(iSym),nBas(iSym))
               Call Square(Ovlp(iOvl),OvlS,1,nBas(iSym),nBas(iSym))
*
               Call DGEMM_('N','N',
     &                     nBas(iSym),nOrb(iSym),nBas(iSym),
     &                     1.0d0,OvlS,nBas(iSym),
     &                           CMO(it,jD),nBas(iSym),
     &                     0.0d0,Aux1,nBas(iSym))
               Call DGEMM_('N','N',
     &                     nBas(iSym),nOrb(iSym),nBas(iSym),
     &                     1.0d0,DnsS,nBas(iSym),
     &                           Aux1,nBas(iSym),
     &                     0.0d0,Aux2,nBas(iSym))
               Call DGEMM_('N','N',
     &                     nBas(iSym),nOrb(iSym),nBas(iSym),
     &                     1.0d0,OvlS,nBas(iSym),
     &                           Aux2,nBas(iSym),
     &                     0.0d0,Aux1,nBas(iSym))
               Call   MxMt(CMO(it,jD),     nBas(iSym),1,
     &                     Aux1,1,nBas(iSym),
     &                     DMoO,
     &                     nOrb(iSym),nBas(iSym))
*
*              Call TriPrt('D(mo)','(8F12.6)',DMoO,nOrb(iSym))
               Do i = nOcc(iSym,jD)+1 , nOrb(iSym)
                  Do j = 1, nOcc(iSym,jD)
                     DMOMax = Max(DMOMax,
     &                        DBLE(nD)*Abs(DMoO(i*(i-1)/2+j)))
                  End Do
               End Do
            End If
*
            it   = it   + iiBO
            id   = id   + iiBT
            iOvl = iOvl + iiBT
         End Do
*
      End Do
*
*---- Deallocate memory
      Call mma_deallocate(Aux2)
      Call mma_deallocate(Aux1)
      Call mma_deallocate(DMoO)
      Call mma_deallocate(OvlS)
      Call mma_deallocate(DnsS)
*
#ifdef _DEBUG_
      Write(6,*)' DMOMax in MODens',DMOMax
#endif
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(14) = TimFld(14) + (Cpu2 - Cpu1)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
