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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
      Subroutine MkEorb(Fock,nFck,CMO,nCMO,Eorb,nEorb,nSym,nBas,nOrb,nD)
      Implicit None
      Integer nFck, nCMO, nEOrb, nD
      Real*8 Fock(nFck,nD)
      Real*8 CMO(nCMO,nD)
      Real*8 EOrb(nEOrb,nD)
      Integer nSym
      Integer nBas(nSym)
      Integer nOrb(nSym)
*
      Integer iD
*
      Do iD = 1, nD
         Call MkEorb_(Fock(1,iD),nFck,CMO(1,iD),nCMO,Eorb(1,iD),nEorb,
     &                nSym,nBas,nOrb)
      End Do
*
      Return
      End
      Subroutine MkEorb_(Fock,nFck,CMO,nCMO,Eorb,nEorb,nSym,nBas,nOrb)
************************************************************************
*                                                                      *
*  This routine calculates the diagonal elements of the MO Fock matrix *
*  (orbital energies).                                                 *
*                                                                      *
*  Input:                                                              *
*    Fock    Fock matrix in AO basis                                   *
*    CMO     Orbitals                                                  *
*                                                                      *
*  Output:                                                             *
*    Eorb    Orbital energies.                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      Implicit None
*     Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
*----------------------------------------------------------------------*
* Dummy arguments.                                                     *
*----------------------------------------------------------------------*
      Integer nFck, nCMO, nEOrb
      Real*8 Fock(nFck)
      Real*8 CMO(nCMO)
      Real*8 EOrb(nEOrb)
      Integer nSym
      Integer nBas(nSym)
      Integer nOrb(nSym)
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Real*8  t
      Real*8, Dimension(:), Allocatable:: FckSqr
      Integer iSym
      Integer iBas
      Integer jBas
      Integer iOrb
      Integer MaxTri
      Integer MaxSqr
      Integer iOffTri
      Integer iOffCMO
      Integer npFckSqr
      Integer indE
      Integer indF
      Integer indx
      Integer jndx
*----------------------------------------------------------------------*
* Some preliminary setup.                                              *
*----------------------------------------------------------------------*
      MaxTri = 0
      MaxSqr = 0
      Do iSym=1,nSym
         MaxTri=Max(MaxTri,nBas(iSym)*(nBas(iSym)+1)/2)
         MaxSqr=Max(MaxSqr,nBas(iSym)*nBas(iSym))
      End Do
*----------------------------------------------------------------------*
* Allocate matrices.                                                   *
*----------------------------------------------------------------------*
      npFckSqr=MaxSqr
      Call mma_allocate(FckSqr,npFckSqr,Label='FckSqr')
*----------------------------------------------------------------------*
* Do compute orbital energies                                          *
*----------------------------------------------------------------------*
      iOffTri=0
      iOffCMO=0
      indE=1
      Do iSym=1,nSym
         If(nOrb(iSym).gt.0) Then
            Call Square(Fock(1+iOffTri),FckSqr,
     &         1,nBas(iSym),nBas(iSym))
*           Call RecPrt('mkeor: Fock matrix','(20F10.4)',
*    &         FckSqr,nBas(iSym),nBas(iSym))
*           Call RecPrt('mkeor: CMO','(20F10.4)',
*    &         CMO(1+iOffCMO),nBas(iSym),nOrb(iSym))
            Do iOrb=1,nOrb(iSym)
               t=0.0d0
               indF=1
               Do iBas=1,nBas(iSym)
                  Do jBas=1,nBas(iSym)
                     indx=iBas+(iOrb-1)*nBas(iSym)+iOffCMO
                     jndx=jBas+(iOrb-1)*nBas(iSym)+iOffCMO
                     t=t+CMO(indx)*CMO(jndx)*FckSqr(indF)
                     indF=indF+1
                  End Do
               End Do
               Eorb(indE)=t
               indE=indE+1
            End Do
         End If
         iOffTri=iOffTri+nBas(iSym)*(nBas(iSym)+1)/2
         iOffCMO=iOffCMO+nBas(iSym)*nOrb(iSym)
      End Do
*----------------------------------------------------------------------*
* Debug print                                                          *
*----------------------------------------------------------------------*
*     indE=1
*     Do iSym=1,nSym
*        If(nOrb(iSym).gt.0) Then
*           Call RecPrt('mkeor: Orbital energies','(20F10.4)',
*    &                  Eorb(indE),1,nOrb(iSym))
*        End If
*        indE=indE+nOrb(iSym)
*     End Do
*----------------------------------------------------------------------*
* Deallocate matrices.                                                 *
*----------------------------------------------------------------------*
      Call mma_deallocate(FckSqr)
*----------------------------------------------------------------------*
* Done.                                                                *
*----------------------------------------------------------------------*
      Return
      End
