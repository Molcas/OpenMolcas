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
      SubRoutine DFroz(Dlt,nDlt,CMO,nCMO,OccNo)
************************************************************************
*                                                                      *
*     purpose: Compute contribution to the density matrix from frozen  *
*              orbitals                                                *
*                                                                      *
*     output:                                                          *
*       Dlt     : density contribution                                 *
*                                                                      *
*     called from: ModFck                                              *
*                                                                      *
*     calls to: DOne                                                   *
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
      Real*8 Dlt(nDlt),CMO(nCMO),OccNo(*)
*
      Real*8, Dimension(:), Allocatable:: NewOcc
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*---- Allocate memory for new occupation numbers
      Call mma_allocate(NewOcc,nnB,Label='NewOcc')
*
*---- Set proper occupation numbers
      iStrtN = 0
      iStrtO = 0
      Do iSym = 1, nSym
         Do iOrb = 1, nOrb(iSym)
            NewOcc(iStrtN + iOrb) = Zero
            If (iOrb.le.nFro(iSym))
     &          NewOcc(iStrtN + iOrb) = OccNo(iStrtO + iOrb)
         End Do
         iStrtN = iStrtN + nOrb(iSym)
         iStrtO = iStrtO + nOrb(iSym)
      End Do
*
*---- Compute density contribution
      Call DOne_SCF_froz(nSym,nBas,nOrb,nFrz,CMO,nCMO,NewOcc,Dlt)
*
*---- Deallocate memory for new occupation numbers
      Call mma_deallocate(NewOcc)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
*----------------------------------------------------------------------*
*
************************************************************************
      SubRoutine DOne_SCF_froz(nSym,nBas,nOrb,nFro,Cff,nCff,Occ,Dlt)
************************************************************************
*                                                                      *
*     purpose: Compute density matrix in AO basis                      *
*                                                                      *
*     input:                                                           *
*       nSym    : number of symmetries                                 *
*       nBas(i) : number of basis functions (i = 1, nSym)              *
*       nOrb(i) : number of orbitals (i = 1, nSym)                     *
*       nFro(i) : number of frozen orbitals (i = 1, nSym)              *
*       Cff     : molecular orbitals                                   *
*       Occ     : occupation numbers                                   *
*                                                                      *
*     output:                                                          *
*       Dlt     : density matrix in triangular storrage                *
*                                                                      *
*     called from: DFroz                                               *
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
*
      Implicit Real*8 (a-h,o-z)
*
      Real*8    Cff(nCff),Occ(*),Dlt(*)
      Integer nBas(nSym),nOrb(nSym),nFro(nSym)
*
#include "real.fh"
*
*---- Statement function for triangular storrage
      Ind(i,j) = i*(i - 1)/2 + j
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*
      Scale=Two
      SScale=One
      ipCff  = 0
      ipDlt  = 0
      ipOcc  = 0
      Do iSym = 1, nSym
*
         nBs = nBas(iSym)
         nOr = nOrb(iSym)
         nFr = nFro(iSym)
*
         nBF = nBs - nFr
         nOF = nOr - nFr
         lth = nBs*(nBs + 1)/2
*
         ipCff = ipCff + nBs*nFr
         Do iRow = 1, nBs
            Sum = Zero
            ij  = -1
            Do i = nFr + 1, nOr
               ij  = ij  + 1
*              If (Occ(ipOcc+i).eq.0.0D0) Go To 100
               Sum = Sum + Occ(ipOcc + i)*
     &         Cff(ipCff + iRow + ij*nBs)*
     &         Cff(ipCff + iRow + ij*nBs)
            End Do
*100        Continue
            Dlt(ipDlt + Ind(iRow,iRow)) = Sum*SScale
*
            Do iCol = 1, iRow - 1
               Sum = Zero
               ij  = -1
               Do i = nFr + 1, nOr
                  ij = ij  + 1
*                 If (Occ(ipOcc+i).eq.0.0D0) Go To 200
                  Sum = Sum + Occ(ipOcc + i)*
     &            Cff(ipCff + iRow + ij*nBs)*
     &            Cff(ipCff + iCol + ij*nBs)
               End Do
*200           Continue
               Dlt(ipDlt + Ind(iRow,iCol)) = Scale*Sum
            End Do
         End Do
*define _DEBUG_
#ifdef _DEBUG_
         Call NrmClc(Dlt(ipDlt),nBs,'DOne_SCF_froz','Dlt(ipDlt)')
#endif
*
         ipCff = ipCff + nBs*nOF
         ipDlt = ipDlt + lth
         ipOcc = ipOcc + nOr
      End Do
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
