!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************
      SubRoutine DFroz(Dlt,nDlt,CMO,nCMO,OccNo)
!***********************************************************************
!                                                                      *
!     purpose: Compute contribution to the density matrix from frozen  *
!              orbitals                                                *
!                                                                      *
!***********************************************************************
      use InfSCF, only: nBas, nnB, nSym, nFro, nOrb
      use Constants, only: Zero, One, Two
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
!
      Integer nDlt, nCMO
      Real*8 Dlt(nDlt),CMO(nCMO)
      Integer OccNo(*)
!
      Integer iOrb, iStrtN, iStrtO, iSym
      Real*8, Dimension(:), Allocatable:: NewOcc
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!---- Allocate memory for new occupation numbers
      Call mma_allocate(NewOcc,nnB,Label='NewOcc')
!
!---- Set proper occupation numbers
      iStrtN = 0
      iStrtO = 0
      Do iSym = 1, nSym
         Do iOrb = 1, nOrb(iSym)
            NewOcc(iStrtN + iOrb) = Zero
            If (iOrb.le.nFro(iSym)) NewOcc(iStrtN + iOrb) = OccNo(iStrtO + iOrb)
         End Do
         iStrtN = iStrtN + nOrb(iSym)
         iStrtO = iStrtO + nOrb(iSym)
      End Do
!
!---- Compute density contribution
      Call DOne_SCF_froz(CMO,nCMO,NewOcc,nnB)
!
!---- Deallocate memory for new occupation numbers
      Call mma_deallocate(NewOcc)
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Return
      Contains
!----------------------------------------------------------------------*
!
!***********************************************************************
      SubRoutine DOne_SCF_froz(Cff, nCff, Occ,nnB)
!***********************************************************************
!                                                                      *
!     purpose: Compute density matrix in AO basis                      *
!                                                                      *
!     input:                                                           *
!       Cff     : molecular orbitals                                   *
!       Occ     : occupation numbers                                   *
!                                                                      *
!     output:                                                          *
!       Dlt     : density matrix in triangular storrage                *
!                                                                      *
!     called from: DFroz                                               *
!                                                                      *
!***********************************************************************
!
      Implicit None
!
      Integer nCff, nnB
      Real*8  Cff(nCff), Occ(nnB)

      Integer i, j, iCol, ij, ipCff, ipDlt, ipOcc, iRow, lth, nBs, nFr, nOr, nOF, iSym
      Integer Ind
      Real*8 Scale, SScale, Sum
!
!---- Statement function for triangular storrage
      Ind(i,j) = i*(i - 1)/2 + j
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!
      Scale=Two
      SScale=One
      ipCff  = 0
      ipDlt  = 0
      ipOcc  = 0
      Do iSym = 1, nSym
!
         nBs = nBas(iSym)
         nOr = nOrb(iSym)
         nFr = nFro(iSym)
!
         nOF = nOr - nFr
         lth = nBs*(nBs + 1)/2
!
         ipCff = ipCff + nBs*nFr
         Do iRow = 1, nBs
            Sum = Zero
            ij  = -1
            Do i = nFr + 1, nOr
               ij  = ij  + 1
!              If (Occ(ipOcc+i).eq.Zero) Go To 100
               Sum = Sum + Occ(ipOcc + i)* Cff(ipCff + iRow + ij*nBs)* Cff(ipCff + iRow + ij*nBs)
            End Do
!100        Continue
            Dlt(ipDlt + Ind(iRow,iRow)) = Sum*SScale
!
            Do iCol = 1, iRow - 1
               Sum = Zero
               ij  = -1
               Do i = nFr + 1, nOr
                  ij = ij  + 1
!                 If (Occ(ipOcc+i).eq.Zero) Go To 200
                  Sum = Sum + Occ(ipOcc + i)* Cff(ipCff + iRow + ij*nBs)* Cff(ipCff + iCol + ij*nBs)
               End Do
!200           Continue
               Dlt(ipDlt + Ind(iRow,iCol)) = Scale*Sum
            End Do
         End Do
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
         Call NrmClc(Dlt(ipDlt),nBs,'DOne_SCF_froz','Dlt(ipDlt)')
#endif
!
         ipCff = ipCff + nBs*nOF
         ipDlt = ipDlt + lth
         ipOcc = ipOcc + nOr
      End Do
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Return
      End SubRoutine DOne_SCF_froz

      End SubRoutine DFroz
