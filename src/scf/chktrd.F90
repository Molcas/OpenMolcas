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
! Copyright (C) 1996, Martin Schuetz                                   *
!***********************************************************************
      SubRoutine ChkTrD(nSym,nBas,nOrb,Occ,nOcc,Dlt,nDlt)
!***********************************************************************
!                                                                      *
!     purpose: Compute trace of density matrix and compare with sum    *
!              over occupation numbers...                              *
!                                                                      *
!     input:                                                           *
!       nSym    : number of symmetries                                 *
!       nBas(i) : number of basis functions (i = 1, nSym)              *
!       Occ     : occupation numbers                                   *
!       Dlt     : density matrix in triangular storrage                *
!                                                                      *
!***********************************************************************
!
      use SCF_Arrays, only: Ovrlp
      use Constants, only: Zero, One
      Implicit None
!
!     declaration of subroutine parameters...
      Integer nSym, nOcc, nDlt
      Real*8 Occ(nOcc), Dlt(nDlt)
      Integer nBas(nSym),nOrb(nSym)
!
!     declaration of some local variables...
      Integer iOr, ipDlt, ipOcc, ipOvl, iSym, lth, nBs, nOr
      Real*8 :: ThrDif=1.0d-7
      Real*8 :: Scale, SumOcc, TrDns
      Real*8, External:: DDot_
!
      ipDlt = 1
      ipOvl = 1
      ipOcc = 0
      Scale= One
      Do iSym = 1, nSym
        nBs = nBas(iSym)
        nOr = nOrb(iSym)
        lth = nBs*(nBs + 1)/2
!       count occupation number...
        SumOcc=Zero
        Do iOr = 1, nOr
          SumOcc=SumOcc+Occ(ipOcc+iOr)*Scale
        End Do
!       do trace of PS for symmetry block...
        TrDns=DDOT_(lth,Dlt(ipDlt),1,Ovrlp(ipOvl),1)
        ipDlt = ipDlt+lth
        ipOvl = ipOvl+lth
        ipOcc = ipOcc+nOr
        If (Abs(SumOcc-TrDns)>ThrDif) Then
        Write (6,*) Abs(SumOcc-TrDns)
!         print Warning...
          Call WarningMessage(1,'WARNING: trace of density is inconsistent with occupation !')
          Write(6,'(A,I1,A,3F12.7)') 'SymBlock: ',iSym,' deviation: ', SumOcc-TrDns,SumOcc,TrDns
        End If
      End Do

      End SubRoutine ChkTrD
