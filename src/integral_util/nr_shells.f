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
! Copyright (C) 1998, Roland Lindh                                     *
!***********************************************************************
      Subroutine Nr_Shells(nSkal)
!***********************************************************************
!                                                                      *
!     Object: to compute the number of unique shells in the input.     *
!                                                                      *
!     Author: Roland Lindh, Chemical Physics, University of Lund,      *
!             Sweden. January '98.                                     *
!***********************************************************************
      use Basis_Info
      use BasisMode
      Implicit None
      Integer, Intent(Out):: nSkal

      Integer iCnttp, nTest, iCnt, iAng, iShll, nExpi, nBasisi
!                                                                      *
!***********************************************************************
!                                                                      *
!     Determine the number of shells
!
      nSkal=0
      If (Basis_Mode.ne.Valence_Mode .and.
     &    Basis_Mode.ne.Auxiliary_Mode .and.
     &    Basis_Mode.ne.Fragment_Mode .and.
     &    Basis_Mode.ne.With_Auxiliary_Mode .and.
     &    Basis_Mode.ne.With_Fragment_Mode .and.
     &    Basis_Mode.ne.All_Mode) Then
         Call WarningMessage(2,'Nr_Shells: illegal Basis_Mode')
         Call Abend()
      End If
!
      Select Case (Atomic)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     Molecular set up                                                 *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Case (.False.)
      Do iCnttp = 1, nCnttp
         nTest = dbsc(iCnttp)%nVal-1
         Do iCnt = 1, dbsc(iCnttp)%nCntr
!
            Do iAng=0, nTest
               iShll = dbsc(iCnttp)%iVal + iAng
               nExpi=Shells(iShll)%nExp
               If (nExpi.eq.0) Cycle
               nBasisi=Shells(iShll)%nBasis
               If (nBasisi.eq.0) Cycle
!
               If (Basis_Mode.eq.Valence_Mode .and.
     &             (Shells(iShll)%Aux.or.Shells(iShll)%Frag)) Cycle
               If (Basis_Mode.eq.Auxiliary_Mode .and.
     &             .Not.Shells(iShll)%Aux) Cycle
               If (Basis_Mode.eq.Fragment_Mode .and.
     &             .Not.Shells(iShll)%Frag) Cycle
               If (Basis_Mode.eq.With_Auxiliary_Mode .and.
     &             Shells(iShll)%Frag) Cycle
               If (Basis_Mode.eq.With_Fragment_Mode .and.
     &             Shells(iShll)%Aux) Cycle
               nSkal = nSkal + 1
              End Do                     ! iAng
         End Do                          ! iCnt
      End Do                             ! iCnttp
!
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     Atomic set up                                                    *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      Case (.True.)
!
      Do iCnttp = kCnttp, lCnttp
      nTest = dbsc(iCnttp)%nVal-1
      Do iAng=0, nTest
         iShll = dbsc(iCnttp)%iVal + iAng
         nExpi=Shells(iShll)%nExp
         If (nExpi.eq.0) Cycle
         nBasisi=Shells(iShll)%nBasis
         If (nBasisi.eq.0) Cycle
!
         If (Shells(iShll)%Frag) Cycle
         nSkal = nSkal + 1
!
      End Do                     ! iAng
      End Do
      If (dbsc(kCnttp)%Aux) nSkal=nSkal+1 ! Add dummy shell
      End Select
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
      End Subroutine Nr_Shells
