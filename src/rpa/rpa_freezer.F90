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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************
      Subroutine RPA_Freezer()
!
!     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
!     Figure out symmetry distribution of frozen occupied orbitals.
!
      Implicit None
#include "rpa_config.fh"
#include "rpa_data.fh"
#include "WrkSpc.fh"

      Integer  RPA_iUHF
      External RPA_iUHF

      Logical Freeze
      Logical Prnt

      Integer iUHF
      Integer iSym
      Integer iSpin
      Integer ip_Fro, l_Fro

      ! set restricted(1)/unrestricted(2)
      iUHF=RPA_iUHF()

      ! freeze orbitals (if requested)
      iSpin=1
      Freeze=nFreeze(iSpin).gt.0
      Do While (.not.Freeze .and. iSpin.lt.iUHF)
         iSpin=iSpin+1
         Freeze=Freeze.or.nFreeze(iSpin).gt.0
      End Do
      If (Freeze) Then
         Prnt=iPrint.ge.0
         l_Fro=nSym
         Call GetMem('OccFrz','Allo','Inte',ip_Fro,l_Fro)
         Do iSpin=1,iUHF
            If (nFreeze(iSpin).gt.0) Then
               Call RPA_Frz(nFreeze(iSpin),Prnt,nSym,                   &
     &                      Work(ip_OccEn(iSpin)),                      &
     &                      nFro(1,iSpin),nOcc(1,iSpin),                &
     &                      iWork(ip_Fro))
               Do iSym=1,nSym
                  nFro(iSym,iSpin)=nFro(iSym,iSpin)+iWork(ip_Fro-1+iSym)
               End Do
            End If
         End Do
         Call GetMem('OccFrz','Free','Inte',ip_Fro,l_Fro)
      End If

      ! correct number of active occupied orbitals
      Do iSpin=1,iUHF
         Do iSym=1,nSym
            nOcc(iSym,iSpin)=nOcc(iSym,iSpin)-nFro(iSym,iSpin)
         End Do
      End Do

      End
