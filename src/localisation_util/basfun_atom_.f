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
! Copyright (C) Yannick Carissan                                       *
!               Thomas Bondo Pedersen                                  *
!               Francesco Aquilante                                    *
!***********************************************************************
      SubRoutine BasFun_Atom_(nBas_per_Atom,nBas_Start,Name,            &
     &                        jBas,nBas,nAtoms,DoPrint)
!
!     Author: Y. Carissan / T. B. Pedersen
!             [adapted to cases with symmetry by F. Aquilante]
!
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      Integer jBas, nBas, nAtoms
      Integer nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
      Character*(LENIN8) Name(nBas)
      Character*(LENIN) AtName(nAtoms)
      Logical DoPrint

      Character*12 SecNam
      Parameter (SecNam = 'BasFun_Atom_')

      Integer iAt, kBas, iCount

      Character*(LENIN)  Lbl
      Character*80 Txt, Formt


!     Counters.
!     ---------

! IFG: To count basis functions per atom, we need a list of atom names,
!      since there is no guarantee all atoms will be present in a give irrep
      Call Get_cArray('Unique Atom Names',AtName,(LENIN)*nAtoms)

      kBas = jBas
      Do iAt = 1, nAtoms
         nBas_per_Atom(iAt) = 0
         Lbl = AtName(iAt)
         Do While ((Name(kBas)(1:LENIN).eq.Lbl) .and. (kBas.le.nBas))
            nBas_per_Atom(iAt) = nBas_per_Atom(iAt) + 1
            kBas = kBas + 1
         End Do
      End Do

!     Offsets.
!     --------

      iCount = 0
      Do iAt = 1,nAtoms
         nBas_Start(iAt) = iCount + 1
         iCount = iCount + nBas_per_Atom(iAt)
      End Do
      jCount = iCount + jBas - 1
      If (jCount .ne. nBas) Then
         Write(Txt,'(A,I9,A,I9)') 'jCount =',jCount,'  nBas =',nBas
         Call SysAbendMsg(SecNam,'jCount.NE.nBas',Txt)
      End If

!     Print.
!     ------

      If (DoPrint) Then
         Write(Formt,'(3(a6,i3,a5))') '(/,a6,',nAtoms,'i5,/,',          &
     &                                '   a6,',nAtoms,'i5,/,',          &
     &                                '   a6,',nAtoms,'i5)'
         Write(6,Formt) 'Atom  ',(iAt,iAt=1,nAtoms),                    &
     &                  'Start ',(nBas_Start(iAt),iAt=1,nAtoms),        &
     &                  'nBas  ',(nBas_per_Atom(iAt),iAt=1,nAtoms)
      End If

      End
