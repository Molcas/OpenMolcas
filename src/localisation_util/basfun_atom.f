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
* Copyright (C) Yannick Carissan                                       *
************************************************************************
      SubRoutine BasFun_Atom(nBas_per_Atom,nBas_Start,Name,
     &                       nBas,nAtoms,DoPrint)
C
C     Author: Y. Carissan [put in separate subroutine by T.B. Pedersen]
C
      Implicit None
#include "Molcas.fh"
      Integer nBas, nAtoms
      Integer nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
      Character*(LENIN4) Name(nBas)
      Logical DoPrint

      Character*11 SecNam
      Parameter (SecNam = 'BasFun_Atom')

      Integer iAt, iAt1, nBasAt, iBas, iCount

      Character*(LENIN)  Lbl, LblOld
      Character*80 Txt, Formt

C     Counters.
C     ---------

      iAt    = 1
      nBasAt = 1
      LblOld = Name(1)(1:LENIN)
      Do iBas = 2,nBas
         Lbl = Name(iBas)(1:LENIN)
         If (Lbl .ne. LblOld) Then
            nBas_per_Atom(iAt) = nBasAt
            iAt    = iAt + 1
            nBasAt = 0
            LblOld = Lbl
         End If
         nBasAt = nBasAt + 1
      End Do
      nBas_per_Atom(iAt) = nBasAt

      If (iAt .ne. nAtoms) Then ! centers without basis functions
         iAt1 = iAt + 1
         Do iAt = iAt1,nAtoms
            nBas_per_Atom(iAt) = 0
         End Do
      End If

C     Offsets.
C     --------

      iCount = 0
      Do iAt = 1,nAtoms
         nBas_Start(iAt) = iCount + 1
         iCount = iCount + nBas_per_Atom(iAt)
      End Do
      If (iCount .ne. nBas) Then
         Write(Txt,'(A,I9,A,I9)') 'iCount =',iCount,'  nBas =',nBas
         Call SysAbendMsg(SecNam,'iCount.NE.nBas',Txt)
      End If

C     Print.
C     ------

      If (DoPrint) Then
         Write(Formt,'(3(a6,i3,a5))') '(/,a6,',nAtoms,'i5,/,',
     &                                '   a6,',nAtoms,'i5,/,',
     &                                '   a6,',nAtoms,'i5)'
         Write(6,Formt) 'Atom  ',(iAt,iAt=1,nAtoms),
     &                  'Start ',(nBas_Start(iAt),iAt=1,nAtoms),
     &                  'nBas  ',(nBas_per_Atom(iAt),iAt=1,nAtoms)
      End If

      End
