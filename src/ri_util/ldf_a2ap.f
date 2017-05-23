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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_SetA2AP()
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Set up atom to atom pair mapping
C              (stored in ldf_a2ap.fh).
C              To unset, call LDF_UnsetA2AP.
C
C     The mapping can be used as follows:
C


C
C     Integer  LDF_nAtom
C     External LDF_nAtom
C
C     Integer A, AB, iAB
C     Integer NumberOfAtomPairsContainingA
C     Integer PointerToListOfAtomPairsContainingA
C
C     Integer i, j
C     Integer A2AP
C     A2AP(i,j)=iWork(ip_A2AP-1+2*(j-1)+i)
C
C     Do A=1,LDF_nAtom()
C        NumberOfAtomPairsContainingA=A2AP(1,A)
C        PointerToListOfAtomPairsContainingA=A2AP(2,A)
C        Do iAB=1,NumberOfAtomPairsContainingA
C           AB=iWork(PointerToListOfAtomPairsContainingA-1+iAB)
C           ....do stuff for atom pair AB....
C        End Do
C     End Do
C
      Implicit None
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_a2ap.fh"

      Character*11 SecNam
      Parameter (SecNam='LDF_SetA2AP')

      Integer  LDF_nAtom
      External LDF_nAtom

#if defined (_DEBUG_)
      Integer n, iAB
#endif
      Integer nAtoms
      Integer ip, l
      Integer AB, A, B

      Character*8 Label

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      If (l_A2AP.gt.0) Then
#if defined (_DEBUG_)
         Call WarningMessage(0,
     &                SecNam//': atom to atom pair map already set up!')
#endif
         Return
      End If

      nAtoms=LDF_nAtom()

      l_A2AP=2*nAtoms
      Call GetMem('LDFA2AP','Allo','Inte',ip_A2AP,l_A2AP)
      Call iZero(iWork(ip_A2AP),l_A2AP)

      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)-1
         B=AP_Atoms(2,AB)-1
         iWork(ip_A2AP+2*A)=iWork(ip_A2AP+2*A)+1
         If (B.ne.A) Then
            iWork(ip_A2AP+2*B)=iWork(ip_A2AP+2*B)+1
         End If
      End Do

      Do A=0,nAtoms-1
         l=iWork(ip_A2AP+2*A)
         If (l.lt.1) Then
            Call WarningMessage(2,
     &                SecNam//': An atom is not part of any atom pair!')
            Call LDF_Quit(1)
         Else
            Write(Label,'(A3,I5.5)') 'A2P',A
            Call GetMem(Label,'Allo','Inte',ip,l)
            iWork(ip_A2AP+2*A+1)=ip
         End If
      End Do

      Do A=0,nAtoms-1
         iWork(ip_A2AP+2*A)=0
      End Do
      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)-1
         B=AP_Atoms(2,AB)-1
         ip=iWork(ip_A2AP+2*A+1)
         iWork(ip+iWork(ip_A2AP+2*A))=AB
         iWork(ip_A2AP+2*A)=iWork(ip_A2AP+2*A)+1
         If (B.ne.A) Then
            ip=iWork(ip_A2AP+2*B+1)
            iWork(ip+iWork(ip_A2AP+2*B))=AB
            iWork(ip_A2AP+2*B)=iWork(ip_A2AP+2*B)+1
         End If
      End Do

#if defined (_DEBUG_)
      n=0
      Do A=1,nAtoms
         l=iWork(ip_A2AP+2*(A-1))
         If (l.lt.1) Then
            Call WarningMessage(1,
     &                SecNam//': an atom is not part of any atom pair!')
            Write(6,'(A,I9)') 'Atom=',A
            Call xFlush(6)
            n=n+1
         Else
            ip=iWork(ip_A2AP+2*(A-1)+1)
            Do iAB=1,l
               AB=iWork(ip-1+iAB)
               If (AP_Atoms(1,AB).ne.A .and. AP_Atoms(2,AB).ne.A) Then
                  n=n+1
               End If
            End Do
         End If
      End Do
      If (n.ne.0) Then
         Call WarningMessage(2,SecNam//': Error detected!')
         Call LDF_Quit(1)
      End If
#endif

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_UnsetA2AP()
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: shut down atom to atom pair mapping.
C
      Implicit None
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_a2ap.fh"

      Character*13 SecNam
      Parameter (SecNam='LDF_UnsetA2AP')

      Integer  LDF_nAtom
      External LDF_nAtom

      Character*8 Label

      Integer A, ip, l

      If (l_A2AP.lt.1) Then
#if defined (_DEBUG_)
         Call WarningMessage(0,
     &                 SecNam//' called, but A2AP has not been set up!')
         Call xFlush(6)
#endif
         Return
      End If

      Do A=0,LDF_nAtom()-1
         l=iWork(ip_A2AP+2*A)
         If (l.lt.1) Then
            Call WarningMessage(2,
     &                SecNam//': an atom is not part of any atom pair!')
            Call LDF_Quit(1)
         Else
            ip=iWork(ip_A2AP+2*A+1)
            Write(Label,'(A3,I5.5)') 'A2P',A
            Call GetMem(Label,'Free','Inte',ip,l)
         End If
      End Do
      Call GetMem('LDFA2AP','Free','Inte',ip_A2AP,l_A2AP)
      ip_A2AP=0
      l_A2AP=0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_PrintA2AP()
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: print atom to atom pair mapping.
C
      Implicit None
#include "WrkSpc.fh"
#include "ldf_a2ap.fh"

#if defined (_DEBUG_)
      Character*13 SecNam
      Parameter (SecNam='LDF_PrintA2AP')
#endif

      Integer  LDF_nAtom
      External LDF_nAtom

      Integer A, n, ip

      Integer i, j
      Integer A2AP
      A2AP(i,j)=iWork(ip_A2AP-1+2*(j-1)+i)

      If (l_A2AP.lt.1) Then
#if defined (_DEBUG_)
         Call WarningMessage(0,
     &         SecNam//' called, but atom to atom pair map not set up!')
         Call xFlush(6)
#endif
         Return
      End If

      Call Cho_Head('Atom to Atom Pair Map','-',80,6)
      Write(6,'(/,A)') '    Atom   #Pairs           List of Pairs'
      Write(6,'(118A1)') ('-',i=1,118)
      Do A=1,LDF_nAtom()
         n=A2AP(1,A)
         ip=A2AP(2,A)
         Call LDF_PrintAtomInfo_(A,n,iWork(ip))
      End Do
      Write(6,'(118A1)') ('-',i=1,118)
      Call xFlush(6)

      End
