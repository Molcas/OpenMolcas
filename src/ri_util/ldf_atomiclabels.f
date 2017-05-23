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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_SetAtomicLabels()
C
C     Thomas Bondo Pedersen, March 2011.
C
C     Purpose: set atomic labels in ldf_atomiclabels.fh
C              To unset (deallocation), call LDF_UnsetAtomicLabels().
C
C     After this routine, call LDF_GetAtomicLabel(A,Label) to retrieve
C     the four-character label Label representing atom A.
C
      Implicit None
#include "WrkSpc.fh"
#include "ldf_atomiclabels.fh"
#include "localdf_bas.fh"

      Integer  LDF_nAtom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nAtom, LDF_nShell_Atom, LDF_lShell_Atom
      Integer  LDF_GetLenIn4

      Integer LenIn4
      Integer ip_Temp, l_Temp
      Integer ip_SB, l_SB
      Integer n
      Integer iShell
      Integer A
      Integer nS
      Integer j
      Integer ipAL
      Integer ipT
#if defined (_DEBUG_)
      Integer ip
      Integer iS
      Integer ii
#endif

      Integer i
      Integer nBasSh
      Integer iSB
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iSB(i)=iWork(ip_SB-1+i)

      If (AtomicLabelsSet) Return
      l_AtomicLabels=4*LDF_nAtom()
      Call GetMem('LDFALbl','Allo','Char',ip_AtomicLabels,
     &                                     l_AtomicLabels)
      LenIn4=LDF_GetLenIn4()
      If (LenIn4.lt.4) Then
         Call WarningMessage(2,'LDF_SetAtomicLabels: LenIn4 < 4')
         Call LDF_Quit(1)
      End If
      l_Temp=LenIn4*nBas_Valence
      Call GetMem('LDFALTmp','Allo','Char',ip_Temp,l_Temp)
      Call Get_cArray('Unique Basis Names',cWork(ip_Temp),l_Temp)
      l_SB=nShell_Valence
      Call GetMem('LDFALSB','Allo','Inte',ip_SB,l_SB)
      n=0
      Do iShell=1,nShell_Valence
         iWork(ip_SB-1+iShell)=n
         n=n+nBasSh(iShell)
      End Do
      If (n.ne.nBas_Valence) Then
         Call WarningMessage(2,'LDF_SetAtomicLabels: n != nBas_Valence')
         Call LDF_Quit(1)
      End If
      Do A=1,LDF_nAtom()
         ipAL=ip_AtomicLabels+4*(A-1)
         nS=LDF_nShell_Atom(A)
         If (nS.gt.0) Then
            iShell=iWork(LDF_lShell_Atom(A))
            ipT=ip_Temp+LenIn4*iSB(iShell)
            Do j=0,3
               cWork(ipAL+j)=cWork(ipT+j)
            End Do
#if defined (_DEBUG_)
            ip=LDF_lShell_Atom(A)-1
            Do iS=1,nS
               iShell=iWork(ip+iS)
               Do ii=0,nBasSh(iShell)-1
                  ipT=ip_Temp+LenIn4*(iSB(iShell)+ii)
                  Do j=0,3
                     If (cWork(ipT+j).ne.cWork(ipAL+j)) Then
                        Call WarningMessage(2,
     &                        'LDF_SetAtomicLabels: Bfn label mismatch')
                        Write(6,'(A,I10)') 'Atom=',A
                        Call LDF_Quit(1)
                     End If
                  End Do
               End Do
            End Do
#endif
         Else
            Call WarningMessage(2,'LDF_SetAtomicLabels: nS < 1')
            Write(6,'(A,I10)') 'Atom=',A
            Call LDF_Quit(1)
         End If
      End Do
      Call GetMem('LDFALSB','Free','Inte',ip_SB,l_SB)
      Call GetMem('LDFALTmp','Free','Char',ip_Temp,l_Temp)
      AtomicLabelsSet=.True.

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_GetAtomicLabel(A,Label)
C
C     Thomas Bondo Pedersen, March 2011.
C
C     Purpose: return label representing atom A.
C              'NONE' is returned if labels are not set.
C
      Implicit None
      Integer A
      Character*4 Label
#include "WrkSpc.fh"
#include "ldf_atomiclabels.fh"

      Integer i

      If (.not.AtomicLabelsSet) Then
         Label='NONE'
      Else
         Write(Label,'(4A1)') (cWork(ip_AtomicLabels+4*(A-1)+i),i=0,3)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_UnsetAtomicLabels()
C
C     Thomas Bondo Pedersen, March 2011.
C
C     Purpose: unset atomic labels in ldf_atomiclabels.fh
C
      Implicit None
#include "ldf_atomiclabels.fh"

      If (.not.AtomicLabelsSet) Return
      Call GetMem('LDFALbl','Free','Char',ip_AtomicLabels,
     &                                     l_AtomicLabels)
      l_AtomicLabels=0
      ip_AtomicLabels=0
      AtomicLabelsSet=.False.

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_PrintAtomicLabels()
      Implicit None
      Integer  LDF_nAtom
      External LDF_nAtom

      Character*4 Label
      Integer A

      Do A=1,LDF_nAtom()
         Call LDF_GetAtomicLabel(A,Label)
         Write(6,'(A,I10,A,A)') 'Atom=',A,' Label=',Label
      End Do
      Call xFlush(6)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Integer Function LDF_GetLenIn4()
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      LDF_GetLenIn4=LENIN4
      End
