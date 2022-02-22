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
C     Purpose: set atomic labels in ldf_atomiclabels
C              To unset (deallocation), call LDF_UnsetAtomicLabels().
C
C     After this routine, call LDF_GetAtomicLabel(A,Label) to retrieve
C     the four-character label Label representing atom A.
C
      Use ldf_atomiclabels, only: AtomicLabels
      Use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
#include "WrkSpc.fh"
#include "localdf_bas.fh"

      Integer  LDF_nAtom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nAtom, LDF_nShell_Atom, LDF_lShell_Atom
      Integer  LDF_GetLenIn8

      Integer LenIn8
      Integer l_Temp
      Integer ip_SB, l_SB
      Integer n
      Integer iShell
      Integer A
      Integer nS
      Integer j
#if defined (_DEBUGPRINT_)
      Integer ip
      Integer iS
      Integer ii
#endif

      Integer i
      Integer nBasSh
      Integer iSB
      Character, Allocatable :: Temp(:)
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iSB(i)=iWork(ip_SB-1+i)

      If (allocated(AtomicLabels)) Return
      call mma_allocate(AtomicLabels,LDF_nAtom())
      LenIn8=LDF_GetLenIn8()
      If (LenIn8.lt.8) Then
         Call WarningMessage(2,'LDF_SetAtomicLabels: LenIn8 < 8')
         Call LDF_Quit(1)
      End If
      l_Temp=LenIn8*nBas_Valence
      Call mma_allocate(Temp,l_Temp,label='LDFALTmp')
      Call Get_cArray('Unique Basis Names',Temp,l_Temp)
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
         nS=LDF_nShell_Atom(A)
         If (nS.gt.0) Then
            iShell=iWork(LDF_lShell_Atom(A))
            Do j=1,4
               AtomicLabels(A)(j:j) = Temp(LenIn8*iSB(iShell)+j)
            End Do
#if defined (_DEBUGPRINT_)
            ip=LDF_lShell_Atom(A)-1
            Do iS=1,nS
               iShell=iWork(ip+iS)
               Do ii=0,nBasSh(iShell)-1
                  Do j=1,4
                     If (Temp(LenIn8*(iSB(iShell)+ii)+j).ne.
     &                   AtomicLabels(A)(j:j)) Then
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
      call mma_deallocate(Temp)

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
      Use ldf_atomiclabels, only: AtomicLabels
      Implicit None
      Integer A
      Character*4 Label
#include "WrkSpc.fh"

      If (.not.allocated(AtomicLabels)) Then
         Label='NONE'
      Else
         Label = AtomicLabels(A)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_UnsetAtomicLabels()
C
C     Thomas Bondo Pedersen, March 2011.
C
C     Purpose: unset atomic labels in ldf_atomiclabels
C
      Use ldf_atomiclabels, only: AtomicLabels
      Use stdalloc, only: mma_deallocate
      Implicit None

      If (allocated(AtomicLabels)) Call mma_deallocate(AtomicLabels)

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
      Integer Function LDF_GetLenIn8()
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      LDF_GetLenIn8=LENIN8
      End
