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
      Subroutine LDF_PrintAtomInfo()
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: Print LDF atom info.
C
      Implicit None
#include "localdf_print.fh"
#include "ldf_atom_info.fh"
#include "WrkSpc.fh"

      Logical  LDF_AtomInfoIsSet, LDF_AtomInfoIsUnset
      External LDF_AtomInfoIsSet, LDF_AtomInfoIsUnset

      Character*4 Label

      Integer  LDF_Coord_Atom
      Integer  LDF_nAtom, LDF_UniqueAtom, LDF_nUniqueAtom
      Integer  LDF_nShell_Atom, LDF_nAuxShell_Atom
      Integer  LDF_lShell_Atom, LDF_lAuxShell_Atom
      External LDF_Coord_Atom
      External LDF_nAtom, LDF_UniqueAtom, LDF_nUniqueAtom
      External LDF_nShell_Atom, LDF_nAuxShell_Atom
      External LDF_lShell_Atom, LDF_lAuxShell_Atom

      Integer i, iAtom, nAtom
      Integer ip, n

      If (.not.LDF_AtomInfoIsSet()) Then
         Write(6,'(A)')
     &   'LDF_PrintAtomInfo called but info not set....'
         If (.not.LDF_AtomInfoIsUnset()) Then
            Write(6,'(A)')
     &      '...but the Unset Flag is not raised !'
            Call WarningMessage(2,
     &           'LDF_PrintAtomInfo: Set and Unset Flags mismatch!')
            Call LDF_Quit(1)
         End If
         Return
      End If

      Call Cho_Head('LDF Atom Info','=',120,6)
      Write(6,'(/,A,I10)')
     & 'Number of atoms................',NumberOfAtoms
      Write(6,'(A,I10)')
     & 'Number of unique atoms.........',LDF_nUniqueAtom()
      If (iPrint.ge.5) Then
         Write(6,'(A,I10,1X,I10)')
     &   'ip_Coord,l_Coord...............',ip_Coord,l_Coord
         Write(6,'(A,I10,1X,I10)')
     &   'ip_A_Unique,l_A_Unique.........',ip_A_Unique,l_A_Unique
         Write(6,'(A,I10,1X,I10)')
     &   'ip_A_Shells,l_A_Shells.........',ip_A_Shells,l_A_Shells
         Write(6,'(A,I10,1X,I10)')
     &   'ip_A_AuxShells,l_A_AuxShells...',ip_A_AuxShells,l_A_AuxShells
      End If

      nAtom=LDF_nAtom()
      If (nAtom.ne.NumberOfAtoms) Then
         Call WarningMessage(2,
     &               'LDF_PrintAtomInfo: Inconsistent number of atoms!')
         Call LDF_Quit(1)
      End If

      Call Cho_Head('Valence Shells','-',120,6)
      Write(6,'(/,A)') '    Atom  #Shells          List of shells'
      Write(6,'(118A1)') ('-',i=1,118)
      Do iAtom=1,nAtom
         n=LDF_nShell_Atom(iAtom)
         If (n.gt.0) Then
            ip=LDF_lShell_Atom(iAtom)
            Call LDF_PrintAtomInfo_(iAtom,n,iWork(ip))
         End If
      End Do
      Write(6,'(118A1)') ('-',i=1,118)

      Call Cho_Head('Auxiliary Shells','-',120,6)
      Write(6,'(/,A)') '    Atom  #Shells          List of shells'
      Write(6,'(118A1)') ('-',i=1,118)
      Do iAtom=1,nAtom
         n=LDF_nAuxShell_Atom(iAtom)
         If (n.gt.0) Then
            ip=LDF_lAuxShell_Atom(iAtom)
            Call LDF_PrintAtomInfo_(iAtom,n,iWork(ip))
         End If
      End Do
      Write(6,'(118A1)') ('-',i=1,118)

      Call LDF_SetAtomicLabels()

      Call Cho_Head('Atomic Coordinates','-',120,6)
      Write(6,'(/,A)')
     & '       Atom            x               y               z'
      Write(6,'(61A1)') ('-',i=1,61)
      Do iAtom=1,nAtom
         Call LDF_GetAtomicLabel(iAtom,Label)
         ip=LDF_Coord_Atom(iAtom)
         Write(6,'(I8,1X,A,3(1X,F15.5))') iAtom,Label,(Work(ip+i),i=0,2)
      End Do
      Write(6,'(61A1)') ('-',i=1,61)

      Call Cho_Head('Unique Atoms','-',120,6)
      Write(6,'(/,A)')
     & '       Atom     Unique'
      Write(6,'(22A1)') ('-',i=1,22)
      Do iAtom=1,nAtom
         Call LDF_GetAtomicLabel(iAtom,Label)
         n=LDF_UniqueAtom(iAtom)
         Write(6,'(I8,1X,A,1X,I8)') iAtom,Label,n
      End Do
      Write(6,'(22A1)') ('-',i=1,22)

      Call LDF_UnsetAtomicLabels()

      Call xFlush(6)

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_PrintAtomInfo_(iAtom,n,List)
      Implicit None
      Integer iAtom, n
      Integer List(n)

      Integer nBatch, iBatch, i1, i2, i, n_this_batch

      If (n.lt.1) Return

      nBatch=(n-1)/11+1
      Do iBatch=1,nBatch
         If (iBatch.eq.nBatch) Then
            n_this_batch=n-11*(nBatch-1)
         Else
            n_this_batch=11
         End If
         i1=11*(iBatch-1)+1
         i2=i1+n_this_batch-1
         If (iBatch.eq.1) Then
            Write(6,'(I8,1X,I8,2X,11(1X,I8))')
     &      iAtom,n,(List(i),i=i1,i2)
         Else
            Write(6,'(19X,11(1X,I8))')
     &      (List(i),i=i1,i2)
         End If
      End Do

      End
