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
      Subroutine LDF_SetAtomInfo(Verbose,irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: Set atom info stored in ldf_atom_info.fh
C
C       NumberOfAtoms = #atoms
C
C       Coord(i,iAtom): x,y,z (i=1,2,3) coordinates of atom iAtom
C
C       A_Unique(iAtom): unique atom corresponding to atom iAtom
C
C       A_Shells(1,iAtom) = #valence shells on atom iAtom
C       A_Shells(2,iAtom) = pointer to list of valence shells on atom
C                           iAtom
C
C       A_AuxShells(1,iAtom)=#auxiliary shells on atom iAtom
C       A_AuxShells(2,iAtom)=pointer to list of auxiliary shells on atom
C                            iAtom
C
C     To unset the information (incl. deallocations), call
C       LDF_UnsetAtomInfo
C
C     To access the data, use the functions described in
C       ldf_atom_info.fh
C
C     Return codes:
C       irc=0: success
C       irc=1: already set
C
      Implicit None
      Logical Verbose
      Integer irc
#include "ldf_atom_info.fh"
#include "localdf_bas.fh"
#include "WrkSpc.fh"

      Character*15 SecNam
      Parameter (SecNam='LDF_SetAtomInfo')

      Logical FirstCall
      Save FirstCall
      Data FirstCall /.True./

      ! Set return code
      irc=0

      ! If already set, print warning message and return
      If (FirstCall) Then
         FirstCall=.False.
      Else
         If (LDF_AtomInfo_Status.eq.LDF_AtomInfo_Set) Then
            If (Verbose) Then
               Call WarningMessage(1,
     &                           SecNam//': LDF Atom Info already set!')
            End If
            irc=1
            Return
         End If
      End If

      ! Get number of atoms
      Call Get_iScalar('Bfn Atoms',NumberOfAtoms)

      ! Allocate and get atomic coordinates
      l_Coord=3*NumberOfAtoms
      Call GetMem('LDF_Coord','Allo','Real',ip_Coord,l_Coord)
      Call Get_dArray('Bfn Coordinates',Work(ip_Coord),l_Coord)

      ! Allocate and get unique atoms
      ! Note: this MUST be done after coordinate fetching...
      l_A_Unique=NumberOfAtoms
      Call GetMem('A_Unique','Allo','Inte',ip_A_Unique,l_A_Unique)
      Call LDF_GetAtomToUniqueAtomMap(iWork(ip_A_Unique),l_A_Unique)

      ! Allocate and set A_Shells and A_AuxShells:
      ! A_Shells(1,iAtom)=#shells on atom iAtom
      ! A_Shells(2,iAtom)=pointer to list of shells on atom iAtom
      ! A_AuxShells(1,iAtom)=#auxiliary shells on atom iAtom
      ! A_AuxShells(2,iAtom)=pointer to list of auxiliary shells on atom
      !                      iAtom
      l_A_Shells=2*NumberOfAtoms
      l_A_AuxShells=l_A_Shells
      Call GetMem('A_Shells','Allo','Inte',ip_A_Shells,l_A_Shells)
      Call GetMem('A_AuxShells','Allo','Inte',ip_A_AuxShells,
     &                                         l_A_AuxShells)
      Call LDF_SetAtomInfo_(nShell_Valence,nShell_Auxiliary,
     &                      NumberOfAtoms,
     &                      iWork(ip_A_Shells),iWork(ip_A_AuxShells))

      ! Mark info set
      LDF_AtomInfo_Status=LDF_AtomInfo_Set

      ! Print
      If (Verbose) Then
         Call LDF_PrintAtomInfo()
      End If

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_SetAtomInfo_(nShell_Valence,nShell_Auxiliary,nAtom,
     &                            A_Shells,A_AuxShells)
      use iSD_data
      Implicit None
      Integer nShell_Valence, nShell_Auxiliary, nAtom
      Integer A_Shells(2,nAtom), A_AuxShells(2,nAtom)
#include "WrkSpc.fh"
#include "nsd.fh"
#include "setup.fh"

      Character*6 Label

      Integer iShell, jShell, iAtom
      Integer ip_Countr, l_Countr
      Integer ip, l

      Integer i
      Integer AtomOfShell
      AtomOfShell(i)=iSD(10,i)

C     Allocate scratch memory.
C     ========================

      l_Countr=nAtom
      Call GetMem('Countr','Allo','Inte',ip_Countr,l_Countr)

C     Set A_Shells
C     ============

      ! Zero A_Shells
      Call iZero(A_Shells,2*nAtom)

      ! Count number of valence shells on each atom
      Do iShell=1,nShell_Valence
         iAtom=AtomOfShell(iShell)
         A_Shells(1,iAtom)=A_Shells(1,iAtom)+1
      End Do

      ! Allocate atomic valence shell lists
      Do iAtom=1,nAtom
         l=A_Shells(1,iAtom)
         If (l.gt.0) Then
            Write(Label,'(A,I4.4)') 'SA',iAtom-1
            Call GetMem(Label,'Allo','Inte',ip,l)
            A_Shells(2,iAtom)=ip
         End If
      End Do

      ! Set atomic valence shell lists
      Call iZero(iWork(ip_Countr),l_Countr)
      Do iShell=1,nShell_Valence
         iAtom=AtomOfShell(iShell)
         iWork(ip_Countr-1+iAtom)=iWork(ip_Countr-1+iAtom)+1
         jShell=iWork(ip_Countr-1+iAtom)
         iWork(A_Shells(2,iAtom)-1+jShell)=iShell
      End Do

C     Set A_AuxShells
C     ===============

      ! Zero A_AuxShells
      Call iZero(A_AuxShells,2*nAtom)

      ! Count number of auxiliary shells on each atom
      Do iShell=nShell_Valence+1,nShell_Valence+nShell_Auxiliary
         iAtom=AtomOfShell(iShell)
         A_AuxShells(1,iAtom)=A_AuxShells(1,iAtom)+1
      End Do

      ! Allocate atomic auxiliary shell lists
      Do iAtom=1,nAtom
         l=A_AuxShells(1,iAtom)
         If (l.gt.0) Then
            Write(Label,'(A,I4.4)') 'AA',iAtom-1
            Call GetMem(Label,'Allo','Inte',ip,l)
            A_AuxShells(2,iAtom)=ip
         End If
      End Do

      ! Set atomic auxiliary shell lists
      Call iZero(iWork(ip_Countr),l_Countr)
      Do iShell=nShell_Valence+1,nShell_Valence+nShell_Auxiliary
         iAtom=AtomOfShell(iShell)
         iWork(ip_Countr-1+iAtom)=iWork(ip_Countr-1+iAtom)+1
         jShell=iWork(ip_Countr-1+iAtom)
         iWork(A_AuxShells(2,iAtom)-1+jShell)=iShell
      End Do

C     Free scratch memory
C     ===================

      Call GetMem('Countr','Free','Inte',ip_Countr,l_Countr)

      End
