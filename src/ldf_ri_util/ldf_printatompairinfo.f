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
      Subroutine LDF_PrintAtomPairInfo()
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: Print atom pair info.
C
      Implicit None
#include "localdf_print.fh"
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

      Character*2 Unt

      Integer iAtomPair, nAtomPair
      Integer iAtom, jAtom
      Integer ip, iD
      Integer l, n1, n2
      Integer n1max, n2max
      Integer n, M, iAddr, LenOfC

      Real*8 R, Rmax_All, Dmax, Dmax_All, x
      Real*8 DBmax, DBmax_All
      Real*8 LenTot, Byte

      Logical  LDF_AtomPairInfoIsSet, LDF_AtomPairInfoIsUnset
      External LDF_AtomPairInfoIsSet, LDF_AtomPairInfoIsUnset

      Integer  LDF_Coord_Atom
      Integer  LDF_DiskAddressOfC, LDF_nBas_Atom
      Integer  LDF_nAtom, LDF_nBasAux_Pair
      Integer  LDF_nAtomPair, LDF_AtomPair_DiagDim
      Integer  LDF_UniqueAtomPair, LDF_nUniqueAtomPair
      External LDF_Coord_Atom
      External LDF_DiskAddressOfC, LDF_nBas_Atom
      External LDF_nAtom, LDF_nBasAux_Pair
      External LDF_nAtomPair, LDF_AtomPair_DiagDim
      External LDF_UniqueAtomPair, LDF_nUniqueAtomPair

      Real*8   LDF_AtomicDistance, LDF_InteractionRange
      External LDF_AtomicDistance, LDF_InteractionRange

      Integer i,j
      Integer ip_D, ip_DB, n1CLinDep, n2CFunctions
      Integer AP_Atoms
      Logical isUnique
      ip_DB(i)=iWork(ip_AP_DiagBak-1+i)
      ip_D(i)=iWork(ip_AP_Diag-1+i)
      n1CLinDep(i)=iWork(ip_AP_1CLinDep+2*(i-1))
      n2CFunctions(i)=iWork(ip_AP_2CFunctions+2*(i-1))
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      isUnique(i)=iWork(ip_AP_Unique-1+i).eq.i

      If (.not.LDF_AtomPairInfoIsSet()) Then
         Write(6,'(A)')
     &   'LDF_PrintAtomPairInfo called but info not set....'
         If (.not.LDF_AtomPairInfoIsUnset()) Then
            Write(6,'(A)')
     &      '...but the Unset Flag is not raised !'
            Call WarningMessage(2,
     &           'LDF_PrintAtomPairInfo: Set and Unset Flags mismatch!')
            Call LDF_Quit(1)
         End If
         Return
      End If

      Call Cho_Head('LDF Atom Pair Info','=',80,6)
      n=NumberOfAtompairs
      x=dble(LDF_nAtom())
      x=x*(x+1.0d0)/2.0d0
      R=1.0d2*dble(n)/x
      Write(6,'(/,A,I10,A,F7.2,A)')
     & 'Number of atom pairs.................',n,' (',R,'%)'
      n=LDF_nUniqueAtomPair()
      R=1.0d2*dble(n)/x
      Write(6,'(A,I10,A,F7.2,A)')
     & 'Number of unique atom pairs..........',n,' (',R,'%)'
      If (iPrint.ge.5) Then
         Write(6,'(A,I10,1X,I10)')
     &   'ip_AP_Atoms, l_AP_Atoms..............',ip_AP_Atoms,l_AP_Atoms
         Write(6,'(A,I10,1X,I10)')
     &   'ip_AP_Unique, l_AP_Unique............',
     &   ip_AP_Unique,l_AP_Unique
         Write(6,'(A,I10,1X,I10)')
     &   'ip_AP_Diag,l_AP_Diag.................',ip_AP_Diag,l_AP_Diag
         Write(6,'(A,I10,1X,I10)')
     &   'ip_AP_DiagBak,l_AP_DiagBak...........',
     &    ip_AP_DiagBak,l_AP_DiagBak
         Write(6,'(A,I10,1X,I10)')
     &   'ip_AP_1CLinDep, l_AP_1CLinDep........',
     &   ip_AP_1CLinDep,l_AP_1CLinDep
         Write(6,'(A,I10,1X,I10)')
     &   'ip_AP_2CFunctions, l_AP_2CFunctions..',
     &   ip_AP_2CFunctions,l_AP_2CFunctions
         Write(6,'(A,I10,1X,I10)')
     &   'ip_AP_DiskC, l_AP_DiskC..............',ip_AP_DiskC,l_AP_DiskC
      End If

      nAtomPair=LDF_nAtomPair()
      If (nAtomPair.ne.NumberOfAtomPairs) Then
         Call WarningMessage(2,
     &      'LDF_PrintAtomPairInfo: Inconsistent number of atom pairs!')
         Call LDF_Quit(1)
      End If

      Call Cho_Head('Atom Pairs','-',80,6)
      Write(6,'(/,A,A)')
     & 'Atom Pair   Atom 1   Atom 2        Distance    Max Diagonal ',
     & ' Aux Bas 1CLinDep 2CFunctions DiskAddress DimensionC'
      Write(6,'(112A1)') ('-',i=1,112)
      Rmax_All=-9.9d9
      Dmax_All=-9.9d9
      DBmax_All=-9.9d9
      n1max=-9999999
      n2max=-9999999
      LenTot=0.0d0
      Do iAtomPair=1,nAtomPair
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         R=LDF_AtomicDistance(iAtom,jAtom)
         Rmax_All=max(Rmax_All,R)
         l=LDF_AtomPair_DiagDim(iAtomPair)
         ip=ip_D(iAtomPair)
         Dmax=Work(ip)
         Do iD=1,l-1
            Dmax=max(Dmax,Work(ip+iD))
         End Do
         Dmax_All=max(Dmax_All,Dmax)
         ip=ip_DB(iAtomPair)
         DBmax=Work(ip)
         Do iD=1,l-1
            DBmax=max(DBmax,Work(ip+iD))
         End Do
         DBmax_All=max(DBmax_All,DBmax)
         n1=n1CLinDep(iAtomPair)
         n1max=max(n1max,n1)
         n2=n2CFunctions(iAtomPair)
         n2max=max(n2max,n2)
         M=LDF_nBasAux_Pair(iAtomPair)
         iAddr=LDF_DiskAddressOfC(iAtomPair)
         LenOfC=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)*M
         Call Cho_Word2Byte(LenOfC,8,Byte,Unt)
         If (isUnique(iAtomPair)) Then
            LenTot=LenTot+dble(LenOfC)
         End If
         Write(6,'(3(1X,I8),2(1X,D15.6),2(1X,I8),4X,I8,I12,F8.3,1X,A)')
     &   iAtomPair,iAtom,jAtom,R,Dmax,M,n1,n2,iAddr,Byte,Unt
      End Do
      Write(6,'(112A1)') ('-',i=1,112)
      Write(6,'(A,D15.6)')
     & 'Max Distance......................',Rmax_All
      Write(6,'(A,D15.6)')
     & 'Max Diagonal......................',Dmax_All
      Write(6,'(A,D15.6)')
     & 'Max Original Diagonal.............',DBmax_All
      Write(6,'(A,I8)')
     & 'Max 1CLinDep......................',n1max
      Write(6,'(A,I8)')
     & 'Max 2CFunctions...................',n2max
      Call Cho_RWord2Byte(LenTot,Byte,Unt)
      Write(6,'(A,F8.3,1X,A)')
     & 'Total Unique Coefficient Storage...',Byte,Unt

      Call Cho_Head('Unique Atom Pairs','-',80,6)
      Write(6,'(/,A,I10)')
     & 'Number of Unique Atom Pairs:',LDF_nUniqueAtomPair()
      Write(6,'(/,A)')
     & 'Atom Pair   Unique'
      Write(6,'(18A1)') ('-',i=1,18)
      Do iAtomPair=1,nAtomPair
         n=LDF_UniqueAtomPair(iAtomPair)
         Write(6,'(1X,I8,1X,I8)') iAtomPair,n
      End Do
      Write(6,'(18A1)') ('-',i=1,18)

      Call LDF_SetA2AP()
      Call LDF_PrintA2AP()
      Call Cho_Head('Interaction Ranges (Bohr)','-',80,6)
      Write(6,'(/,A,A)')
     & '    Atom          x               y               z',
     & '              Range'
      Write(6,'(72A1)') ('-',i=1,72)
      Do iAtom=1,LDF_nAtom()
         ip=LDF_Coord_Atom(iAtom)
         Write(6,'(I8,3(1X,F15.5),1X,1P,D15.5)')
     &   iAtom,(Work(ip+i),i=0,2),LDF_InteractionRange(iAtom)
      End Do
      Write(6,'(72A1)') ('-',i=1,72)
      Call LDF_UnsetA2AP()

      Call xFlush(6)

      End
