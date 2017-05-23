************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
C#define _DEBUG_
C======================================================================C
C This file collects the subroutines and functions used to get info    C
C from the LDF common blocks.                                          C
C======================================================================C
C >>> Logical Function LDF_With2CF()
      Logical Function LDF_With2CF()
      Implicit None
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#if defined (_DEBUG_)
      Logical  LDF_AtomPairInfoIsSet
      External LDF_AtomPairInfoIsSet
#endif
      Integer iAtomPair
      Integer i
      Integer n2CFunctions
      n2CFunctions(i)=iWork(ip_AP_2CFunctions+2*(i-1))
#if defined (_DEBUG_)
      If (.not.LDF_AtomPairInfoIsSet()) Then
         Call WarningMessage(2,'LDF_With2CF: atom pair info unset!')
         Call LDF_Quit(1)
      End If
#endif
      LDF_With2CF=.False.
      iAtomPair=0
      Do While (iAtomPair.lt.NumberOfAtomPairs .and. .not.LDF_With2CF)
         iAtomPair=iAtomPair+1
         LDF_With2CF=n2CFunctions(iAtomPair).gt.0
      End Do
      End
************************************************************************
C >>> Integer Function LDF_UniqueAtomPair(iAtomPair)
      Integer Function LDF_UniqueAtomPair(iAtomPair)
      Implicit None
      Integer iAtomPair
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"
#if defined (_DEBUG_)
      Logical  LDF_AtomPairInfoIsSet
      External LDF_AtomPairInfoIsSet
      Integer  LDF_nAtomPair
      External LDF_nAtomPair
      If (.not.LDF_AtomPairInfoIsSet()) Then
         Call WarningMessage(2,'LDF_UniqueAtomPair: info not set!')
         Call LDF_Quit(1)
      End If
      If (iAtomPair.lt.1 .or. iAtomPair.gt.LDF_nAtomPair()) Then
         Call WarningMessage(2,
     &                   'LDF_UniqueAtomPair: iAtomPair out of bounds!')
         Call LDF_Quit(1)
      End If
#endif
      LDF_UniqueAtomPair=iWork(ip_AP_Unique-1+iAtomPair)
      End
************************************************************************
C >>> Integer Function LDF_UniqueAtom(iAtom)
      Integer Function LDF_UniqueAtom(iAtom)
      Implicit None
      Integer iAtom
#include "ldf_atom_info.fh"
#include "WrkSpc.fh"
#if defined (_DEBUG_)
      Logical  LDF_AtomInfoIsSet
      External LDF_AtomInfoIsSet
      Integer  LDF_nAtom
      External LDF_nAtom
      If (.not.LDF_AtomInfoIsSet()) Then
         Call WarningMessage(2,'LDF_UniqueAtom: info not set!')
         Call LDF_Quit(1)
      End If
      If (iAtom.lt.1 .or. iAtom.gt.LDF_nAtom()) Then
         Call WarningMessage(2,'LDF_UniqueAtom: iAtom out of bounds!')
         Call LDF_Quit(1)
      End If
#endif
      LDF_UniqueAtom=iWork(ip_A_Unique-1+iAtom)
      End
************************************************************************
C >>> Integer Function LDF_nUniqueAtomPair()
      Integer Function LDF_nUniqueAtomPair()
      Implicit None
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"
      Integer iAtomPair
      Integer AP2UAP, i
      AP2UAP(i)=iWork(ip_AP_Unique-1+i)
      LDF_nUniqueAtomPair=0
      Do iAtomPair=1,NumberOfAtomPairs
         If (AP2UAP(iAtomPair).eq.iAtomPair)
     &      LDF_nUniqueAtomPair=LDF_nUniqueAtomPair+1
      End Do
      End
************************************************************************
C >>> Integer Function LDF_nUniqueAtom()
      Integer Function LDF_nUniqueAtom()
      Implicit None
#include "ldf_atom_info.fh"
#include "WrkSpc.fh"
      Integer iAtom
      Integer A2UA, i
      A2UA(i)=iWork(ip_A_Unique-1+i)
      LDF_nUniqueAtom=0
      Do iAtom=1,NumberOfAtoms
         If (A2UA(iAtom).eq.iAtom) LDF_nUniqueAtom=LDF_nUniqueAtom+1
      End Do
      End
************************************************************************
C >>> Integer Function LDF_NumberOfValenceShellsOnAtom(iAtom)
      Integer Function LDF_NumberOfValenceShellsOnAtom(iAtom)
      Implicit None
      Integer iAtom
#include "ldf_atom_info.fh"
#include "WrkSpc.fh"
#if defined (_DEBUG_)
      Logical  LDF_AtomInfoIsSet
      External LDF_AtomInfoIsSet

      If (.not.LDF_AtomInfoIsSet()) Then
         Call SysAbendMsg('LDF_NumberOfValenceShellsOnAtom',
     &                    'LDF Atom Info is not defined!',
     &                    'Unable to proceed')
      End If
      If (iAtom.lt.1 .or. iAtom.gt.NumberOfAtoms) Then
         Call SysAbendMsg('LDF_NumberOfValenceShellsOnAtom',
     &                    'Input variable iAtom out of bounds',
     &                    'iAtom < 1 or iAtom > NumberOfAtoms')
      End If
#endif
      LDF_NumberOfValenceShellsOnAtom=iWork(ip_A_Shells+2*(iAtom-1))
      End
************************************************************************
C >>> Integer Function LDF_nShell_Atom(iAtom)
      Integer Function LDF_nShell_Atom(iAtom)
      Implicit None
      Integer iAtom
      Integer LDF_NumberOfValenceShellsOnAtom
      LDF_nShell_Atom=LDF_NumberOfValenceShellsOnAtom(iAtom)
      End
************************************************************************
C >>> Integer Function LDF_NumberOfAuxiliaryShellsOnAtom(iAtom)
      Integer Function LDF_NumberOfAuxiliaryShellsOnAtom(iAtom)
      Implicit None
      Integer iAtom
#include "ldf_atom_info.fh"
#include "WrkSpc.fh"
#if defined (_DEBUG_)
      Logical  LDF_AtomInfoIsSet
      External LDF_AtomInfoIsSet

      If (.not.LDF_AtomInfoIsSet()) Then
         Call SysAbendMsg('LDF_NumberOfAuxiliaryShellsOnAtom',
     &                    'LDF Atom Info is not defined!',
     &                    'Unable to proceed')
      End If
      If (iAtom.lt.1 .or. iAtom.gt.NumberOfAtoms) Then
         Call SysAbendMsg('LDF_NumberOfAuxiliaryShellsOnAtom',
     &                    'Input variable iAtom out of bounds!',
     &                    'iAtom < 1 or iAtom > NumberOfAtoms')
      End If
#endif
      LDF_NumberOfAuxiliaryShellsOnAtom=iWork(ip_A_AuxShells
     &                                        +2*(iAtom-1))
      End
************************************************************************
C >>> Integer Function LDF_nAuxShell_Atom(iAtom)
      Integer Function LDF_nAuxShell_Atom(iAtom)
      Implicit None
      Integer iAtom
      Integer LDF_NumberOfAuxiliaryShellsOnAtom
      LDF_nAuxShell_Atom=LDF_NumberOfAuxiliaryShellsOnAtom(iAtom)
      End
************************************************************************
C >>> Integer Function LDF_NumberOfAtoms()
      Integer Function LDF_NumberOfAtoms()
      Implicit None
#include "ldf_atom_info.fh"
#if defined (_DEBUG_)
      Logical  LDF_AtomInfoIsSet
      External LDF_AtomInfoIsSet

      If (.not.LDF_AtomInfoIsSet()) Then
         Call SysAbendMsg('LDF_NumberOfAtoms',
     &                    'LDF Atom Info is not defined!',
     &                    'Unable to proceed')
      End If
      If (NumberOfAtoms.lt.1) Then
         Call SysAbendMsg('LDF_NumberOfAtoms',
     &                    'Number of atoms zero or negative!',
     &                    'Unable to proceed')
      End If
#endif
      LDF_NumberOfAtoms=NumberOfAtoms
      End
************************************************************************
C >>> Integer Function LDF_nAtom()
      Integer Function LDF_nAtom()
      Implicit None
      Integer LDF_NumberOfAtoms
      LDF_nAtom=LDF_NumberOfAtoms()
      End
************************************************************************
C >>> Integer Function LDF_nShell()
      Integer Function LDF_nShell()
      Implicit None
#include "localdf_bas.fh"
      LDF_nShell=nShell_Valence
      End
************************************************************************
C >>> Integer Function LDF_nBasSh(iShell)
      Integer Function LDF_nBasSh(iShell)
      Implicit None
      Integer iShell
#include "localdf_bas.fh"
#include "WrkSpc.fh"
#if defined (_DEBUG_)
      If (iShell.lt.1 .or.
     &    iShell.gt.(nShell_Valence+nShell_Auxiliary+1)) Then
         Call WarningMessage(2,'LDF_nBasSh_Atom: iAtom out of bounds')
         Call LDF_Quit(1)
      End If
#endif
      LDF_nBasSh=iWork(ip_nBasSh-1+iShell)
      End
************************************************************************
C >>> Integer Function LDF_nBasSh_Atom(iShell,iAtom)
      Integer Function LDF_nBasSh_Atom(iShell,iAtom)
      Implicit None
      Integer iShell, iAtom
#include "WrkSpc.fh"
      Integer jShell
      Integer  LDF_lShell_Atom, LDF_nBasSh
      External LDF_lShell_Atom, LDF_nBasSh
#if defined (_DEBUG_)
      Logical  LDF_AtomInfoIsSet
      External LDF_AtomInfoIsSet
      Integer  LDF_nAtom, LDF_nShell_Atom
      External LDF_nAtom, LDF_nShell_Atom
      If (.not.LDF_AtomInfoIsSet()) Then
         Call WarningMessage(2,'LDF_nBasSh_Atom: LDF Atom Info not set')
         Call LDF_Quit(1)
      End If
      If (iAtom.lt.1 .or. iAtom.gt.LDF_nAtom()) Then
         Call WarningMessage(2,'LDF_nBasSh_Atom: iAtom out of bounds')
         Call LDF_Quit(1)
      End If
      If (iShell.lt.1 .or. iShell.gt.LDF_nShell_Atom(iAtom)) Then
         Call WarningMessage(2,'LDF_nBasSh_Atom: iShell out of bounds')
         Call LDF_Quit(1)
      End If
#endif
      jShell=iWork(LDF_lShell_Atom(iAtom)-1+iShell)
      LDF_nBasSh_Atom=LDF_nBasSh(jShell)
      End
************************************************************************
C >>> Integer Function LDF_nBasAuxSh_Atom(iAuxShell,iAtom)
      Integer Function LDF_nBasAuxSh_Atom(iAuxShell,iAtom)
      Implicit None
      Integer iAuxShell, iAtom
#include "WrkSpc.fh"
#include "localdf_bas.fh"
      Integer jAuxShell
      Integer  LDF_lAuxShell_Atom
      External LDF_lAuxShell_Atom
#if defined (_DEBUG_)
      Logical  LDF_AtomInfoIsSet
      External LDF_AtomInfoIsSet
      Integer  LDF_nAtom, LDF_nAuxShell_Atom
      External LDF_nAtom, LDF_nAuxShell_Atom
#endif
      Integer i
      Integer nBasSh
      nBasSh(i)=iWork(ip_nBasSh-1+i)
#if defined (_DEBUG_)
      If (.not.LDF_AtomInfoIsSet()) Then
         Call SysWarnMsg('LDF_nBasAuxSh_Atom',
     &                      'LDF Atom Info not set!',
     &                      '')
         Call LDF_Quit(1)
      End If
      If (iAtom.lt.1 .or. iAtom.gt.LDF_nAtom()) Then
         Call SysWarnMsg('LDF_nBasAuxSh_Atom',
     &                      'iAtom out of bounds!',
     &                      'iAtom<1 or iAtom>nAtom')
         Call LDF_Quit(1)
      End If
      If (iAuxShell.lt.1 .or.
     &    iAuxShell.gt.LDF_nAuxShell_Atom(iAtom)) Then
         Call SysWarnMsg('LDF_nBasAuxSh_Atom',
     &                      'iAuxShell out of bounds!',
     &                      'iAuxShell<1 or iAuxShell>nAuxShell_Atom')
         Call LDF_Quit(1)
      End If
#endif
      jAuxShell=iWork(LDF_lAuxShell_Atom(iAtom)+iAuxShell-1)
      LDF_nBasAuxSh_Atom=nBasSh(jAuxShell)
      End
************************************************************************
C >>> Integer Function LDF_nBasAux_Pair(iAtomPair)
      Integer Function LDF_nBasAux_Pair(iAtomPair)
      Implicit None
      Integer iAtomPair
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"
      Integer iAtom, jAtom
      Integer  LDF_nBasAux_Atom
      External LDF_nBasAux_Atom
#if defined (_DEBUG_)
      Integer  nAtom
      Integer  LDF_nAtom
      External LDF_nAtom
      Logical  LDF_AtomPairInfoIsSet
      External LDF_AtomPairInfoIsSet
#endif
      Integer i, j
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
#if defined (_DEBUG_)
      If (.not.LDF_AtomPairInfoIsSet()) Then
         Call WarningMessage(2,
     &                      'LDF_nBasAux_Pair: Atom Pair Info not set!')
         Call LDF_Quit(1)
      End If
      If (iAtomPair.lt.1 .or. iAtomPair.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,
     &                    'LDF_nBasAux_Pair: iAtomPair out of bounds')
         Call LDF_Quit(1)
      End If
#endif
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)
#if defined (_DEBUG_)
      nAtom=LDF_nAtom()
      If (iAtom.lt.1 .or. iAtom.gt.nAtom .or.
     &    jAtom.lt.1 .or. jAtom.gt.nAtom) Then
         Call WarningMessage(2,
     &                    'LDF_nBasAux_Pair: iAtom/jAtom out of bounds')
         Call LDF_Quit(1)
      End If
#endif
      LDF_nBasAux_Pair=LDF_nBasAux_Atom(iAtom)
      If (jAtom.ne.iAtom) Then
         LDF_nBasAux_Pair=LDF_nBasAux_Pair+LDF_nBasAux_Atom(jAtom)
      End If
      LDF_nBasAux_Pair=LDF_nBasAux_Pair-AP_1CLinDep(1,iAtomPair)
      LDF_nBasAux_Pair=LDF_nBasAux_Pair+AP_2CFunctions(1,iAtomPair)
      End
************************************************************************
C >>> Integer Function LDF_nBasAux_Pair_wLD(iAtomPair)
      Integer Function LDF_nBasAux_Pair_wLD(iAtomPair)
      Implicit None
      Integer iAtomPair
      Integer LDF_nBasAux_Pair_WithLinDep
      LDF_nBasAux_Pair_wLD=LDF_nBasAux_Pair_WithLinDep(iAtomPair)
      End
************************************************************************
C >>> Integer Function LDF_nBasAux_Pair_WithLinDep(iAtomPair)
      Integer Function LDF_nBasAux_Pair_WithLinDep(iAtomPair)
      Implicit None
      Integer iAtomPair
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"
      Integer iAtom, jAtom
      Integer  LDF_nBasAux_Atom
      External LDF_nBasAux_Atom
#if defined (_DEBUG_)
      Integer  nAtom
      Integer  LDF_nAtom
      External LDF_nAtom
      Logical  LDF_AtomPairInfoIsSet
      External LDF_AtomPairInfoIsSet
#endif
      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
#if defined (_DEBUG_)
      If (.not.LDF_AtomPairInfoIsSet()) Then
         Call WarningMessage(2,
     &           'LDF_nBasAux_Pair_WithLinDep: Atom Pair Info not set!')
         Call LDF_Quit(1)
      End If
      If (iAtomPair.lt.1 .or. iAtomPair.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,
     &           'LDF_nBasAux_Pair_WithLinDep: iAtomPair out of bounds')
         Call LDF_Quit(1)
      End If
#endif
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)
#if defined (_DEBUG_)
      nAtom=LDF_nAtom()
      If (iAtom.lt.1 .or. iAtom.gt.nAtom .or.
     &    jAtom.lt.1 .or. jAtom.gt.nAtom) Then
         Call WarningMessage(2,
     &         'LDF_nBasAux_Pair_WithLinDep: iAtom/jAtom out of bounds')
         Call LDF_Quit(1)
      End If
#endif
      LDF_nBasAux_Pair_WithLinDep=LDF_nBasAux_Atom(iAtom)
      If (jAtom.ne.iAtom) Then
         LDF_nBasAux_Pair_WithLinDep=LDF_nBasAux_Pair_WithLinDep
     &                              +LDF_nBasAux_Atom(jAtom)
      End If
      LDF_nBasAux_Pair_WithLinDep=LDF_nBasAux_Pair_WithLinDep
     &                           +AP_2CFunctions(1,iAtomPair)
      End
************************************************************************
C >>> Integer Function LDF_nBasAux_Atom(iAtom)
      Integer Function LDF_nBasAux_Atom(iAtom)
      Implicit None
      Integer iAtom
      Integer i
      Integer  LDF_nAuxShell_Atom, LDF_nBasAuxSh_Atom
      External LDF_nAuxShell_Atom, LDF_nBasAuxSh_Atom
#if defined (_DEBUG_)
      Integer  LDF_nAtom
      External LDF_nAtom
      If (iAtom.lt.1 .or. iAtom.gt.LDF_nAtom()) Then
         Call WarningMessage(2,'LDF_nBasAux_Atom: iAtom out of bounds')
         Call LDF_Quit(1)
      End If
#endif
      LDF_nBasAux_Atom=0
      Do i=1,LDF_nAuxShell_Atom(iAtom)
         LDF_nBasAux_Atom=LDF_nBasAux_Atom+LDF_nBasAuxSh_Atom(i,iAtom)
      End Do
      End
************************************************************************
C >>> Integer Function LDF_nBas_Atom(iAtom)
      Integer Function LDF_nBas_Atom(iAtom)
      Implicit None
      Integer iAtom
      Integer i
      Integer  LDF_nShell_Atom, LDF_nBasSh_Atom
      External LDF_nShell_Atom, LDF_nBasSh_Atom
#if defined (_DEBUG_)
      Integer  LDF_nAtom
      External LDF_nAtom
      If (iAtom.lt.1 .or. iAtom.gt.LDF_nAtom()) Then
         Call WarningMessage(2,'LDF_nBas_Atom: iAtom out of bounds')
         Call LDF_Quit(1)
      End If
#endif
      LDF_nBas_Atom=0
      Do i=1,LDF_nShell_Atom(iAtom)
         LDF_nBas_Atom=LDF_nBas_Atom+LDF_nBasSh_Atom(i,iAtom)
      End Do
      End
************************************************************************
C >>> Integer Function LDF_nAuxShell()
      Integer Function LDF_nAuxShell()
      Implicit None
#include "localdf_bas.fh"
      LDF_nAuxShell=nShell_Auxiliary
      End
************************************************************************
C >>> Integer Function LDF_nAtomPair()
      Integer Function LDF_nAtomPair()
      Implicit None
#include "ldf_atom_pair_info.fh"
      LDF_nAtomPair=NumberOfAtomPairs
      End
************************************************************************
C >>> Integer Function LDF_ListOfValenceShellsOnAtom(iAtom)
      Integer Function LDF_ListOfValenceShellsOnAtom(iAtom)
      Implicit None
      Integer iAtom
#include "ldf_atom_info.fh"
#include "WrkSpc.fh"
#if defined (_DEBUG_)
      Logical  LDF_AtomInfoIsSet
      External LDF_AtomInfoIsSet
      If (.not.LDF_AtomInfoIsSet()) Then
         Call SysAbendMsg('LDF_ListOfValenceShellsOnAtom',
     &                    'LDF Atom Info is not defined!',
     &                    'Unable to proceed')
      End If
      If (iAtom.lt.1 .or. iAtom.gt.NumberOfAtoms) Then
         Call SysAbendMsg('LDF_ListOfValenceShellsOnAtom',
     &                    'Input variable iAtom out of bounds',
     &                    'iAtom < 1 or iAtom > NumberOfAtoms')
      End If
#endif
      LDF_ListOfValenceShellsOnAtom=iWork(ip_A_Shells+2*(iAtom-1)+1)
      End
************************************************************************
C >>> Integer Function LDF_lShell_Atom(iAtom)
      Integer Function LDF_lShell_Atom(iAtom)
      Implicit None
      Integer iAtom
      Integer LDF_ListOfValenceShellsOnAtom
      LDF_lShell_Atom=LDF_ListOfValenceShellsOnAtom(iAtom)
      End
************************************************************************
C >>> Integer Function LDF_ListOfAuxiliaryShellsOnAtom(iAtom)
      Integer Function LDF_ListOfAuxiliaryShellsOnAtom(iAtom)
      Implicit None
      Integer iAtom
#include "ldf_atom_info.fh"
#include "WrkSpc.fh"
#if defined (_DEBUG_)
      Logical  LDF_AtomInfoIsSet
      External LDF_AtomInfoIsSet
      If (.not.LDF_AtomInfoIsSet()) Then
         Call SysAbendMsg('LDF_NumberOfAuxiliaryShellsOnAtom',
     &                    'LDF Atom Info is not defined!',
     &                    'Unable to proceed')
      End If
      If (iAtom.lt.1 .or. iAtom.gt.NumberOfAtoms) Then
         Call SysAbendMsg('LDF_ListOfAuxiliaryShellsOnAtom',
     &                    'Input variable iAtom out of bounds!',
     &                    'iAtom < 1 or iAtom > NumberOfAtoms')
      End If
#endif
      LDF_ListOfAuxiliaryShellsOnAtom=iWork(ip_A_AuxShells+2*(iAtom-1)
     &                                      +1)
      End
************************************************************************
C >>> Integer Function LDF_lAuxShell_Atom(iAtom)
      Integer Function LDF_lAuxShell_Atom(iAtom)
      Implicit None
      Integer iAtom
      Integer LDF_ListOfAuxiliaryShellsOnAtom
      LDF_lAuxShell_Atom=LDF_ListOfAuxiliaryShellsOnAtom(iAtom)
      End
************************************************************************
C >>> Integer Function LDF_DiskAddressOfC(iAtomPair)
      Integer Function LDF_DiskAddressOfC(iAtomPair)
C Note: returns -1 if fitting coefficients are not available on disk.
      Implicit None
      Integer iAtomPair
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#if defined (_DEBUG_)
      Logical  LDF_AtomPairInfoIsSet
      External LDF_AtomPairInfoIsSet
      If (.not.LDF_AtomPairInfoIsSet()) Then
         Call WarningMessage(2,'LDF_DiskAddressOfC: info not set')
         Call LDF_Quit(1)
      End If
      If (iAtomPair.lt.1 .or. iAtomPair.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,
     &                    'LDF_DiskAddressOfC: iAtomPair out of bounds')
         Call LDF_Quit(1)
      End If
#endif
      LDF_DiskAddressOfC=iWork(ip_AP_DiskC-1+iAtomPair)
      End
************************************************************************
C >>> Integer Function LDF_Coord_Atom(iAtom)
      Integer Function LDF_Coord_Atom(iAtom)
      Implicit None
      Integer iAtom
#include "ldf_atom_info.fh"
#if defined (_DEBUG_)
      If (l_Coord.lt.1) Then
         Call WarningMessage(2,'LDF_Coord_Atom: Info not set!')
         Call LDF_Quit(1)
      End If
      If (iAtom.lt.1 .or. iAtom.gt.NumberOfAtoms) Then
         Call WarningMessage(2,'LDF_Coord_Atom: iAtom out of bounds!')
         Call LDF_Quit(1)
      End If
#endif
      LDF_Coord_Atom=ip_Coord+3*(iAtom-1)
      End
************************************************************************
C >>> Integer Function LDF_AtomWithCoordinates(R)
      Integer Function LDF_AtomWithCoordinates(R)
      Implicit None
      Real*8 R(3)
#include "ldf_atom_info.fh"
#include "WrkSpc.fh"
      Integer iAtom, ip_i
      Real*8  d
      Integer  LDF_Coord_Atom
      External LDF_Coord_Atom
#if defined (_DEBUG_)
      If (l_Coord.lt.1) Then
         Call WarningMessage(2,
     &             'LDF_AtomWithCoordinates: LDF Atom Info is not set!')
         Call LDF_Quit(1)
      End If
#endif
      LDF_AtomWithCoordinates=0
      Do iAtom=1,NumberOfAtoms
         ip_i=LDF_Coord_Atom(iAtom)
         d=sqrt((Work(ip_i)-R(1))**2
     &         +(Work(ip_i+1)-R(2))**2
     &         +(Work(ip_i+2)-R(3))**2)
         If (d.lt.1.0d-12) Then
            LDF_AtomWithCoordinates=iAtom
            Return
         End If
      End Do
      End
************************************************************************
C >>> Logical Function LDF_AtomPairInfoIsSet()
      Logical Function LDF_AtomPairInfoIsSet()
      Implicit None
#include "ldf_atom_pair_info.fh"
      LDF_AtomPairInfoIsSet=
     &                   LDF_AtomPairInfo_Status.eq.LDF_AtomPairInfo_Set
      End
************************************************************************
C >>> Logical Function LDF_AtomPairInfoIsUnset()
      Logical Function LDF_AtomPairInfoIsUnset()
      Implicit None
#include "ldf_atom_pair_info.fh"
      LDF_AtomPairInfoIsUnset=
     &                 LDF_AtomPairInfo_Status.eq.LDF_AtomPairInfo_Unset
      End
************************************************************************
C >>> Integer Function LDF_AtomPair_DiagDim(iAtomPair)
      Integer Function LDF_AtomPair_DiagDim(iAtomPair)
      Implicit None
      Integer iAtomPair
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"
      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom
#if defined (_DEBUG_)
      Logical  LDF_AtomPairInfoIsSet
      External LDF_AtomPairInfoIsSet
#endif
      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
#if defined (_DEBUG_)
      If (.not.LDF_AtomPairInfoIsSet()) Then
         Call WarningMessage(2,
     &                 'LDF_AtomPair_DiagDim: Atom Pair Info not set!')
         Call LDF_Quit(1)
      End If
      If (iAtomPair.lt.1 .or. iAtomPair.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,
     &                 'LDF_AtomPair_DiagDim: iAtomPair out of bounds!')
         Call LDF_Quit(1)
      End If
#endif
      LDF_AtomPair_DiagDim=LDF_nBas_Atom(AP_Atoms(1,iAtomPair))
     &                    *LDF_nBas_Atom(AP_Atoms(2,iAtomPair))
      End
************************************************************************
C >>> Logical Function LDF_AtomInfoIsSet()
      Logical Function LDF_AtomInfoIsSet()
      Implicit None
#include "ldf_atom_info.fh"
      LDF_AtomInfoIsSet=LDF_AtomInfo_Status.eq.LDF_AtomInfo_Set
      End
************************************************************************
C >>> Logical Function LDF_AtomInfoIsUnset()
      Logical Function LDF_AtomInfoIsUnset()
      Implicit None
#include "ldf_atom_info.fh"
      LDF_AtomInfoIsUnset=LDF_AtomInfo_Status.eq.LDF_AtomInfo_Unset
      End
************************************************************************
C >>> Real*8 Function LDF_AtomicDistance(iAtom,jAtom)
      Real*8 Function LDF_AtomicDistance(iAtom,jAtom)
      Implicit None
      Integer iAtom, jAtom
#include "WrkSpc.fh"
      Integer  LDF_Coord_Atom
      External LDF_Coord_Atom
      Integer ip_i, ip_j
      ip_i=LDF_Coord_Atom(iAtom)
      ip_j=LDF_Coord_Atom(jAtom)
      LDF_AtomicDistance=sqrt((Work(ip_i)-Work(ip_j))**2
     &                       +(Work(ip_i+1)-Work(ip_j+1))**2
     &                       +(Work(ip_i+2)-Work(ip_j+2))**2)
      End
