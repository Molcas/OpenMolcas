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
      Subroutine LDF_PrintAuxBasInfo(iAtomPair)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Print auxiliary basis info for atom pair iAtomPair.
C
      Implicit None
      Integer iAtomPair
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

#if defined (_DEBUG_)
      Logical  LDF_AtomPairInfoIsSet
      External LDF_AtomPairInfoIsSet
#endif
      Integer  LDF_nBasAux_Pair
      Integer  LDF_nBas_Atom, LDF_nBasAux_Atom
      External LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Atom

      Integer iAtom, jAtom
      Integer n,nProducts

      Integer i, j
      Integer AP_Atoms, AP_1CLinDep, AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

#if defined (_DEBUG_)
      If (.not.LDF_AtomPairInfoIsSet()) Then
         Call WarningMessage(2,
     &                   'LDF_PrintAuxBasInfo: Atom Pair Info not set!')
         Call LDF_Quit(1)
      End If
      If (iAtomPair.lt.1 .or. iAtomPair.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,
     &                  'LDF_PrintAuxBasInfo: iAtomPair out of bounds!')
         Call LDF_Quit(1)
      End If
#endif

      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)
      n=LDF_nBas_Atom(iAtom)
      If (iAtom.eq.jAtom) Then
         nProducts=n*(n+1)/2
      Else
         nProducts=n*LDF_nBas_Atom(jAtom)
      End If

      Write(6,'(/,A,1X,I9)')
     & 'Atom pair....................................',iAtomPair
      Write(6,'(A,2(1X,I9))')
     & 'Atoms........................................',iAtom,jAtom
      Write(6,'(A,1X,I9)')
     & 'Number of product functions to fit...........',nProducts
      Write(6,'(A,1X,I9)')
     & 'Total number of auxiliary basis functions....',
     & LDF_nBasAux_Pair(iAtomPair)
      If (iAtom.eq.jAtom) Then
         Write(6,'(A,1X,I9)')
     &   'Total number of one-center aux functions.....',
     &   LDF_nBasAux_Atom(iAtom)
      Else
         Write(6,'(A,2(1X,I9))')
     &   'Total number of one-center aux functions.....',
     &   LDF_nBasAux_Atom(iAtom)+LDF_nBasAux_Atom(jAtom)
      End If
      Write(6,'(A,1X,I9)')
     & 'Linearly dependent one-center aux functions..',
     & AP_1CLinDep(1,iAtomPair)
      Write(6,'(A,1X,I9)')
     & 'Number of two-center auxiliary functions.....',
     & AP_2CFunctions(1,iAtomPair)

      End
