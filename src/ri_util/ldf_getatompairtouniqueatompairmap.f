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
      Subroutine LDF_GetAtomPairToUniqueAtomPairMap(AP2UAP,nAP)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Get mapping from atom pair to unique atom pair.
C
C TODO/FIXME: only unique atoms are used at the moment.
C
      Implicit None
      Integer nAP
      Integer AP2UAP(nAP)
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

      Integer ip0
      Integer iAtomPair
      Integer iAtom, jAtom
      Integer ip_A2AP, l_A2AP

      Integer  LDF_nAtom, LDF_UniqueAtom
      External LDF_nAtom, LDF_UniqueAtom

      Integer i, j
      Integer A2AP
      Integer AP_Atoms
      A2AP(i)=iWork(ip_A2AP-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

#if defined (_DEBUG_)
      If (nAP.lt.NumberOfAtomPairs) Then
         Call WarningMessage(2,
     &           'LDF_GetAtomPairToUniqueAtomPairMap: dimension error!')
         Call LDF_Quit(1)
      End If
#endif

      l_A2AP=LDF_nAtom()
      Call Getmem('A2AP','Allo','Inte',ip_A2AP,l_A2AP)

      ip0=ip_A2AP-1
      Do iAtomPair=1,NumberOfAtomPairs
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         If (iAtom.eq.jAtom) Then
            iWork(ip0+iAtom)=iAtomPair
         End If
      End Do

      Do iAtomPair=1,NumberOfAtomPairs
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         If (iAtom.eq.jAtom) Then
            AP2UAP(iAtomPair)=A2AP(LDF_UniqueAtom(iAtom))
         Else ! TODO/FIXME: all off-diagonal pairs are considered unique
            AP2UAP(iAtomPair)=iAtomPair
         End If
      End Do

      Call Getmem('A2AP','Free','Inte',ip_A2AP,l_A2AP)

      End
