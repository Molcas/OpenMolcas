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
      Subroutine LDF_ZeroAuxBasVector(ip_V)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Zero aux bas vector.
C
      Implicit None
      Integer ip_V
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nAtom, LDF_nBasAux_Atom
      External LDF_nAtom, LDF_nBasAux_Atom

      Integer iAtom, nAtom, iAtomPair
      Integer ip0, ip, l

      Integer i, j
      Integer AP_2CFunctions
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      nAtom=LDF_nAtom()
      ip0=ip_V-1
      Do iAtom=1,nAtom
         l=LDF_nBasAux_Atom(iAtom)
         ip=iWork(ip0+iAtom)
         Call Cho_dZero(Work(ip),l)
      End Do
      ip0=ip0+nAtom
      Do iAtomPair=1,NumberOfAtomPairs
         l=AP_2CFunctions(1,iAtomPair)
         ip=iWork(ip0+iAtomPair)
         Call Cho_dZero(Work(ip),l)
      End Do

      End
