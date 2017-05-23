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
      Subroutine LDF_AuxBasVectorNorm(ip_V,ip_Norm)
C
C     Thomas Bondo Pedersen, December 2010.
C
C     Purpose: compute Frobenius norm of aux bas vector.
C
      Implicit None
      Integer ip_V
      Integer ip_Norm
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nAtom, LDF_nBasAux_Atom
      External LDF_nAtom, LDF_nBasAux_Atom

      real*8 ddot_
      external ddot_

      Integer iAtom, nAtom, iAtomPair
      Integer ip0, ipN, ip, l

      Integer i, j
      Integer AP_2CFunctions
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      nAtom=LDF_nAtom()
      ip0=ip_V-1
      ipN=ip_Norm-1
      Do iAtom=1,nAtom
         l=LDF_nBasAux_Atom(iAtom)
         ip=iWork(ip0+iAtom)
         Work(ipN+iAtom)=sqrt(dDot_(l,Work(ip),1,Work(ip),1))
      End Do
      ip0=ip0+nAtom
      ipN=ipN+nAtom
      Do iAtomPair=1,NumberOfAtomPairs
         l=AP_2CFunctions(1,iAtomPair)
         ip=iWork(ip0+iAtomPair)
         Work(ipN+iAtomPair)=sqrt(dDot_(l,Work(ip),1,Work(ip),1))
      End Do

      End
