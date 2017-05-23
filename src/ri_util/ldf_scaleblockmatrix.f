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
      Subroutine LDF_ScaleBlockMatrix(ip_Blocks,Factor)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Scale (in place) matrix blocks.
C
      Implicit None
      Integer ip_Blocks
      Real*8  Factor
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Real*8 Tol
      Parameter (Tol=1.0d-15)

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer iAtomPair
      Integer ip, l

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      If (abs(Factor-1.0d0).lt.Tol) Return

      Do iAtomPair=1,NumberOfAtomPairs
         l=LDF_nBas_Atom(AP_Atoms(1,iAtomPair))
     &    *LDF_nBas_Atom(AP_Atoms(2,iAtomPair))
         ip=iWork(ip_Blocks-1+iAtomPair)
         Call dScal_(l,Factor,Work(ip),1)
      End Do

      End
