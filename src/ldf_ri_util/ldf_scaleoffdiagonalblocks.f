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
      Subroutine LDF_ScaleOffdiagonalMatrixBlocks(ip_Blocks,Factor)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Scale (in place) off-diagonal blocks.
C
      Implicit None
      Integer ip_Blocks
      Real*8  Factor
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer iAtomPair
      Integer iAtom, jAtom
      Integer ip, l

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      Do iAtomPair=1,NumberOfAtomPairs
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         If (iAtom.ne.jAtom) Then
            l=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
            ip=iWork(ip_Blocks-1+iAtomPair)
            Call dScal_(l,Factor,Work(ip),1)
         End If
      End Do

      End
