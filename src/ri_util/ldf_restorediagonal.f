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
      Subroutine LDF_RestoreDiagonal(iAtomPair)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: restore original diagonal for atom pair iAtomPair.
C
      Implicit None
      Integer iAtomPair
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_AtomPair_DiagDim
      External LDF_AtomPair_DiagDim

      Integer l, ip1, ip2

      l=LDF_AtomPair_DiagDim(iAtomPair)
      ip1=iWork(ip_AP_DiagBak-1+iAtomPair)
      ip2=iWork(ip_AP_Diag-1+iAtomPair)
      Call dCopy_(l,Work(ip1),1,Work(ip2),1)

      End
