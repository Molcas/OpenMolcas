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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Real*8 Function LDF_Charge(PackedD,ip_D)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: compute charge from density using only atom pair blocks
C              as defined by LDF.
C
      Implicit None
      Logical PackedD
      Integer ip_D
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      real*8 ddot_
      external ddot_

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer ip_DBlocks
      Integer ip_SBlocks
      Integer AB
      Integer A, B
      Integer nAB
      Integer ipD
      Integer ipS

      Real*8 Factor

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      LDF_Charge=0.0d0
      Call LDF_AllocateBlockMatrix('Den',ip_DBlocks)
      Call LDF_Full2Blocked(Work(ip_D),PackedD,ip_DBlocks)
      Call LDF_AllocateBlockMatrix('Ovl',ip_SBlocks)
      Call LDF_GetBlockedOverlapMatrix(0,ip_SBlocks)
      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
         ipD=iWork(ip_DBlocks-1+AB)
         ipS=iWork(ip_SBlocks-1+AB)
         Factor=1.0d0+dble(min(abs(A-B),1))
         LDF_Charge=LDF_Charge+Factor*dDot_(nAB,Work(ipD),1,Work(ipS),1)
      End Do
      Call LDF_DeallocateBlockMatrix('Ovl',ip_SBlocks)
      Call LDF_DeallocateBlockMatrix('Den',ip_DBlocks)

      End
