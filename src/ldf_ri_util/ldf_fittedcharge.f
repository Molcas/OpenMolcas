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
      Real*8 Function LDF_FittedCharge(PackedD,ip_D)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: compute charge from density using only atom pair blocks
C              as defined by LDF and the fitted overlap matrix,
C              S(uv)=sum_J C(uv,J)*n(J),  n(J)=\int J(r)dr.
C
      Implicit None
      Logical PackedD
      Integer ip_D
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      real*8 ddot_
      external ddot_

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair_wLD
      External LDF_nBas_Atom, LDF_nBasAux_Pair_wLD

      Integer ip_DBlocks
      Integer ip_AuxIntVec
      Integer AB
      Integer A, B
      Integer nAB, nABmax
      Integer ipD
      Integer ip_S, l_S
      Integer ip_C, l_C, l_Cmax

      Real*8 Factor

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      LDF_FittedCharge=0.0d0
      l_Cmax=0
      nABmax=0
      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
         l_Cmax=max(l_Cmax,nAB*LDF_nBasAux_Pair_wLD(AB))
         nABmax=max(nABmax,nAB)
      End Do
      If (l_Cmax.lt.1) Return
      l_S=nABmax
      Call GetMem('FQS','Allo','Real',ip_S,l_S)
      Call GetMem('FQC','Allo','Real',ip_C,l_Cmax)
      Call LDF_AllocateBlockMatrix('Den',ip_DBlocks)
      Call LDF_Full2Blocked(Work(ip_D),PackedD,ip_DBlocks)
      Call LDF_AllocateAuxBasVector('Int',ip_AuxIntVec)
      Call LDF_ComputeAuxInt(ip_AuxIntVec)
      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
         l_C=nAB*LDF_nBasAux_Pair_wLD(AB)
         If (l_C.gt.0) Then
            Call LDF_CIO_ReadC_wLD(AB,Work(ip_C),l_C)
            Call LDF_ComputeOverlapFromAuxInt(AB,
     &                                1.0d0,l_C,Work(ip_C),ip_AuxIntVec,
     &                                             0.0d0,nAB,Work(ip_S))
            ipD=iWork(ip_DBlocks-1+AB)
            Factor=1.0d0+dble(min(abs(A-B),1))
            LDF_FittedCharge=LDF_FittedCharge
     &                      +Factor*dDot_(nAB,Work(ipD),1,Work(ip_S),1)
         End If
      End Do
      Call LDF_DeallocateAuxBasVector('Int',ip_AuxIntVec)
      Call LDF_DeallocateBlockMatrix('Den',ip_DBlocks)
      Call GetMem('FQS','Free','Real',ip_S,l_S)
      Call GetMem('FQC','Free','Real',ip_C,l_Cmax)

      End
