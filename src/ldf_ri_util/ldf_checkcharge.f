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
      Subroutine LDF_CheckCharge(doPrint,PackedD,ip_D,MAE,AB_MAE,
     &                           Q,deltaQ)
C
C     Thomas Bondo Pedersen, February 2011.
C     - based on LDF_CheckAllOverlapIntegrals by T.B. Pedersen.
C
C     Purpose: check charge
C
C              Q = sum_uv D(uv) * S(uv)
C              deltaQ = sum_uv D(uv) * [sum_J C(uv,J)*\int J(r)dr
C                                      -S(uv)]
C
C     Input:
C        doPrint  -- if true, print info
C        PackedD  -- true if LT storage is used for input density
C        ip_D     -- pointer to input density in Work array
C     Output:
C        MAE      -- max abs error
C        AB_MAE   -- atom pair corresponding to MAE
C        Q        -- total charge (as defined above)
C        deltaQ   -- total charge difference (as defined above)
C                    so that Q(LDF)=Q+deltaQ.
C
      Implicit None
      Logical doPrint
      Logical PackedD
      Integer ip_D
      Real*8  MAE
      Integer AB_MAE
      Real*8  Q
      Real*8  deltaQ
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*15 SecNam
      Parameter (SecNam='LDF_CheckCharge')

      Real*8   dDot_, LDF_AtomicDistance
      external ddot_, LDF_AtomicDistance

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair_wLD, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair_wLD, LDF_nBasAux_Pair

      Integer ip_DBlocks, ip_SBlocks
      Integer ip_AuxIntVec
      Integer AB
      Integer A, B
      Integer nAB
      Integer ip_C, l_C, l_Cmax
      Integer ipD, ipS
      Integer ip_dQ, l_dQ
      Integer ip_Stat, l_Stat

      Real*8  QAB, dQ

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      MAE=-9.9d9
      AB_MAE=-1
      Q=0.0d0
      deltaQ=0.0d0
      l_Cmax=0
      Do AB=1,NumberOfAtomPairs
         l_Cmax=max(l_Cmax,LDF_nBas_Atom(AP_Atoms(1,AB))
     &                    *LDF_nBas_Atom(AP_Atoms(2,AB))
     &                    *LDF_nBasAux_Pair_wLD(AB))
      End Do
      If (l_Cmax.lt.1) Return
      If (doPrint) Then
         l_dQ=NumberOfAtomPairs
         Call GetMem('dQ','Allo','Real',ip_dQ,l_dQ)
         Call Cho_dZero(Work(ip_dQ),l_dQ)
      End If
      Call GetMem('Coeff','Allo','Real',ip_C,l_Cmax)
      Call LDF_AllocateBlockMatrix('Den',ip_DBlocks)
      Call LDF_Full2Blocked(Work(ip_D),PackedD,ip_DBlocks)
      Call LDF_ScaleOffdiagonalMatrixBlocks(ip_DBlocks,2.0d0)
      Call LDF_AllocateBlockMatrix('Ovl',ip_SBlocks)
      Call LDF_GetBlockedOverlapMatrix(0,ip_SBlocks)
      Call LDF_AllocateAuxBasVector('Int',ip_AuxIntVec)
      Call LDF_ComputeAuxInt(ip_AuxIntVec)
      If (doPrint) Call Cho_Head('LDF Charge Check','-',80,6)
      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
         l_C=nAB*LDF_nBasAux_Pair_wLD(AB)
         If (l_C.gt.0) Then
            Call LDF_CIO_ReadC_wLD(AB,Work(ip_C),l_C)
            ipD=iWork(ip_DBlocks-1+AB)
            ipS=iWork(ip_SBlocks-1+AB)
            QAB=dDot_(nAB,Work(ipD),1,Work(ipS),1)
            Q=Q+QAB
            Call LDF_ComputeOverlapFromAuxInt(AB,
     &                                1.0d0,l_C,Work(ip_C),ip_AuxIntVec,
     &                                   -1.0d0,nAB,Work(ipS))
            dQ=dDot_(nAB,Work(ipD),1,Work(ipS),1)
            deltaQ=deltaQ+dQ
            If (abs(dQ).gt.MAE) Then
               MAE=abs(dQ)
               AB_MAE=AB
            End If
            If (doPrint) Then
               Work(ip_dQ-1+AB)=dQ
               Write(6,'(/,2X,A,10X,I10,2X,A,2I10,2X,A,1P,D20.10)')
     &         'Atom pair..............',AB,
     &         'Atoms..................',A,B,
     &         'Atomic distance........',LDF_AtomicDistance(A,B)
               Write(6,'(2X,A,5X,I15,2X,A,5X,I15,A)')
     &         'Dimension..............',nAB,
     &         'Auxiliary basis........',
     &         LDF_nBasAux_Pair(AB),' (w/o LinDep)'
               Write(6,'(1P,3(2X,A,D20.10))')
     &         'Charge.................',QAB,
     &         'LDF charge.............',QAB+dQ,
     &         'Error..................',dQ
               Write(6,'(1P,3(2X,A,D20.10))')
     &         'Accumulated charge.....',Q,
     &         'Acccumulated LDF charge',Q+deltaQ,
     &         'Accumulated error......',deltaQ
               Call xFlush(6)
            End If
         End If
      End Do
      Call LDF_DeallocateAuxBasVector('Int',ip_AuxIntVec)
      Call LDF_DeallocateBlockMatrix('Ovl',ip_SBlocks)
      Call LDF_DeallocateBlockMatrix('Den',ip_DBlocks)
      Call GetMem('Coeff','Free','Real',ip_C,l_Cmax)

      If (doPrint) Then
         l_Stat=7
         Call GetMem('Stat','Allo','Real',ip_Stat,l_Stat)
         Call Cho_Head('LDF Charge Error Statistics','-',80,6)
         Call Statistics(Work(ip_dQ),NumberOfAtomPairs,Work(ip_Stat),
     &                   1,2,3,4,5,6,7)
         Write(6,*)
         Write(6,'(1P,3(2X,A,D20.10))')
     &   'Total charge......',Q,
     &   'Total LDF charge..',Q+deltaQ,
     &   'Total LDF error...',deltaQ
         Write(6,'(1P,3(2X,A,D20.10))')
     &   'Average error.....',Work(ip_Stat),
     &   'Standard deviation',Work(ip_Stat+5),
     &   'Abs average error.',Work(ip_Stat+1)
         Write(6,'(1P,2(2X,A,D20.10))')
     &   'Minimum error.....',Work(ip_Stat+2),
     &   'Maximum error.....',Work(ip_Stat+3)
         If (AB_MAE.gt.0) Then
            Write(6,'(/,2X,A,1P,D20.10,1X,A,I10,2X,A,D20.10)')
     &      'Max abs charge error...',MAE,'@AB=',AB_MAE,'Distance=',
     &      LDF_AtomicDistance(AP_Atoms(1,AB_MAE),AP_Atoms(2,AB_MAE))
         End If
         Call xFlush(6)
         Call GetMem('Stat','Free','Real',ip_Stat,l_Stat)
         Call GetMem('dQ','Free','Real',ip_dQ,l_dQ)
      End If

      If (Q.lt.0.0d0 .or. (Q+deltaQ).lt.0.0d0) Then
         Write(6,'(1P,2(2X,A,D20.10))')
     &   'Q=',Q,'Q_LDF=',Q+deltaQ
         Call WarningMessage(2,SecNam//': this is unphysical....')
         Call LDF_Quit(1)
      End If

      End
