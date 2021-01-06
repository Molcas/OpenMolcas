************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_Distrib_Vec(Jin,Jfi,iDV,nV)
C
C     Unless you know exactly what you are doing,
C     do NOT call this routine directly; use Cho_P_Distrib_Vec instead!
C
      Use Para_Info, Only: MyRank, nProcs
      Implicit None
      Integer  Jin, Jfi, nV
      Integer  iDV(*)

      Integer i, iNode

      nV=0
      Do i=Jin,Jfi
         iNode=MOD(i-1,nProcs)
         If (myRank .eq. iNode) Then
            nV = nV + 1
            iDV(nV) = i
         End If
      End Do

      End
