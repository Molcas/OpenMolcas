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
      Logical Function LDF_CheckIntegrals_JK_2P(Tol,Silent)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: check symmetry of integrals (J_AB | K_CD).
C              Returns .True. if symmetric (debug code).
C
      Implicit None
      Real*8  Tol
      Logical Silent
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair

      Logical S
      Logical LDF_CheckIntegrals_JK_2P_

      Integer AB, CD
      Integer MAB, MCD
      Integer ip_ABCD, l_ABCD
      Integer ip_CDAB, l_CDAB
      Integer iCount

      iCount=0
      Do CD=1,NumberOfAtomPairs
         Do AB=CD,NumberOfAtomPairs
            ! Get auxiliary basis dimensions
            MAB=LDF_nBasAux_Pair(AB)
            MCD=LDF_nBasAux_Pair(CD)
            ! Allocate integral arrays
            l_ABCD=MAB*MCD
            l_CDAB=l_ABCD
            Call GetMem('CIJKABCD','Allo','Real',ip_ABCD,l_ABCD)
            Call GetMem('CIJKCDAB','Allo','Real',ip_CDAB,l_CDAB)
            ! Compute integrals
            Call LDF_ComputeIntegrals_JK_2P(AB,CD,l_ABCD,Work(ip_ABCD))
            Call LDF_ComputeIntegrals_JK_2P(CD,AB,l_CDAB,Work(ip_CDAB))
            ! Check symmetry
            S=LDF_CheckIntegrals_JK_2P_(MAB,MCD,
     &                                  Work(ip_ABCD),Work(ip_CDAB),
     &                                  Tol)
            If (.not.S) iCount=iCount+1
            If (.not.Silent) Then
               If (S) Then
                  Write(6,'(A,I9,1X,I9)')
     &            '(J|K) = (K|J) for atom pairs',AB,CD
               Else
                  Write(6,'(A,I9,1X,I9,A,I9,A)')
     &            '(J|K) != (K|J) for atom pairs',AB,CD,
     &            '(Error',iCount,')'
               End If
            End If
            ! Deallocate
            Call GetMem('CIJKCDAB','Free','Real',ip_CDAB,l_CDAB)
            Call GetMem('CIJKABCD','Free','Real',ip_ABCD,l_ABCD)
         End Do
      End Do

      ! Set return value
      LDF_CheckIntegrals_JK_2P=iCount.eq.0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Logical Function LDF_CheckIntegrals_JK_2P_(N,M,X,Y,Tol)
      Implicit None
      Integer N
      Integer M
      Real*8  X(N,M)
      Real*8  Y(M,N)
      Real*8  Tol

      Integer I, J

      LDF_CheckIntegrals_JK_2P_=.True.
      Do J=1,M
         Do I=1,N
            If (abs(X(I,J)-Y(J,I)).gt.Tol) Then
               LDF_CheckIntegrals_JK_2P_=.False.
               Return
            End If
         End Do
      End Do

      End
