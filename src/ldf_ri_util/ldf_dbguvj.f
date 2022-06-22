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
      Logical Function LDF_DbguvJ(Tol,Silent)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Compare integrals of type (uv|J) computed by two
C              different routines (debug code).
C
      Implicit None
      Real*8  Tol
      Logical Silent
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      real*8 ddot_
      external ddot_

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Integer AB, A, B, nuv, M, iCount
      Integer ip_Int1, l_Int1, ip_Int2, l_Int2
      Real*8  x

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      iCount=0
      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nuv=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
         M=LDF_nBasAux_Pair(AB)
         l_Int1=nuv*M
         If (l_Int1.gt.0) Then
            l_Int2=l_Int1
            Call GetMem('Int1','Allo','Real',ip_Int1,l_Int1)
            Call GetMem('Int2','Allo','Real',ip_Int2,l_Int2)
            Call LDF_SetIndxG(AB)
            Call LDF_ComputeIntegrals_uvJ(AB,l_Int1,Work(ip_Int1))
            Call LDF_UnsetIndxG()
            Call LDF_ComputeIntegrals_uvJ_2P(AB,AB,l_Int2,Work(ip_Int2))
            Call dAXPY_(l_Int1,-1.0d0,Work(ip_Int2),1,Work(ip_Int1),1)
            x=sqrt(dDot_(l_Int1,Work(ip_Int1),1,
     &                         Work(ip_Int1),1))/dble(l_Int1)
            Call GetMem('Int2','Free','Real',ip_Int2,l_Int2)
            Call GetMem('Int1','Free','Real',ip_Int1,l_Int1)
         Else
            x=0.0d0
         End If
         If (.not.Silent) Then
            Write(6,'(A,I9,A,I9,A,1P,D15.6)')
     &      'Atom pair',AB,'   Dimension: ',nuv*M,
     &      '   Normalized diff. norm uvJ-uvJ_2P=',x
         End If
         If (x.gt.Tol) Then
            iCount=iCount+1
         End If
      End Do

      LDF_DbguvJ=iCount.eq.0

      End
