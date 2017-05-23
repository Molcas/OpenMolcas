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
      Subroutine LDF_UpdateDiagonalFromC(Constraint,AB,l_C,C,irc)
      Implicit None
      Integer Constraint
      Integer AB
      Integer l_C
      Real*8  C(l_C)
      Integer irc
#include "localdf.fh"
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*23 SecNam
      Parameter (SecNam='LDF_UpdateDiagonalFromC')

      Integer  LDF_AtomPair_DiagDim, LDF_nBasAux_Pair
      External LDF_AtomPair_DiagDim, LDF_nBasAux_Pair

      Integer nAB
      Integer M
      Integer ip_G, l_G
      Integer ip_X, l_X
      Integer i, J, iJ
      Integer ipD
      Integer ipX0, ipX

      nAB=LDF_AtomPair_DiagDim(AB)
      M=LDF_nBasAux_Pair(AB)
      If (nAB.lt.1 .or. M.lt.1) Return
      If (l_C.lt.(nAB*M)) Then
         Call WarningMessage(2,SecNam//': insufficient array dimension')
         Call LDF_Quit(1)
      End If

      Call LDF_SetIndxG(AB)
      l_G=M**2
      Call GetMem('UDFCG','Allo','Real',ip_G,l_G)
      Call LDF_ComputeGMat(AB,M,Work(ip_G))
      l_X=nAB*M
      Call GetMem('UDFCX','Allo','Real',ip_X,l_X)
      Call LDF_ComputeIntegrals_uvJ(AB,l_X,Work(ip_X))
      Call dGeMM_('N','N',nAB,M,M,
     &            -1.0d0,C,nAB,Work(ip_G),M,
     &             2.0d0,Work(ip_X),nAB)
      ipD=iWork(ip_AP_Diag-1+AB)-1
      ipX0=ip_X-1
      Do J=0,M-1
         iJ=nAB*J
         ipX=ipX0+iJ
         Do i=1,nAB
            Work(ipD+i)=Work(ipD+i)-C(iJ+i)*Work(ipX+i)
         End Do
      End Do
      Call GetMem('UDFCX','Free','Real',ip_X,l_X)
      Call GetMem('UDFCG','Free','Real',ip_G,l_G)
      Call LDF_UnsetIndxG()

      irc=0
      Do i=1,nAB
         If (Work(ipD+i).lt.TooNegative) Then
            irc=irc+1
         End If
      End Do

c Avoid unused argument warnings
      If (.False.) Call Unused_integer(Constraint)
      End
