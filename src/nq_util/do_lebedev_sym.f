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
* Copyright (C) 2019, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine Do_Lebedev_Sym(L_Eff,mPt,ipR)
      Implicit None
#include "WrkSpc.fh"
      Integer, Intent(In) :: L_Eff
      Integer, Intent(Out) :: mPt, ipR
      Integer :: mPt_, ipR_, i, j, ii, jj
      Real*8, Parameter :: Thr = 1.0D-16
      Call Do_Lebedev(L_Eff,mPt_,ipR_)
      mPt=0
      outer: Do i=1,mPt_
        ii = ipR_+(i-1)*4
        Do j=1,i-1
          jj = ipR_+(j-1)*4
          If (All(Abs(Work(jj:jj+2)+Work(ii:ii+2)).lt.Thr)) Then
            Work(ii+3)=0.0D0
            Cycle outer
          End If
        End Do
        mPt=mPt+1
      End Do outer
      Call GetMem('AngRW','Allo','Real',ipR,4*mPt)
      j=1
      Do i=1,mPt_
        ii = ipR_+(i-1)*4
        If (Work(ii+3).ne.0.0D0) Then
          Call DCopy_(4,Work(ii),1,Work(ipR+(j-1)*4),1)
          j=j+1
        End if
      End Do
      Call GetMem('AngRW','Free','Real',ipR_,4*mPt_)
      End Subroutine Do_Lebedev_Sym
