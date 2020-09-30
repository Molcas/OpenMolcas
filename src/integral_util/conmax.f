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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine ConMax(A,nPrim,mPrim,B,nCont,C,mCont)
************************************************************************
*                                                                      *
* Object: to find the largest element in the contraction matrix  for   *
*         each primitive index.                                        *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             July '91                                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 A(nPrim,mPrim), B(nPrim,nCont), C(mPrim,mCont)
#include "print.fh"
#include "real.fh"
*
      iRout = 230
      iPrint = nPrint(iRout)
*
      Do iPrim = 1, nPrim
         Temp= Zero
*        Do iCont = 1, nCont
*           If (Abs(B(iPrim,iCont)).gt.Temp)
*    &         Temp = Abs(B(iPrim,iCont))
*        End Do
         Temp=DDot_(nCont,B(iPrim,1),nPrim,B(iPrim,1),nPrim)
         Do jPrim = 1, mPrim
            A(iPrim,jPrim) = Temp
         End Do
      End Do
*
      Do jPrim = 1, mPrim
*        Temp = Zero
*        Do jCont = 1, mCont
*           If (Abs(C(jPrim,jCont)).gt.Temp)
*    &         Temp = Abs(C(jPrim,jCont))
*        End Do
         Temp=DDot_(mCont,C(jPrim,1),mPrim,C(jPrim,1),mPrim)
         Do iPrim = 1, nPrim
*           A(iPrim,jPrim) = A(iPrim,jPrim)*Temp
            A(iPrim,jPrim) = Sqrt(A(iPrim,jPrim)*Temp)
         End Do
      End Do
*
      Return
      End
