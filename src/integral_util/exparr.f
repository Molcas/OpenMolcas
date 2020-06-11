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
      SubRoutine ExpArr(Array,Ind,nArray,lArray)
************************************************************************
*                                                                      *
* Object: to expand arrays according to an index array.                *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             Augusti '91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Array(lArray,nArray)
      Integer Ind(nArray)
*
      iRout = 235
      iPrint = nPrint(iRout)
      Call qEnter('ExpArr')
*
      Do 100 iArray = nArray, 1, -1
         jArray = Ind(iArray)
         If (jArray.le.0) Then
*           Set column iArray to zero
            call dcopy_(lArray,[Zero],0,Array(1,iArray),1)
         Else If (jArray.lt.iArray) Then
*           Copy row jArray to position iArray
            call dcopy_(lArray,Array(1,jArray),1,Array(1,iArray),1)
         End If
 100  Continue
      Call qExit('ExpArr')
      Return
      End
