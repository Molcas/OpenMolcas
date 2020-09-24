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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine HrrCtl(Arr1,nArr1,Arr2,nArr2,
     &                  la,lb,lc,ld,nabMax,ncdMax,nTR,
     &                  A,B,C,D,IfGrad)
************************************************************************
*                                                                      *
* Object: to act as a shell towards the HRR subroutines.               *
*                                                                      *
* Called from: Rysg1                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              Hrr2Da                                                  *
*              Hrr2Db                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Arr1(nTR,3*nArr1), Arr2(nTR,3*nArr2),
     &       A(3), B(3), C(3), D(3)
      Logical IfGrad(3,4)
*
      iRout = 233
      iPrint = nPrint(iRout)
*
      Call Hrr2Da(Arr1,nTR,nabMax,ncdMax,Arr2,A,B,la,lb,lc,ld,IfGrad)
*
      Call Hrr2Db(Arr2,nTR,       ncdMax,Arr1,C,D,la,lb,lc,ld,IfGrad)
*
      Return
      End
