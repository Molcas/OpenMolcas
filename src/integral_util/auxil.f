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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SUBROUTINE Auxil(T,nT,Fm,mHigh)
************************************************************************
*                                                                      *
*     Object: to compute the auxiliary functions in quadruple precision*
*             for a number of arguments.                               *
*                                                                      *
* Called from: RtsWgh                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              Fm                                                      *
*              RecPrt                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      REAL*8 Fm(nT,0:mHigh), T(nT)
*
      iRout = 53
      iPrint = nPrint(iRout)
      Call qEnter('Auxil')
*
      Call HighFm(Fm(1,mHigh),T,mHigh,nT)
*
*     Now use recusion formula for Fm, 0<=m<mHigh
*
      Do 30  i = 1, nT
         Ti=T(i)
         Do 31 m = mHigh-1, 0, -1
            Fm(i,m) = (Two*Ti*Fm(i,m+1)+Exp(-Ti))/DBLE(2*m+1)
 31      Continue
 30   Continue
*     Call RecPrt(' Fm',' ',Fm,nT,mHigh+1)
*     Call GetMem('Auxil','CHECK','REAL',iDum,iDum)
*
      Call qExit('Auxil')
      Return
      End
