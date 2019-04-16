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
      SubRoutine Assmbl(Rnxyz,Axyz,la,Rxyz,lr,Bxyz,lb,nZeta,HerW,nHer)
************************************************************************
*                                                                      *
* Object: to assemble the cartesian components of the multipole moment *
*         matrix within the framwork of the Gauss-Hermite quadrature.  *
*                                                                      *
* Called from: PrpInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Rnxyz(nZeta*3,0:la,0:lb,0:lr), HerW(nHer),
     &       Axyz(nZeta*3,nHer,0:la),
     &       Rxyz(nZeta*3,nHer,0:lr),
     &       Bxyz(nZeta*3,nHer,0:lb)
      Character*80 Label
*
      iRout = 123
      iPrint = nPrint(iRout)
*     Call qEnter('Assmbl ')
      If (iPrint.ge.99) Then
         Call RecPrt(' In Assmbl:HerW',' ',HerW,1,nHer)
         Call RecPrt(' In Assmbl:Axyz',' ',Axyz,nZeta*3,nHer*(la+1))
         Call RecPrt(' In Assmbl:Bxyz',' ',Bxyz,nZeta*3,nHer*(lb+1))
         Call RecPrt(' In Assmbl:Rxyz',' ',Rxyz,nZeta*3,nHer*(lr+1))
      End If
*
*
      call dcopy_(3*nZeta*(la+1)*(lb+1)*(lr+1),[Zero],0,
     &           Rnxyz,1)
      Do 100 ia = 0, la
         Do 110 ib = 0, lb
            Do 120 ir = 0, lr
*
*     Generate the cartesian components of the multipole moment
*     matrix as a sum of the value of the integrand, evaluated
*     at a root, times a weight.
*
               Do 30 iHer = 1, nHer
                  Do 10 iZCar = 1, 3*nZeta
c vv. splitted in order to make GNU compiler on Mac more happy.
#ifdef _DARWIN_
                  vv1=Axyz(iZCar,iHer,ia)*Rxyz(iZCar,iHer,ir)
                       vv2=Bxyz(iZCar,iHer,ib)*HerW(iHer)
                  Rnxyz(iZCar,ia,ib,ir)= Rnxyz(iZCar,ia,ib,ir)+vv1*vv2
#else
                     Rnxyz(iZCar,ia,ib,ir) = Rnxyz(iZCar,ia,ib,ir) +
     &                             Axyz(iZCar,iHer,ia)*
     &                             Rxyz(iZCar,iHer,ir)*
     &                             Bxyz(iZCar,iHer,ib)*
     &                             HerW(iHer)
#endif
 10               Continue
 30            Continue
*
               If (iPrint.ge.99) Then
                  Write (Label,'(A,I2,A,I2,A,I2,A)')
     &            ' In Assmbl: Rnxyz(',ia,',',ib,',',ir,')'
                  Call RecPrt(Label,' ',Rnxyz(1,ia,ib,ir),nZeta,3)
               End If
 120        Continue
 110     Continue
 100  Continue
*
*     Call GetMem(' Exit Assmbl ','LIST','REAL',iDum,iDum)
*     Call qExit('Assmbl ')
      Return
      End
