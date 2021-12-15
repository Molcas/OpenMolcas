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
      SubRoutine vAssmbl(Rnxyz,Axyz,la,Rxyz,lr,Bxyz,lb,nZeta,HerW,nHer,
     &                  Temp)
************************************************************************
*                                                                      *
* Object: to assemble the cartesian components of the multipole moment *
*         matrix within the framework of the Gauss-Hermite quadrature. *
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
     &       Bxyz(nZeta*3,nHer,0:lb), Temp(nZeta*3,nHer)
      Character*80 Label
*
      iRout = 123
      iPrint = nPrint(iRout)
*     iPrint = 99
      If (iPrint.ge.99) Then
         Call RecPrt(' In vAssmbl:HerW',' ',HerW,1,nHer)
         Call RecPrt(' In vAssmbl:Axyz',' ',Axyz,nZeta*3,nHer*(la+1))
         Call RecPrt(' In vAssmbl:Bxyz',' ',Bxyz,nZeta*3,nHer*(lb+1))
         Call RecPrt(' In vAssmbl:Rxyz',' ',Rxyz,nZeta*3,nHer*(lr+1))
      End If
*
*
      call dcopy_(3*nZeta*(la+1)*(lb+1)*(lr+1),[Zero],0,
     &           Rnxyz,1)
      Do 100 ia = 0, la
         Do 110 ib = 0, lb
            Do 111 iHer = 1, nHer
               Do 112 iZCar = 1, 3*nZeta
                  Temp(iZCar,iHer) =  Axyz(iZCar,iHer,ia)*
     &                                Bxyz(iZCar,iHer,ib)*HerW(iHer)
 112           Continue
 111        Continue
            Do 120 ir = 0, lr
*
*     Generate the cartesian components of the multipole moment
*     matrix as a sum of the value of the integrand, evaluated
*     at a root, times a weight.
*
               Do 30 iHer = 1, nHer
                  Do 10 iZCar = 1, 3*nZeta
                     Rnxyz(iZCar,ia,ib,ir) = Rnxyz(iZCar,ia,ib,ir) +
     &                             Temp(iZCar,iHer)*
     &                             Rxyz(iZCar,iHer,ir)
 10               Continue
 30            Continue
*
               If (iPrint.ge.99) Then
                  Write (Label,'(A,I2,A,I2,A,I2,A)')
     &            ' In vAssmbl: Rnxyz(',ia,',',ib,',',ir,')'
                  Call RecPrt(Label,' ',Rnxyz(1,ia,ib,ir),nZeta,3)
               End If
 120        Continue
 110     Continue
 100  Continue
*
*     Call GetMem(' Exit vAssmbl ','LIST','REAL',iDum,iDum)
      Return
      End
