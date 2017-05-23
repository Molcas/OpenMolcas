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
      SubRoutine CrtCmp(Zeta,P,nZeta,A,Axyz,na,HerR,nHer,ABeq)
************************************************************************
*                                                                      *
* Object: to compile the value of the angular part of a basis function *
*         evaluated at a root of the quadrature.                       *
*                                                                      *
* Called from: PrpInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy (ESSL)                                            *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Zeta(nZeta), P(nZeta,3), A(3), HerR(nHer),
     &       Axyz(nZeta,3,nHer,0:na)
*define _DEBUG_
#ifdef _DEBUG_
      Character*80 Label
#endif
      Logical ABeq(3)
*
C     Call qEnter('CrtCmp')
*
      If (na.lt.0) Then
         Call WarningMessage(2,'CrtCmp: na.lt.0')
         Call Abend()
      End If
#ifdef _DEBUG_
      Call RecPrt(' In CrtCmp: HerR',' ',HerR,1,nHer)
      Call RecPrt(' In CrtCmp: Zeta',' ',Zeta,nZeta,1)
      Call RecPrt(' In CrtCmp: A   ',' ',A   ,1    ,3)
      Call RecPrt(' In CrtCmp: P   ',' ',P   ,nZeta,3)
#endif
      call dcopy_(nZeta*3*nHer,One,0,Axyz(1,1,1,0),1)
      If (na.eq.0) then
C        Call qExit('CrtCmp')
         Return
      End If
*
      Do 10 iHer = 1, nHer
         Do 20 iCar = 1, 3
*
            If (ABeq(iCar)) Then
               Do 31 iZeta = 1, nZeta
                  Axyz(iZeta,iCar,iHer,1) =
     &                   HerR(iHer)*1/Sqrt(Zeta(iZeta))
*    &                   HerR(iHer)*Zeta(iZeta)**(-Half)
 31            Continue
            Else
               Do 30 iZeta = 1, nZeta
                  Axyz(iZeta,iCar,iHer,1) =
     &                   HerR(iHer)*1/Sqrt(Zeta(iZeta)) +
*    &                   HerR(iHer)*Zeta(iZeta)**(-Half) +
     &                   P(iZeta,iCar) - A(iCar)
 30            Continue
            End If
*
            Do 40 ia = 2, na
               Do 50 iZeta = 1, nZeta
                  Axyz(iZeta,iCar,iHer,ia) = Axyz(iZeta,iCar,iHer,1) *
     &                                       Axyz(iZeta,iCar,iHer,ia-1)
 50            Continue
 40         Continue
*
 20      Continue
 10   Continue
*
#ifdef _DEBUG_
      Write (Label,'(A)') ' In CrtCmp: Axyz '
      Call RecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1))
#endif
C     Call qExit('CrtCmp')
      Return
      End
