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
      SubRoutine vCrtCmp(Zeta12,P,nZeta,A,Axyz,na,HerR,nHer,ABeq)
************************************************************************
*                                                                      *
* Object: to compile the value of the angular part of a basis function *
*         evaluated at a root of the quadrature.                       *
*                                                                      *
* Called from: PrpInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy                                                   *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Zeta12(nZeta), P(nZeta,3), A(3), HerR(nHer),
     &       Axyz(nZeta,3,nHer,0:na)
      Logical ABeq(3)
*
      iRout = 116
      iPrint = nPrint(iRout)
*     iPrint = 99
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In vCrtCmp: HerR',' ',HerR,1,nHer)
         Call RecPrt(' In vCrtCmp: Zeta',' ',Zeta12,nZeta,1)
         Call RecPrt(' In vCrtCmp: A   ',' ',A   ,1    ,3)
         Call RecPrt(' In vCrtCmp: P   ',' ',P   ,nZeta,3)
      End If
      call dcopy_(nZeta*3*nHer,[One],0,Axyz(1,1,1,0),1)
      If (na.eq.0) Go To 999
*
      Do 10 iHer = 1, nHer
         Do 20 iCar = 1, 3
*
            If (ABeq(iCar)) Then
               Do 31 iZeta = 1, nZeta
                  Axyz(iZeta,iCar,iHer,1) =
     &                   HerR(iHer)*Zeta12(iZeta)
 31            Continue
            Else
               Do 30 iZeta = 1, nZeta
                  Axyz(iZeta,iCar,iHer,1) =
     &                   HerR(iHer)*Zeta12(iZeta) +
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
 999  Continue
*     If (iPrint.ge.99) Then
*        Do 100 ia = 0, na
*           Do 110 iHer = 1, nHer
*              Write (Label,'(A,I2,A,I2,A)') ' In vCrtCmp: Axyz (iHer=',
*    &               iHer,',ia=',ia,')'
*              Call RecPrt(Label,' ',Axyz(1,1,iHer,ia),nZeta,3)
*110        Continue
*100     Continue
*     End If
*     Call GetMem(' Exit vCrtCmp','CHECK','REAL',iDum,iDum)
      Return
      End
