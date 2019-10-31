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
* Copyright (C) 1990,2015, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine CCrtCmp(Zeta,P,nZeta,A,Axyz,na,HerR,nHer,ABeq,KVector)
************************************************************************
*                                                                      *
* Object: to compile the value of the angular part of a basis function *
*         evaluated at a root of the quadrature.                       *
*                                                                      *
* Called from: PrpInt                                                  *
*                                                                      *
* Calling    :                                                         *
*              RecPrt                                                  *
*              DCopy (ESSL)                                            *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*                                                                      *
*             Roland Lindh, Uppsala Universitet, Uppsala Sweden        *
*             December 2015.                                           *
*             Modification to wave vectors and complex representation. *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Zeta(nZeta), P(nZeta,3), A(3), HerR(nHer), KVector(3)
      Complex*16 Axyz(nZeta,3,nHer,0:na)
      Character*80 Label
      Logical ABeq(3)
*
      iRout = 116
      iPrint = nPrint(iRout)
*
      If (na.lt.0) Then
         Call WarningMessage(2,'CCrtCmp: na.lt.0')
         Call Abend()
      End If
      If (iPrint.ge.99) Then
         Call RecPrt(' In CCrtCmp: HerR',' ',HerR,1,nHer)
         Call RecPrt(' In CCrtCmp: Zeta',' ',Zeta,nZeta,1)
         Call RecPrt(' In CCrtCmp: A   ',' ',A   ,1    ,3)
         Call RecPrt(' In CCrtCmp: P   ',' ',P   ,nZeta,3)
         Call RecPrt(' In CCrtCmp: KVec',' ',KVector,1,3)
      End If
      Do iHer = 1, nHer
         Do iCar = 1, 3
            Do iZeta = 1, nZeta
               Axyz(iZeta,iCar,iHer,0) = DCMPLX(One,Zero)
            End Do
         End Do
      End Do
      If (na.eq.0) Go to 99
*
      Do iHer = 1, nHer
         Do iCar = 1, 3
*
            Do iZeta = 1, nZeta
               Axyz(iZeta,iCar,iHer,1) =
     &                DCMPLX(
     &                HerR(iHer)*1/Sqrt(Zeta(iZeta)) +
     &                P(iZeta,iCar) - A(iCar),
     &                KVector(iCar)/(Two*Zeta(iZeta)))
            End Do
*
            Do ia = 2, na
               Do iZeta = 1, nZeta
                  Axyz(iZeta,iCar,iHer,ia) = Axyz(iZeta,iCar,iHer,1) *
     &                                       Axyz(iZeta,iCar,iHer,ia-1)
               End Do
            End Do
*
         End Do
      End Do
 99   Continue
*
      If (iPrint.ge.99) Then
         Write (Label,'(A)') ' In CCrtCmp: Axyz '
         Call CRecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1),'R')
         Call CRecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1),'I')
      End If
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical_array(ABeq)
      End
