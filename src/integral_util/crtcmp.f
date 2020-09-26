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
* Copyright (C) 1990,2020, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine CrtCmp(Zeta,P,nZeta,A,Axyz,na,HerR,nHer,ABeq)
************************************************************************
*                                                                      *
* Object: to compile the value of the angular part of a basis function *
*         evaluated at a root of the quadrature.                       *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Implicit None
#include "real.fh"
      Integer, Intent(In):: nZeta, nHer, na
      Integer iHer, iCar, ia
      Real*8, Intent(In)::  Zeta(nZeta), P(nZeta,3), A(3), HerR(nHer)
      Real*8, Intent(InOut):: Axyz(nZeta,3,nHer,0:na)
!#define _DEBUG_
#ifdef _DEBUG_
      Character*80 Label
#endif
      Logical ABeq(3)
*
*
      If (na.lt.0) Then
         Call WarningMessage(2,'CrtCmp: na.lt.0')
         Call Abend()
      End If
#ifdef _DEBUG_
      Write (Label,'(A)') ' In CrtCmp: Axyz(in)'
      Call RecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1))
      Call RecPrt(' In CrtCmp: HerR',' ',HerR,1,nHer)
      Call RecPrt(' In CrtCmp: Zeta',' ',Zeta,nZeta,1)
      Call RecPrt(' In CrtCmp: A   ',' ',A   ,1    ,3)
      Call RecPrt(' In CrtCmp: P   ',' ',P   ,nZeta,3)
#endif
      Axyz(:,:,:,0)=One
      If (na.ne.0) then
*
         Do iHer = 1, nHer
            Do iCar = 1, 3
*
               If (ABeq(iCar)) Then
                  Axyz(:,iCar,iHer,1) =
     &                   HerR(iHer)*1/Sqrt(Zeta(:))
               Else
                  Axyz(:,iCar,iHer,1) =
     &                      HerR(iHer)*1/Sqrt(Zeta(:)) +
     &                      P(:,iCar) - A(iCar)
               End If
*
               Do ia = 2, na
                  Axyz(:,iCar,iHer,ia) = Axyz(:,iCar,iHer,1) *
     &                                      Axyz(:,iCar,iHer,ia-1)
               End Do
*
            End Do
         End Do
*
      End If
*
#ifdef _DEBUG_
      Write (Label,'(A)') ' In CrtCmp: Axyz(out) '
      Call RecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1))
#endif
      Return
      End
