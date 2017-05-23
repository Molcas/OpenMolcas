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
      SubRoutine Kntc_GIAO(Txyz,Rxyz,Wxyz,na,nb,nr,Alpha,Beta,nZeta)
************************************************************************
*                                                                      *
* Object: to assemble the cartesian components of the kinetic energy   *
*         integral from the cartesian components of the overlap        *
*         integral.                                                    *
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
      Real*8 Txyz(nZeta,3,0:na  ,0:nb,  0:1),
     &       Rxyz(nZeta,3,0:na+1,0:nb+1,0:1),
     &       Wxyz(nZeta,3,0:na  ,0:nb  ,  2),
     &       Alpha(nZeta), Beta(nZeta)
      Character*80 Label
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 115
      iPrint = nPrint(iRout)
*     Call qEnter('Kntc_GIAO')
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.99) Then
         Call RecPrt(' In Kntc: Alpha',' ',Alpha,nZeta,1)
         Call RecPrt(' In Kntc: Beta ',' ',Beta ,nZeta,1)
         Do ia = 0, na+1
            Do ib = 0, nb+1
               Write (Label,'(A,I2,A,I2,A)')
     &               ' In Kntc: Rxyz(',ia,',',ib,',0)'
               Call RecPrt(Label,' ',Rxyz(1,1,ia,ib,0),nZeta,3)
               Write (Label,'(A,I2,A,I2,A)')
     &               ' In Kntc: Rxyz(',ia,',',ib,',1)'
               Call RecPrt(Label,' ',Rxyz(1,1,ia,ib,1),nZeta,3)
            End Do
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Do ia = 0, na
         Do ib = 0, nb
            If (ia.eq.0 .and. ib.eq.0) Then
               Do iCar = 1, 3
                  Do iZeta = 1, nZeta
                     Txyz(iZeta,iCar,ia,ib,0) =
     &                   Two * Alpha(iZeta) * Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib+1,0)
                     Txyz(iZeta,iCar,ia,ib,1) =
     &                   Two * Alpha(iZeta) * Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib+1,1)
                     Wxyz(iZeta,iCar,ia,ib,1) =
     &                  -Two * Alpha(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib  ,0)
                     Wxyz(iZeta,iCar,ia,ib,2) =
     &                  -Two *  Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia  ,ib+1,0)
                  End Do
               End Do
            Else If (ia.eq.0) Then
               Do iCar = 1, 3
                  Do iZeta = 1, nZeta
                     Txyz(iZeta,iCar,ia,ib,0) =
     &                   Two * Alpha(iZeta) * Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib+1,0)
     &                - Alpha(iZeta) * ib * Rxyz(iZeta,iCar,ia+1,ib-1,0)
                     Txyz(iZeta,iCar,ia,ib,1) =
     &                   Two * Alpha(iZeta) * Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib+1,1)
     &                - Alpha(iZeta) * ib * Rxyz(iZeta,iCar,ia+1,ib-1,1)
                     Wxyz(iZeta,iCar,ia,ib,1) =
     &                  -Two * Alpha(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib  ,0)
                     Wxyz(iZeta,iCar,ia,ib,2) =
     &                  -Two *  Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia  ,ib+1,0)
     &                  +              ib * Rxyz(iZeta,iCar,ia  ,ib-1,0)
                  End Do
               End Do
            Else If (ib.eq.0) Then
               Do iCar = 1, 3
                  Do iZeta = 1, nZeta
                     Txyz(iZeta,iCar,ia,ib,0) =
     &                   Two * Alpha(iZeta) * Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib+1,0)
     &                - Beta(iZeta)  * ia * Rxyz(iZeta,iCar,ia-1,ib+1,0)
                     Txyz(iZeta,iCar,ia,ib,1) =
     &                   Two * Alpha(iZeta) * Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib+1,1)
     &                - Beta(iZeta)  * ia * Rxyz(iZeta,iCar,ia-1,ib+1,1)
                     Wxyz(iZeta,iCar,ia,ib,1) =
     &                  -Two * Alpha(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib  ,0)
     &                  +              ia * Rxyz(iZeta,iCar,ia-1,ib  ,0)
                     Wxyz(iZeta,iCar,ia,ib,2) =
     &                  -Two *  Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia  ,ib+1,0)
                  End Do
               End Do
            Else
               Do iCar = 1, 3
                  Do iZeta = 1, nZeta
                     Txyz(iZeta,iCar,ia,ib,0) =
     &                     Half * ia * ib * Rxyz(iZeta,iCar,ia-1,ib-1,0)
     &                - Beta(iZeta)  * ia * Rxyz(iZeta,iCar,ia-1,ib+1,0)
     &                - Alpha(iZeta) * ib * Rxyz(iZeta,iCar,ia+1,ib-1,0)
     &                + Two * Alpha(iZeta) * Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib+1,0)
                     Txyz(iZeta,iCar,ia,ib,1) =
     &                     Half * ia * ib * Rxyz(iZeta,iCar,ia-1,ib-1,1)
     &                - Beta(iZeta)  * ia * Rxyz(iZeta,iCar,ia-1,ib+1,1)
     &                - Alpha(iZeta) * ib * Rxyz(iZeta,iCar,ia+1,ib-1,1)
     &                + Two * Alpha(iZeta) * Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib+1,1)
                     Wxyz(iZeta,iCar,ia,ib,1) =
     &                  -Two * Alpha(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia+1,ib  ,0)
     &                  +              ia * Rxyz(iZeta,iCar,ia-1,ib  ,0)
                     Wxyz(iZeta,iCar,ia,ib,2) =
     &                  -Two *  Beta(iZeta) *
     &                                      Rxyz(iZeta,iCar,ia  ,ib+1,0)
     &                  +              ib * Rxyz(iZeta,iCar,ia  ,ib-1,0)
                  End Do
               End Do
            End If
*                                                                      *
************************************************************************
*                                                                      *
            If (iPrint.ge.99) Then
               Write (Label,'(A,I2,A,I2,A)') ' In Kntc: Txyz(',ia,',',
     &                ib,',0)'
               Call RecPrt(Label,' ',Txyz(1,1,ia,ib,0),nZeta,3)
               Write (Label,'(A,I2,A,I2,A)') ' In Kntc: Txyz(',ia,',',
     &                ib,',1)'
               Call RecPrt(Label,' ',Txyz(1,1,ia,ib,1),nZeta,3)
               Write (Label,'(A,I2,A,I2,A)') ' In Kntc: Wxyz(',ia,',',
     &                ib,',1)'
               Call RecPrt(Label,' ',Wxyz(1,1,ia,ib,1),nZeta,3)
               Write (Label,'(A,I2,A,I2,A)') ' In Kntc: Wxyz(',ia,',',
     &                ib,',2)'
               Call RecPrt(Label,' ',Wxyz(1,1,ia,ib,2),nZeta,3)
            End If
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Call qExit('Kntc_GIAO')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nr)
      End
