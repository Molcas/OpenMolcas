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
* Copyright (C) 1990,1998, Roland Lindh                                *
*               1998, Samuel Mikes                                     *
************************************************************************
      SubRoutine MVe(rV2Int,rV4Int,Sxyz,na,nb,Alpha,Beta,nZeta)
************************************************************************
*                                                                      *
* Object: to compute intermediate integrals for the evaluation of the  *
*          mass velocity integrals.                                    *
*                                                                      *
* Called from: MVeInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*                                                                      *
*             Correct out of bound reference, February '98, Samuel     *
*             Mikes and Roland Lindh.                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 rV2Int(nZeta,3,0:na,0:nb,2),
     &       rV4Int(nZeta,3,0:na,0:nb),
     &       Sxyz(nZeta,3,0:na+2,0:nb+2),
     &       Alpha(nZeta), Beta(nZeta)
      Character*80 Label
*
      iRout = 192
      iPrint = nPrint(iRout)
      Call qEnter('MVe')
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In MVe: Alpha',' ',Alpha,nZeta,1)
         Call RecPrt(' In MVe: Beta ',' ',Beta ,nZeta,1)
         Do ib = 0, nb+2
            Do ia = 0, na+2
               Write (Label,'(A,I2,A,I2,A)')
     &               ' In MVe: Sxyz(',ia,',',ib,')'
               Call RecPrt(Label,' ',Sxyz(1,1,ia,ib),nZeta,3)
            End Do
         End Do
*
      End If
      Do ib = 0, nb
         Do ia = 0, na
               Do iCar = 1, 3
                  Do iZeta = 1, nZeta
*
                     rV2Int(iZeta,iCar,ia,ib,1) =
     &                  Four*Alpha(iZeta)**2*Sxyz(iZeta,iCar,ia+2,ib)
     &                 -Two*Alpha(iZeta)*(Two*Dble(ia)+One)*
     &                                       Sxyz(iZeta,iCar,ia,ib)
                     If (ia.ge.2) Then
                        rV2Int(iZeta,iCar,ia,ib,1) =
     &                     rV2Int(iZeta,iCar,ia,ib,1)
     &                     +Dble(ia*(ia-1))*  Sxyz(iZeta,iCar,ia-2,ib)
                     End If
*
                     rV2Int(iZeta,iCar,ia,ib,2) =
     &                  Four* Beta(iZeta)**2*Sxyz(iZeta,iCar,ia,ib+2)
     &                 -Two* Beta(iZeta)*(Two*Dble(ib)+One)*
     &                                       Sxyz(iZeta,iCar,ia,ib)
                     If (ib.ge.2) Then
                        rV2Int(iZeta,iCar,ia,ib,2) =
     &                     rV2Int(iZeta,iCar,ia,ib,2)
     &                     +Dble(ib*(ib-1))*  Sxyz(iZeta,iCar,ia,ib-2)
                     End If
*
                     rV4Int(iZeta,iCar,ia,ib) =
     &                  Four*Alpha(iZeta)**2*
     &                  Four* Beta(iZeta)**2*Sxyz(iZeta,iCar,ia+2,ib+2)
     &                 -Four*Alpha(iZeta)**2*
     &                  Two* Beta(iZeta)*(Two*Dble(ib)+One)*
     &                                       Sxyz(iZeta,iCar,ia+2,ib)
     &                 -Four* Beta(iZeta)**2*
     &                  Two*Alpha(iZeta)*(Two*Dble(ia)+One)*
     &                                        Sxyz(iZeta,iCar,ia,ib+2)
     &                 +Two*Alpha(iZeta)*(Two*Dble(ia)+One)*
     &                  Two* Beta(iZeta)*(Two*Dble(ib)+One)*
     &                                       Sxyz(iZeta,iCar,ia,ib)
                     If (ia.ge.2) Then
                        rV4Int(iZeta,iCar,ia,ib) =
     &                     rV4Int(iZeta,iCar,ia,ib) + Dble(ia*(ia-1))* (
     &                     Four* Beta(iZeta)**2 *
     &                                       Sxyz(iZeta,iCar,ia-2,ib+2)
     &                    -Two* Beta(iZeta)*(Two*Dble(ib)+One)*
     &                                       Sxyz(iZeta,iCar,ia-2,ib) )
                     End If
                     If (ib.ge.2) Then
                        rV4Int(iZeta,iCar,ia,ib) =
     &                     rV4Int(iZeta,iCar,ia,ib) + Dble(ib*(ib-1))* (
     &                     Four*Alpha(iZeta)**2 *
     &                                       Sxyz(iZeta,iCar,ia+2,ib-2)
     &                    -Two*Alpha(iZeta)*(Two*Dble(ia)+One)*
     &                                       Sxyz(iZeta,iCar,ia,ib-2) )
                     End If
                     If (ia.ge.2.and.ib.ge.2) Then
                        rV4Int(iZeta,iCar,ia,ib) =
     &                     rV4Int(iZeta,iCar,ia,ib)
     &                    +Dble(ia*(ia-1)*ib*(ib-1)) *
     &                                       Sxyz(iZeta,iCar,ia-2,ib-2)
                     End If
*
                  End Do
               End Do
         End Do
      End Do
*
      If (iPrint.ge.99) Then
         Do ib = 0, nb
            Do ia = 0, na
               Write (Label,'(A,I2,A,I2,A)')
     &            'In MVe: rV2Int(',ia,',',ib,',1)'
               Call RecPrt(Label,' ',rV2Int(1,1,ia,ib,1),nZeta,3)
               Write (Label,'(A,I2,A,I2,A)')
     &            'In MVe: rV2Int(',ia,',',ib,',2)'
               Call RecPrt(Label,' ',rV2Int(1,1,ia,ib,2),nZeta,3)
               Write (Label,'(A,I2,A,I2,A)')
     &            'In MVe: rV4Int(',ia,',',ib,')'
               Call RecPrt(Label,' ',rV4Int(1,1,ia,ib),nZeta,3)
            End Do
         End Do
      End If
*
*     Call GetMem(' Exit MVe   ','CHECK','REAL',iDum,iDum)
      Call qExit('MVe ')
      Return
      End
