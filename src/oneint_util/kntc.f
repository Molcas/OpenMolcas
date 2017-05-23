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
      SubRoutine Kntc(Txyz,Sxyz,na,nb,Alpha,Beta,nZeta)
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
      Real*8 Txyz(nZeta,3,0:na,0:nb), Sxyz(nZeta,3,0:na+1,0:nb+1),
     &       Alpha(nZeta), Beta(nZeta)
      Character*80 Label
*
      iRout = 115
      iPrint = nPrint(iRout)
*     Call qEnter('Kntc')
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In Kntc: Alpha',' ',Alpha,nZeta,1)
         Call RecPrt(' In Kntc: Beta ',' ',Beta ,nZeta,1)
         Do 1 ia = 0, na+1
            Do 2 ib = 0, nb+1
               Write (Label,'(A,I2,A,I2,A)')
     &               ' In Kntc: Sxyz(',ia,',',ib,')'
               Call RecPrt(Label,' ',Sxyz(1,1,ia,ib),nZeta,3)
 2          Continue
 1       Continue
      End If
      Do 10 ia = 0, na
         Do 20 ib = 0, nb
            If (ia.eq.0 .and. ib.eq.0) Then
               Do 31 iCar = 1, 3
                  Do 41 iZeta = 1, nZeta
                     Txyz(iZeta,iCar,ia,ib) =
     &                   Two * Alpha(iZeta) * Beta(iZeta) *
     &                                       Sxyz(iZeta,iCar,ia+1,ib+1)
 41               Continue
 31            Continue
            Else If (ia.eq.0) Then
               Do 32 iCar = 1, 3
                  Do 42 iZeta = 1, nZeta
                     Txyz(iZeta,iCar,ia,ib) =
     &                   Two * Alpha(iZeta) * Beta(iZeta) *
     &                                       Sxyz(iZeta,iCar,ia+1,ib+1)
     &           - Alpha(iZeta) * DBLE(ib) * Sxyz(iZeta,iCar,ia+1,ib-1)
 42               Continue
 32            Continue
            Else If (ib.eq.0) Then
               Do 33 iCar = 1, 3
                  Do 43 iZeta = 1, nZeta
                     Txyz(iZeta,iCar,ia,ib) =
     &                   Two * Alpha(iZeta) * Beta(iZeta) *
     &                                       Sxyz(iZeta,iCar,ia+1,ib+1)
     &           - Beta(iZeta)  * DBLE(ia) * Sxyz(iZeta,iCar,ia-1,ib+1)
 43               Continue
 33            Continue
            Else
               Do 30 iCar = 1, 3
                  Do 40 iZeta = 1, nZeta
                     Txyz(iZeta,iCar,ia,ib) =
     &                Half * DBLE(ia * ib) * Sxyz(iZeta,iCar,ia-1,ib-1)
     &           - Beta(iZeta)  * DBLE(ia) * Sxyz(iZeta,iCar,ia-1,ib+1)
     &           - Alpha(iZeta) * DBLE(ib) * Sxyz(iZeta,iCar,ia+1,ib-1)
     &                 + Two * Alpha(iZeta) * Beta(iZeta) *
     &                                       Sxyz(iZeta,iCar,ia+1,ib+1)
 40               Continue
 30            Continue
            End If
            If (iPrint.ge.99) Then
               Write (Label,'(A,I2,A,I2,A)') ' In Kntc: Txyz(',ia,',',
     &                ib,')'
               Call RecPrt(Label,' ',Txyz(1,1,ia,ib),nZeta,3)
            End If
 20      Continue
 10   Continue
*
*     Call GetMem(' Exit Kntc  ','CHECK','REAL',iDum,iDum)
*     Call qExit('Kntc')
      Return
      End
