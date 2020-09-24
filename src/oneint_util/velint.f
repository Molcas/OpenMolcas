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
      SubRoutine VelInt(Vxyz,Sxyz,na,nb,Beta,nZeta)
************************************************************************
*                                                                      *
* Object: to assemble the cartesian components of the velocity inte-   *
*         grals from the cartesian components of the overlap integals. *
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
#include "real.fh"
#include "print.fh"
      Real*8 Vxyz(nZeta,3,0:na,0:nb), Sxyz(nZeta,3,0:na,0:nb+1),
     &       Beta(nZeta)
      Character*80 Label
      iRout=160
      iPrint = nPrint(iRout)

      If (iPrint.ge.99) Then
         Call RecPrt(' In VelInt: Beta ',' ',Beta ,nZeta,1)
      End If
      Do 10 ia = 0, na
         Do 20 ib = 0, nb
            If (ib.eq.0) Then
               Do 33 iCar = 1, 3
                  Do 43 iZeta = 1, nZeta
                     Vxyz(iZeta,iCar,ia,ib) =
     &                 - Beta(iZeta) * Two * Sxyz(iZeta,iCar,ia,ib+1)
 43               Continue
 33            Continue
            Else
               Do 30 iCar = 1, 3
                  Do 40 iZeta = 1, nZeta
                     Vxyz(iZeta,iCar,ia,ib) =
     &                    Dble(ib) * Sxyz(iZeta,iCar,ia,ib-1)
     &                 - Beta(iZeta) * Two * Sxyz(iZeta,iCar,ia,ib+1)
 40               Continue
 30            Continue
            End If
*
            If (iPrint.ge.99) Then
               Write (Label,'(A,I2,A,I2,A)') ' In VelInt: Vxyz(',ia,',',
     &                ib,')'
               Call RecPrt(Label,' ',Vxyz(1,1,ia,ib),nZeta,3)
            End If
 20      Continue
 10   Continue

*     Call GetMem(' Exit VelInt  ','CHECK','REAL',iDum,iDum)
      Return
      End
