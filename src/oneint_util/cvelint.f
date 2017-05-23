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
      SubRoutine CVelInt(Vxyz,Sxyz,na,nb,Alpha,Beta,nZeta)
************************************************************************
*                                                                      *
* Object: to assemble the cartesian components of the velocity inte-   *
*         grals from the cartesian components of the overlap integals. *
*                                                                      *
* Called from: PrpInt                                                  *
*                                                                      *
* Calling    : CRecPrt                                                 *
*               RecPrt                                                 *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Complex*16 Vxyz(nZeta,3,0:na,0:nb,2), Sxyz(nZeta,3,0:na+1,0:nb+1)
      Real*8 Alpha(nZeta), Beta(nZeta)
      Character*80 Label
      iRout=160
      iPrint = nPrint(iRout)

      If (iPrint.ge.99) Then
         Call RecPrt(' In CVelInt: Beta ',' ',Beta ,nZeta,1)
      End If
      Do 10 ia = 0, na
         Do 20 ib = 0, nb
            If (ia.ne.0 .and. ib.ne.0) Then
               Do iCar = 1, 3
                  Do iZeta = 1, nZeta
                     Vxyz(iZeta,iCar,ia,ib,1) =
     &                    Dble(ia) * Sxyz(iZeta,iCar,ia-1,ib)
     &                 - Alpha(iZeta) * Two * Sxyz(iZeta,iCar,ia+1,ib)
                     Vxyz(iZeta,iCar,ia,ib,2) =
     &                    Dble(ib) * Sxyz(iZeta,iCar,ia,ib-1)
     &                 - Beta(iZeta) * Two * Sxyz(iZeta,iCar,ia,ib+1)
                  End Do
               End Do
            Else If (ia.eq.0 .and. ib.ne.0) Then
               Do iCar = 1, 3
                  Do iZeta = 1, nZeta
                     Vxyz(iZeta,iCar,ia,ib,1) =
     &                 - Alpha(iZeta) * Two * Sxyz(iZeta,iCar,ia+1,ib)
                     Vxyz(iZeta,iCar,ia,ib,2) =
     &                    Dble(ib) * Sxyz(iZeta,iCar,ia,ib-1)
     &                 - Beta(iZeta) * Two * Sxyz(iZeta,iCar,ia,ib+1)
                  End Do
               End Do
            Else If (ia.ne.0 .and. ib.eq.0) Then
               Do iCar = 1, 3
                  Do iZeta = 1, nZeta
                     Vxyz(iZeta,iCar,ia,ib,1) =
     &                    Dble(ia) * Sxyz(iZeta,iCar,ia-1,ib)
     &                 - Alpha(iZeta) * Two * Sxyz(iZeta,iCar,ia+1,ib)
                     Vxyz(iZeta,iCar,ia,ib,2) =
     &                 - Beta(iZeta) * Two * Sxyz(iZeta,iCar,ia,ib+1)
                  End Do
               End Do
            Else
               Do iCar = 1, 3
                  Do iZeta = 1, nZeta
                     Vxyz(iZeta,iCar,ia,ib,1) =
     &                 - Alpha(iZeta) * Two * Sxyz(iZeta,iCar,ia+1,ib)
                     Vxyz(iZeta,iCar,ia,ib,2) =
     &                 - Beta(iZeta) * Two * Sxyz(iZeta,iCar,ia,ib+1)
                  End Do
               End Do
            End If
*
            If (iPrint.ge.99) Then
               Write (Label,'(A,I2,A,I2,A)') ' In CVelInt: Vxyz(',
     &                ia,',',ib,',1)'
               Call CRecPrt(Label,' ',Vxyz(1,1,ia,ib,1),nZeta,3,'R')
               Call CRecPrt(Label,' ',Vxyz(1,1,ia,ib,1),nZeta,3,'I')
               Write (Label,'(A,I2,A,I2,A)') ' In CVelInt: Vxyz(',
     &                ia,',',ib,',2)'
               Call CRecPrt(Label,' ',Vxyz(1,1,ia,ib,2),nZeta,3,'R')
               Call CRecPrt(Label,' ',Vxyz(1,1,ia,ib,2),nZeta,3,'I')
            End If
 20      Continue
 10   Continue

      Return
      End
