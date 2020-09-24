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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine Assemble_dTdmu(nZeta,Final,la,lb,Elalbp,Elalbm,Beta)
************************************************************************
*                                                                      *
* Object: to assemble the diamagnetic shielding integrals from         *
*         electric field integrals.                                    *
*                                                                      *
* Called from: DMSInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             February '91                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8  Final (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3),
     &        Elalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,3),
     &        Elalbm(nZeta,(la+1)*(la+2)/2,(lb  )*(lb+1)/2,3),
     &        Beta(nZeta)
      Character*80 Label
*
*     Statement function for cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
      nElem(ix) = (ix+1)*(ix+2)/2
*
      iRout = 231
      iPrint = nPrint(iRout)
*
*     Fact = -1.D6 * One2C2
      If (iPrint.ge.99) Then
          Write (6,*) ' In Assemble_dTdmu la,lb=',la,lb
          Do ia = 1, nElem(la)
             Do ib = 1, nElem(lb+1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalbp(',ia,',',ib,',x)'
                Call RecPrt(Label,' ',Elalbp(1,ia,ib,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalbp(',ia,',',ib,',y)'
                Call RecPrt(Label,' ',Elalbp(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalbp(',ia,',',ib,',z)'
                Call RecPrt(Label,' ',Elalbp(1,ia,ib,3),nZeta,1)
             End Do
          End Do
          Do ia = 1, nElem(la)
             ib_max=nElem(lb-1)
             If (lb.eq.0) ib_max=0
             Do ib = 1, ib_max
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalbm(',ia,',',ib,',x)'
                Call RecPrt(Label,' ',Elalbm(1,ia,ib,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalbm(',ia,',',ib,',y)'
                Call RecPrt(Label,' ',Elalbm(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Elalbm(',ia,',',ib,',z)'
                Call RecPrt(Label,' ',Elalbm(1,ia,ib,3),nZeta,1)
             End Do
          End Do
      End If
*
      Do ixa = la, 0, -1
         Do iya = la-ixa, 0, -1
            iza = la-ixa-iya
            ipa = Ind(la,ixa,iza)
*
      Do ixb = lb, 0, -1
         Do iyb = lb-ixb, 0, -1
            izb = lb-ixb-iyb
            ipb = Ind(lb,ixb,izb)
*
            Do iZeta = 1, nZeta
               xyTmp=-Two*Beta(nzeta)*Elalbp(iZeta,ipa,
     &                Ind(lb+1,ixb  ,izb  ),1)
               yxTmp=-Two*Beta(nzeta)*Elalbp(iZeta,ipa,
     &                Ind(lb+1,ixb+1,izb  ),2)
               yzTmp=-Two*Beta(nzeta)*Elalbp(iZeta,ipa,
     &                Ind(lb+1,ixb  ,izb+1),2)
               zyTmp=-Two*Beta(nzeta)*Elalbp(iZeta,ipa,
     &                Ind(lb+1,ixb  ,izb  ),3)
               zxTmp=-Two*Beta(nzeta)*Elalbp(iZeta,ipa,
     &                Ind(lb+1,ixb+1,izb  ),3)
               xzTmp=-Two*Beta(nzeta)*Elalbp(iZeta,ipa,
     &                Ind(lb+1,ixb  ,izb+1),1)
               If (ixb.ge.1) Then
                  yxTmp = yxTmp + Dble(ixb)*Elalbm(iZeta,ipa,
     &                                      Ind(lb-1,ixb-1,izb  ),2)
                  zxTmp = zxTmp + Dble(ixb)*Elalbm(iZeta,ipa,
     &                                      Ind(lb-1,ixb-1,izb  ),3)
               End If
               If (iyb.ge.1) Then
                  xyTmp = xyTmp + Dble(iyb)*Elalbm(iZeta,ipa,
     &                                      Ind(lb-1,ixb  ,izb  ),1)
                  zyTmp = xyTmp + Dble(iyb)*Elalbm(iZeta,ipa,
     &                                      Ind(lb-1,ixb  ,izb  ),3)
               End If
               If (izb.ge.1) Then
                  xzTmp = xzTmp + Dble(izb)*Elalbm(iZeta,ipa,
     &                                      Ind(lb-1,ixb  ,izb-1),1)
                  yzTmp = yzTmp + Dble(izb)*Elalbm(iZeta,ipa,
     &                                      Ind(lb-1,ixb  ,izb-1),2)
               End If
               Final(iZeta,ipa,ipb,1) = -(xyTmp - yxTmp)
               Final(iZeta,ipa,ipb,2) = -(yzTmp - zyTmp)
               Final(iZeta,ipa,ipb,3) = -(zxTmp - xzTmp)
            End Do
*
         End Do
      End Do
*
         End Do
      End Do
*
      If (iPrint.ge.49) Then
          Do iComp = 1, 3
             Write (Label,'(A,I2,A)') ' Final (',iComp,') '
             Call RecPrt(Label,' ',Final(1,1,1,iComp),nZeta,
     &                   nElem(la)*nELem(lb))
          End Do
      End If
*
      Return
      End
