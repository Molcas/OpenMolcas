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
* Copyright (C) 1994, Bernd Artur Hess                                 *
************************************************************************
      SubRoutine Ass_pXp(Beta,nZeta,Final,la,lb,Slalbp,Slalbm,nComp)
************************************************************************
*                                                                      *
* Object: to assemble the pVp integrals
*                                                                      *
* Called from: PVPInt                                                  *
*                                                                      *
* Calling    : qEnter                                                  *
*              qExit                                                   *
*                                                                      *
*     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
*             Chemie, University of Bonn, Germany, August 1994         *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8  Final (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     *        Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,3,nComp),
     *        Slalbm(nZeta,(la+1)*(la+2)/2, lb   *(lb+1)/2,3,nComp),
     *        Beta(nZeta)
      Character*80 Label
*
*     Statement function for cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
      nElem(ix) = (ix+1)*(ix+2)/2
*
      iRout = 211
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Write (6,*)
         Write (6,*) ' In Ass_pXp la,lb,nComp,=',la,lb,nComp
         Write (6,*)
         Call RecPrt('Beta','(10G15.8)',Beta,nZeta,1)
         Do iComp = 1, nComp
            Write (6,*) 'iComp=',iComp
            Write (Label,'(A,I2,A)')
     &            ' Ass_pXp: Slalbp(1,iComp=',iComp,')'
            Call RecPrt(Label,'(10G15.8)',Slalbp(1,1,1,1,iComp),
     &                   nZeta,nElem(la)*nElem(lb+1))
            Write (Label,'(A,I2,A)')
     &            ' Ass_pXp: Slalbp(2,iComp=',iComp,')'
            Call RecPrt(Label,'(10G15.8)',Slalbp(1,1,1,2,iComp),
     &                   nZeta,nElem(la)*nElem(lb+1))
            Write (Label,'(A,I2,A)')
     &            ' Ass_pXp: Slalbp(3,iComp=',iComp,')'
            Call RecPrt(Label,'(10G15.8)',Slalbp(1,1,1,3,iComp),
     &                   nZeta,nElem(la)*nElem(lb+1))
         If (lb.gt.0) Then
            Write (Label,'(A,I2,A)')
     &             'Ass_pXp: Slalbm(1,iComp=',iComp,')'
            Call RecPrt(Label,'(10G15.8)',Slalbm(1,1,1,1,iComp),
     &                  nZeta,nElem(la)*nElem(lb-1))
            Write (Label,'(A,I2,A)')
     &             'Ass_pXp: Slalbm(2,iComp=',iComp,')'
            Call RecPrt(Label,'(10G15.8)',Slalbm(1,1,1,2,iComp),
     &                  nZeta,nElem(la)*nElem(lb-1))
            Write (Label,'(A,I2,A)')
     &             'Ass_pXp: Slalbm(3,iComp=',iComp,')'
            Call RecPrt(Label,'(10G15.8)',Slalbm(1,1,1,3,iComp),
     &                  nZeta,nElem(la)*nElem(lb-1))
         End If
         End Do
      End If
*
      Do iComp = 1, nComp
*
      Do 10 ixa = la, 0, -1
         Do 11 iya = la-ixa, 0, -1
            iza = la-ixa-iya
            ipa = Ind(la,ixa,iza)
*
      Do 20 ixb = lb, 0, -1
         Do 21 iyb = lb-ixb, 0, -1
            izb = lb-ixb-iyb
            ipb = Ind(lb,ixb,izb)
*
            Do 30 iZeta = 1, nZeta
               Final(iZeta,ipa,ipb,iComp) =
     &                  Two*Beta(iZeta) *
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb  ),1,iComp)
     &                 +Two*Beta(iZeta) *
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb  ,izb  ),2,iComp)
     &                 +Two*Beta(iZeta) *
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb  ,izb+1),3,iComp)
30          Continue
*
            If (ixb.gt.0) Then
               Do 31 iZeta = 1, nZeta
                 Final(iZeta,ipa,ipb,iComp)=Final(iZeta,ipa,ipb,iComp)
     &          -Dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),1,iComp)
31             Continue
            End If
*
            If (iyb.gt.0) Then
               Do 32 iZeta = 1, nZeta
                 Final(iZeta,ipa,ipb,iComp)=Final(iZeta,ipa,ipb,iComp)
     &          -Dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),2,iComp)
32             Continue
            End If
*
            If (izb.gt.0) Then
               Do 33 iZeta = 1, nZeta
                 Final(iZeta,ipa,ipb,iComp)=Final(iZeta,ipa,ipb,iComp)
     &          -Dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),3,iComp)
33             Continue
            End If
*
21       Continue
20    Continue
*
11       Continue
10    Continue
*
      End Do ! iComp
*
      If (iPrint.ge.49) Then
         Do iComp = 1, nComp
            Write (Label,'(A,I2,A,I2,A,I2,A)')
     &            ' Ass_pXp: pXp(iComp=',iComp,')'
            Call RecPrt(Label,'(10G15.8)',Final(1,1,1,iComp),
     &                              nZeta,nElem(la)*nElem(lb))
         End Do
      End If
*
      Return
      End
