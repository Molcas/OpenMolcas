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
      SubRoutine Ass_pX(Alpha,nZeta,Final,la,lb,Slaplb,Slamlb,nComp)
************************************************************************
*                                                                      *
* Object: to assemble the pV integrals from                            *
*         derivative integrals of the electric potential.              *
*                                                                      *
* Called from: pvInt                                                   *
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
      Real*8  Final (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3,nComp),
     *        Slaplb(nZeta,(la+2)*(la+3)/2,(lb+1)*(lb+2)/2,nComp),
     *        Slamlb(nZeta, la   *(la+1)/2,(lb+1)*(lb+2)/2,nComp),
     *        Alpha(nZeta)
      Character*80 Label
*
*     Statement function for cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
      nElem(ix) = (ix+1)*(ix+2)/2
*
      iRout = 203
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
          Write (6,*)
          Write (6,*) ' In Ass_pX la,lb,nComp=',la,lb,nComp
          Write (6,*)
          Call RecPrt('Alpha','(10G15.8)',Alpha,nZeta,1)
          Do iComp = 1, nComp
             Write (6,*)
             Write (6,*) 'iComp=',iComp
             Write (6,*)
             Write (Label,'(A,I2,A)')
     &             'Ass_pX:  Slaplb(iComp=',iComp,')'
             Call RecPrt(Label,'(10f15.8)',Slaplb(1,1,1,iComp),
     &                   nZeta,nElem(la+1)*nElem(lb))
             If (la.gt.0) Then
                Write (Label,'(A,I2,A)')
     &                'Ass_pX: Slamlb(iComp=,',iComp,')'
                Call RecPrt(Label,'(10G15.8)',Slamlb(1,1,1,iComp),
     &                      nZeta,nElem(la-1)*nElem(lb))
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
*                                                                      *
************************************************************************
*                                                                      *
            If (ixa.EQ.0) Then
               ixp=Ind(la+1,ixa+1,iza)
               Do iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,1,iComp) =
     &               Two*Alpha(iZeta)*Slaplb(iZeta,ixp,ipb,iComp)
               End Do
            Else
               ixp=Ind(la+1,ixa+1,iza)
               ixm=Ind(la-1,ixa-1,iza)
               Do iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,1,iComp) =
     &               Two*Alpha(iZeta)*Slaplb(iZeta,ixp,ipb,iComp)
     &              -Dble(ixa)*       Slamlb(iZeta,ixm,ipb,iComp)
               End Do
            End If
*
            If (iya.EQ.0) Then
               iyp=Ind(la+1,ixa,iza)
               Do iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,2,iComp) =
     &               Two*Alpha(iZeta)*Slaplb(iZeta,iyp,ipb,iComp)
               End Do
            Else
               iyp=Ind(la+1,ixa,iza)
               iym=Ind(la-1,ixa,iza)
               Do iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,2,iComp) =
     &               Two*Alpha(iZeta)*Slaplb(iZeta,iyp,ipb,iComp)
     &              -Dble(iya)*       Slamlb(iZeta,iym,ipb,iComp)
               End Do
            End If
*
            If (iza.EQ.0) Then
               izp=Ind(la+1,ixa,iza+1)
               Do iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,3,iComp) =
     &               Two*Alpha(iZeta)*Slaplb(iZeta,izp,ipb,iComp)
               End Do
            Else
               izp=Ind(la+1,ixa,iza+1)
               izm=Ind(la-1,ixa,iza-1)
               Do iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,3,iComp) =
     &               Two*Alpha(iZeta)*Slaplb(iZeta,izp,ipb,iComp)
     &              -Dble(iza)*       Slamlb(iZeta,izm,ipb,iComp)
               End Do
            End If
*                                                                      *
************************************************************************
*                                                                      *
21       Continue
20    Continue
*
11       Continue
10    Continue
*
      End Do
*
      If (iPrint.ge.49) Then
          Write (6,*) ' In Ass_pX la,lb,nComp=',la,lb,nComp
          Do iComp = 1, nComp
             Write (6,*)
             Write (6,*) 'iComp=',iComp
             Write (6,*)
*
             Write (Label,'(A,I2,A)')
     &          ' Ass_pX: pX( 1,iComp=',iComp,')'
             Call RecPrt(Label,'(10G15.8)',
     &                   Final(1,1,1,1,iComp),nZeta,
     &                   nElem(la)*nElem(lb))
*
             Write (Label,'(A,I2,A)')
     &          ' Ass_pX: pX( 2,iComp=',iComp,')'
             Call RecPrt(Label,'(10G15.8)',
     &                   Final(1,1,1,2,iComp),nZeta,
     &                   nElem(la)*nElem(lb))
*
             Write (Label,'(A,I2,A)')
     &          ' Ass_pX: pX( 3,iComp=',iComp,')'
             Call RecPrt(Label,'(10G15.8)',
     &                   Final(1,1,1,3,iComp),nZeta,
     &                   nElem(la)*nElem(lb))
         End Do
      End If
*
      Return
      End
