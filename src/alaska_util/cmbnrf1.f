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
* Copyright (C) 1991,1992,1995, Roland Lindh                           *
************************************************************************
      SubRoutine CmbnRF1(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,Final,nComp,
     &                   Fact,Temp,Alpha,Beta,Grad,nGrad,DAO,IfGrad,
     &                   IndGrd,iStab,jStab,kOp,EF)
************************************************************************
*                                                                      *
* Object: to compute gradient integrals for SC Reaction Fields         *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             Modified for reaction field calculations July '92        *
*             Modified for gradient calculations May '95               *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChBas
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Integer IndGrd(3,2), kOp(2)
      Logical IfGrad(3,2)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp,6),
     &       Zeta(nZeta), rKappa(nZeta), Fact(nZeta), Temp(nZeta),
     &       Rnxyz(nZeta,3,0:la+1,0:lb+1,0:lr), Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2),
     &       Alpha(nZeta), Beta(nZeta), EF(nComp)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
      iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
      nElem(i) = (i+1)*(i+2)/2
*
      iRout = 134
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' In CmbnRF1: EF',' ',EF,nComp,1)
      End If
*
      tTwo=Two
*
      Do 130 iZeta = 1, nZeta
         Fact(iZeta) = rKappa(iZeta) * Zeta(iZeta)**(-Three/Two)
130   Continue
*
*---- Loop over angular components of the basis set
*
      Do 10 ixa = 0, la
         iyaMax=la-ixa
      Do 11 ixb = 0, lb
         iybMax=lb-ixb
         Do 20 iya = 0, iyaMax
            iza = la-ixa-iya
            ipa= Ind(la,ixa,iza)
         Do 21 iyb = 0, iybMax
            izb = lb-ixb-iyb
            ipb= Ind(lb,ixb,izb)
*
*           Combine multipole moment integrals
*
            If (IfGrad(1,1)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     If (ixa.gt.0) Then
                        xa = -ixa
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   (tTwo*Alpha(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa+1,ixb,ix) +
     &                                xa*Rnxyz(iZeta,1,ixa-1,ixb,ix)) *
     &                                   Rnxyz(iZeta,2,iya,iyb,iy)
                        End Do
                     Else
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   tTwo*Alpha(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa+1,ixb,ix) *
     &                                   Rnxyz(iZeta,2,iya,iyb,iy)
                        End Do
                     End If
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)
                        Do iZeta = 1, nZeta
                           Final(iZeta,ipa,ipb,iComp,1) = Temp(iZeta)*
     &                             Rnxyz(iZeta,3,iza,izb,iz)
                        End Do
                     End Do
                  End Do
               End Do
            End If
            If (IfGrad(1,2)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     If (ixb.gt.0) Then
                        xb = -ixb
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   (tTwo*Beta(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb+1,ix) +
     &                                xb*Rnxyz(iZeta,1,ixa,ixb-1,ix)) *
     &                                   Rnxyz(iZeta,2,iya,iyb,iy)
                        End Do
                     Else
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   tTwo*Beta(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb+1,ix) *
     &                                   Rnxyz(iZeta,2,iya,iyb,iy)
                        End Do
                     End If
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)
                        Do iZeta = 1, nZeta
                           Final(iZeta,ipa,ipb,iComp,4) = Temp(iZeta)*
     &                             Rnxyz(iZeta,3,iza,izb,iz)
                        End Do
                     End Do
                  End Do
               End Do
            End If
            If (IfGrad(2,1)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     If (iya.gt.0) Then
                        ya = -iya
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb,ix) *
     &                                   (tTwo*Alpha(iZeta) *
     &                                   Rnxyz(iZeta,2,iya+1,iyb,iy) +
     &                                ya*Rnxyz(iZeta,2,iya-1,iyb,iy))
                        End Do
                     Else
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb,ix) *
     &                                   tTwo*Alpha(iZeta) *
     &                                   Rnxyz(iZeta,2,iya+1,iyb,iy)
                        End Do
                     End If
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)
                        Do iZeta = 1, nZeta
                           Final(iZeta,ipa,ipb,iComp,2) = Temp(iZeta)*
     &                             Rnxyz(iZeta,3,iza,izb,iz)
                        End Do
                     End Do
                  End Do
               End Do
            End If
            If (IfGrad(2,2)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     If (iyb.gt.0) Then
                        yb = -iyb
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb,ix) *
     &                                   (tTwo*Beta(iZeta) *
     &                                   Rnxyz(iZeta,2,iya,iyb+1,iy) +
     &                                yb*Rnxyz(iZeta,2,iya,iyb-1,iy))
                        End Do
                     Else
                        Do iZeta = 1, nZeta
                           Temp(iZeta) = Fact(iZeta) *
     &                                   Rnxyz(iZeta,1,ixa,ixb,ix) *
     &                                   tTwo*Beta(iZeta) *
     &                                   Rnxyz(iZeta,2,iya,iyb+1,iy)
                        End Do
                     End If
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)
                        Do iZeta = 1, nZeta
                           Final(iZeta,ipa,ipb,iComp,5) = Temp(iZeta)*
     &                             Rnxyz(iZeta,3,iza,izb,iz)
                        End Do
                     End Do
                  End Do
               End Do
            End If
            If (IfGrad(3,1)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     Do iZeta = 1, nZeta
                        Temp(iZeta) = Fact(iZeta) *
     &                                Rnxyz(iZeta,1,ixa,ixb,ix)*
     &                                Rnxyz(iZeta,2,iya,iyb,iy)
                     End Do
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)
                        If (iza.gt.0) Then
                           za = -iza
                           Do iZeta = 1, nZeta
                              Final(iZeta,ipa,ipb,iComp,3)= Temp(iZeta)*
     &                                (tTwo*Alpha(iZeta) *
     &                                Rnxyz(iZeta,3,iza+1,izb,iz) +
     &                             za*Rnxyz(iZeta,3,iza-1,izb,iz))
                           End Do
                        Else
                           Do iZeta = 1, nZeta
                              Final(iZeta,ipa,ipb,iComp,3)= Temp(iZeta)*
     &                                tTwo*Alpha(iZeta) *
     &                                Rnxyz(iZeta,3,iza+1,izb,iz)
                           End Do
                        End If
                     End Do
                  End Do
               End Do
            End If
            If (IfGrad(3,2)) Then
               Do ix = 0, lr
                  Do iy = 0, lr-ix
                     Do iZeta = 1, nZeta
                        Temp(iZeta) = Fact(iZeta) *
     &                                Rnxyz(iZeta,1,ixa,ixb,ix)*
     &                                Rnxyz(iZeta,2,iya,iyb,iy)
                     End Do
*
                     Do ir = ix+iy, lr
                        iz = ir-ix-iy
                        iComp=Ind(ir,ix,iz)+iOff(ir)
                        If (izb.gt.0) Then
                           zb = -izb
                           Do iZeta = 1, nZeta
                              Final(iZeta,ipa,ipb,iComp,6)= Temp(iZeta)*
     &                                (tTwo* Beta(iZeta) *
     &                                Rnxyz(iZeta,3,iza,izb+1,iz) +
     &                             zb*Rnxyz(iZeta,3,iza,izb-1,iz))
                           End Do
                        Else
                           Do iZeta = 1, nZeta
                              Final(iZeta,ipa,ipb,iComp,6)= Temp(iZeta)*
     &                                tTwo* Beta(iZeta) *
     &                                Rnxyz(iZeta,3,iza,izb+1,iz)
                           End Do
                        End If
                     End Do
                  End Do
               End Do
            End If
*
 21      Continue
 20      Continue
 11   Continue
 10   Continue
      If (iPrint.ge.99) Then
         Call RecPrt('In CmbnRF1: DAO',' ',DAO,
     &               nZeta,nElem(la)*nElem(lb))
         Call RecPrt('In CmbnRF1: Final',' ',Final,
     &               nZeta*nElem(la)*nElem(lb)*nComp,6)
      End If
*
*-----Trace the gradient integrals
*
      nDAO = nZeta * (la+1)*(la+2)/2 * (lb+1)*(lb+2)/2
      Do iEF = 1, nComp
         Do iCn = 1, 2
            Do iCar = 1, 3
               If (IndGrd(iCar,iCn).ne.0) Then
*                 Accumulate contibution to the gradient
                  iGrad = Abs(IndGrd(iCar,iCn))
                  If (iCn.eq.1) Then
                     i1 = iCar
                     i2 = iCar + 3
                     ps = DBLE( iPrmt( kOp(1), iChBas(1+iCar) ) )
                     Fct = DBLE(iStab)/DBLE(nIrrep)
                  Else
                     i1 = iCar + 3
                     i2 = iCar
                     ps = DBLE( iPrmt( kOp(2), iChBas(1+iCar) ) )
                     Fct = ps * DBLE(jStab)/DBLE(nIrrep)
                  End If
                  If (IndGrd(iCar,iCn).lt.0) Then
*                    Gradient via the translational invariance.
                     Grad(iGrad) = Grad(iGrad) - Fct * EF(iEF) *
     &                   DDot_(nDAO,DAO,1,Final(1,1,1,iEF,i2),1)
                  Else
                     Grad(iGrad) = Grad(iGrad) + Fct * EF(iEF) *
     &                   DDot_(nDAO,DAO,1,Final(1,1,1,iEF,i1),1)
                  End If
               End If
            End Do  ! End loop over cartesian components, iCar
         End Do     ! End loop over centers, iCn
      End Do        ! End loop over EF components, iEF
*
      Return
      End
