!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      Function EstI(Zeta,rKapAB,nAlpha,nBeta,Coeff1,niBas,Coeff2,njBas,
     &              xab,nab,Scrt,nScrt,IndZ)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************
      use Constants, only: Zero
      Implicit None
      Real*8 EstI
      Integer nAlpha,nBeta,niBas,njBas,nab,nScrt
      Real*8 Zeta(nAlpha*nBeta), rKapAB(nAlpha,nBeta),
     &       Coeff1(nAlpha,niBas), Coeff2(nBeta,njBas),
     &       xab(nAlpha*nBeta), Scrt(nScrt)
      Integer IndZ(nAlpha*nBeta+1)

      Integer mZeta, iZeta, iAlpha, iBeta, iEta, iDelta, iGamma, iBas,
     &        jBas, ijBas, iHigh
      Integer, External:: iDAMax_
      Real*8 rab, rcd, AInt
!
#ifdef _DEBUGPRINT_
      Write (6,*) 'Esti:mZeta=',IndZ(nAlpha*nBeta)
      Call RecPrt('Esti:xab',' ',xab,1,nAlpha*nBeta)
      Call RecPrt('Esti:Coeff1',' ',Coeff1,nAlpha,niBas)
      Call RecPrt('Esti:Coeff2',' ',Coeff2,nBeta ,njBas)
#endif
!
      mZeta=IndZ(nAlpha*nBeta+1)
      call dcopy_(niBas*njBas,[Zero],0,Scrt,1)
      Do iZeta = 1, mZeta
         iBeta = (IndZ(iZeta)-1)/nAlpha + 1
         iAlpha = IndZ(iZeta) - (iBeta-1)*nAlpha
         rab = xab(iZeta)
         Do iEta = 1, mZeta
            iDelta = (IndZ(iEta)-1)/nAlpha + 1
            iGamma = IndZ(iEta) - (iDelta-1)*nAlpha
            rcd = xab(iEta)
            AInt=  rab*rcd
            Do iBas = 1, niBas
               Do jBas = 1, njBas
                  ijBas = (jBas-1)*niBas + iBas
                  Scrt(ijBas) = Scrt(ijBas) +
     &             Abs(Coeff1(iAlpha,iBas) * Coeff2(iBeta, jBas))*
     &             Abs(Coeff1(iGamma,iBas) * Coeff2(iDelta,jBas))
     &             * AInt
               End Do
            End Do
         End Do
      End Do
      iHigh = iDAMax_(niBas*njBas,Scrt,1)
      EstI = Sqrt(Scrt(iHigh))
!
      Return
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Zeta)
         Call Unused_real_array(rKapAB)
         Call Unused_integer(nab)
      End If
      End Function EstI
