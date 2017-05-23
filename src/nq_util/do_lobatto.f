************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Do_Lobatto(L_Eff,nPoints,ipR)
************************************************************************
*                                                                      *
*     Computes datas useful for the angular quadrature.                *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "nq_info.fh"
#include "real.fh"
#include "WrkSpc.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Statement Function
*
      Pax(i) = Work(ip_O-1+i)
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate angular grid a la Lobatto
*
      nPoints=0
      nTheta = (L_Eff+3)/2
      Do iTheta = 1, nTheta
         nPhi = L_Eff
         If(iTheta.eq.1.or.iTheta.eq.nTheta) nPhi = 1
         If(iTheta.eq.nTheta/2+1.and.nTheta/2*2-nTheta.eq.-1.and.
     &                         nTheta.gt.3)  nPhi = L_Eff+4
*
         nPoints = nPoints + nPhi
*
      End Do
*
      Call GetMem('AngRW','Allo','Real',ipR,4*nPoints)
*
      nTheta = (L_Eff+3)/2
      nLabatto=3*(nTheta+2)*(nTheta+3)/2
      Call GetMem('Labatto','Allo','Real',ipLabatto,nLabatto)
      Call Lobatto(ntheta,Work(ipLabatto))
*
      mTheta=nTheta-1
      iOffT=ipLabatto + 3*mTheta*(mTheta+1)/2
      iOff   = ipR
      Do iTheta = 1, nTheta
*
         Cos_Theta=Work(iOffT)
         Sin_Theta=Sqrt(One-Cos_Theta**2)
         w_Theta  =Work(iOffT+1)
         iOffT = iOffT + 3
*
         nPhi = L_Eff
         if(iTheta.eq.1.or.iTheta.eq.nTheta) nPhi = 1
         if(iTheta.eq.nTheta/2+1.and.nTheta/2*2-nTheta.eq.-1.and.
     &                         nTheta.gt.3)  nPhi = L_Eff+4
*
         Do iPhi = 1, nPhi
            Call Phi_point(iPhi,nPhi,Cos_Phi,Sin_Phi,w_Phi)
*
            x = Sin_Theta*Cos_Phi
            y = Sin_Theta*Sin_Phi
            z = Cos_Theta
            Work(iOff  )=Pax(1)*x+Pax(4)*y+Pax(7)*z
            Work(iOff+1)=Pax(2)*x+Pax(5)*y+Pax(8)*z
            Work(iOff+2)=Pax(3)*x+Pax(6)*y+Pax(9)*z
            Work(iOff+3)=w_Theta*w_Phi
            iOff = iOff + 4
*
         End Do  ! iPhi
*
      End Do     ! iTheta
      Call GetMem('Labatto','Free','Real',ipLabatto,nLabatto)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
