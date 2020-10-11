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
* Copyright (C) 1990,1992, Roland Lindh                                *
************************************************************************
      SubRoutine TERISq(Zeta,Eta,P,Q,rKapab,rKapcd,T,Fact,ZEInv,nT,
     &                  IsChi,ChiI2)
************************************************************************
*                                                                      *
* Object: to entities for the two-electron integrals which are used in *
*         in the Rys quadrature to evaluate these integrals.           *
*                                                                      *
*         OBSERVE that the prefactor is only partial!!!                *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             May '92 Modified to fit 2nd order differential scheme for*
*             the gradient estimates.                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Zeta(nT), Eta(nT), P(nT,3), Q(nT,3),
     &       rKapab(nT), rKapcd(nT),
     &       T(nT), Fact(nT), ZEInv(nT)
*
      iRout = 56
      iPrint = nPrint(iRout)
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Call RecPrt(' Zeta in TERISq',' ',Zeta,nT,1)
         Call RecPrt(' P in TERISq',' ',P,nT,3)
         Call RecPrt(' Q in TERISq',' ',Q,nT,3)
         Call RecPrt(' Kab in TERISq',' ',rKapab,nT,1)
         Call RecPrt(' Kcd in TERISq',' ',rKapcd,nT,1)
      End If
#endif
*
      Do iT = 1, nT
         tmp = 1.0D0/(Zeta(iT)+Zeta(iT)
     >            +(Zeta(iT)*Zeta(iT)*ChiI2)*Dble(IsChi))
         ZEInv(iT) = tmp
         Rho = Zeta(iT)*Zeta(iT)*tmp
         PQ2 = (P(iT,1)-Q(iT,1))**2
     &       + (P(iT,2)-Q(iT,2))**2
     &       + (P(iT,3)-Q(iT,3))**2
         T(iT) = Rho*PQ2
         Fact(iT) =  rKapab(iT) * rKapcd(iT)
      End Do
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Call RecPrt('Tvalue',' ',T,nT,1)
         Call RecPrt('Fact  ',' ',Fact,nT,1)
      End If
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Eta)
      End
