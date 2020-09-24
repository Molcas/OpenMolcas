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
      SubRoutine TERI(Zeta,Eta,P,Q,rKapab,rKapcd,T,Fact,ZEInv,nT,IsChi,
     &                ChiI2)
************************************************************************
*                                                                      *
* Object: to entities for the two-electron integrals which are used in *
*         in the Rys quadrature to evaluate these integrals.           *
*                                                                      *
* Called from: Rys                                                     *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Intrinsic Sqrt
#include "print.fh"
#include "real.fh"
      Real*8 Zeta(nT), Eta(nT), P(nT,3), Q(nT,3),
     &       rKapab(nT), rKapcd(nT),
     &       T(nT), Fact(nT), ZEInv(nT)
*
      iRout = 56
*
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt(' Zeta in TERI',' ',Zeta,1,nT)
      Call RecPrt(' Eta in TERI',' ',Eta,1,nT)
      Call RecPrt(' P in TERI',' ',P,nT,3)
      Call RecPrt(' Q in TERI',' ',Q,nT,3)
      Call RecPrt(' Kab in TERI',' ',rKapab,1,nT)
      Call RecPrt(' Kcd in TERI',' ',rKapcd,1,nT)
#endif
*
      Do iT = 1, nT
         tmp = 1.0D0/(Zeta(iT)+Eta(iT)
     >       +(Eta(iT)*Zeta(iT)*ChiI2)*Dble(IsChi))
         ZEInv(iT) = tmp
         Rho = Zeta(iT)*Eta(iT)*tmp
         PQ2 = (P(iT,1)-Q(iT,1))**2
     &       + (P(iT,2)-Q(iT,2))**2
     &       + (P(iT,3)-Q(iT,3))**2
         T(iT) = Rho*PQ2
         Fact(iT) =  rKapab(iT) * rKapcd(iT) * Sqrt(tmp)
       End Do
#ifdef _DEBUGPRINT_
       Call RecPrt('Tvalue',' ',T,1,nT)
       Call RecPrt('Fact  ',' ',Fact,1,nT)
#endif
      Return
      End
