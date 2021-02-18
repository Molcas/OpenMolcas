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
      SubRoutine TNAI(Zeta,Eta,P,Q,rKapab,rKapcd,T,Fact,ZEInv,nT,IsChi,
     &                ChiI2)
************************************************************************
*                                                                      *
* Object: to entities for the nucelar attraction integrals which are   *
*         used in the Rys quadrature to evaluate these integrals.      *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Zeta(nT), Eta(nT), P(nT,3), Q(nT,3),
     &       rKapab(nT), rKapcd(nT), ZEInv(nT), T(nT), Fact(nT)
*
#ifdef _DEBUGPRINT_
      iRout = 57
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' Zeta in TNAI',' ',Zeta,nT,1)
         Call RecPrt(' Eta in TNAI',' ',Eta,nT,1)
         Call RecPrt(' P in TNAI',' ',P,nT,3)
         Call RecPrt(' Q in TNAI',' ',Q,nT,3)
         Call RecPrt(' Kab in TNAI',' ',rKapab,nT,1)
         Call RecPrt(' Kcd in TNAI',' ',rKapcd,nT,1)
         Write (6,*) ' In TNAI: ABeqCD=',ABeqCD
      End If
#endif
      Do iT = 1, nT
         PQ2 = (P(iT,1)-Q(iT,1))**2
     &       + (P(iT,2)-Q(iT,2))**2
     &       + (P(iT,3)-Q(iT,3))**2
         T(iT) = Zeta(iT)*PQ2
         ZEInv(iT) = 1.0D0/Zeta(iT)
         Fact(iT) =  2.0D0*rKapab(iT)*Pi/Zeta(iT)
      End Do
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Call RecPrt('Tvalue',' ',T,nT,1)
         Call RecPrt('Fact  ',' ',Fact,nT,1)
      End If
#endif
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Eta)
         Call Unused_real_array(rKapcd)
         Call Unused_integer(IsChi)
         Call Unused_real(ChiI2)
      End If

      Return
      End
