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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine TERIS(Zeta,Eta,P,Q,rKapab,rKapcd,T,Fact,ZEInv,nT,IsChi,
     &                 ChiI2)
************************************************************************
*                                                                      *
* Object: compute the arguments for the reduced list of integrals which*
*         are used in prescreening.                                    *
*                                                                      *
* Called from: Rys                                                     *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             June '91, modified for k2 loop.                          *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Zeta(nT), P(nT,3), rKapab(nT), T(nT), Fact(nT),
     &       ZEInv(nT)
*
      iRout = 244
      iPrint = nPrint(iRout)
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Call RecPrt(' Zeta in TERIS',' ',Zeta,nT,1)
         Call RecPrt(' P in TERIS',' ',P,nT,3)
         Call RecPrt(' Kab in TERIS',' ',rKapab,nT,1)
      End If
#endif
*
      Do iT = 1, nT
         T(iT) = 0.0D0
         tmp = 1.0D0/(Zeta(iT)+Zeta(iT)
     >         +(Zeta(iT)*Zeta(iT)*ChiI2)*Dble(IsChi))
         ZEInv(iT) = tmp
         Fact(iT) =  rKapab(iT) **2 * Sqrt(tmp)
      End Do
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Call RecPrt('In TERIS: Tvalue',' ',T,nT,1)
         Call RecPrt('In TERIS: Fact  ',' ',Fact,nT,1)
      End If
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(P)
#endif
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(Eta)
         Call Unused_real(Q)
         Call Unused_real(rKapcd)
      End If
      End
