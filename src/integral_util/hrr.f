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
      SubRoutine HRR(la,lb,A,B,Target,nPrim,nTrgt,ipIn)
************************************************************************
*                                                                      *
* Object: to generate the contracted integrals (one or two-electron)   *
*         from (a+b,0| or (0,a+b| type of integrals. The approach      *
*         implemented here is the same as that one described by Rys,   *
*         Dupuis and King. In this they implemented the method for the *
*         2D-integrals. The method was later implemented to be applied *
*         on the contracted integrals( Head-Gordon & Pople).           *
*                                                                      *
* Called from: OneEl                                                   *
*              TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              OCHRR                                                   *
*              HRR1                                                    *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             February '90                                             *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden.                                         *
*             Modified to not use pointers June '91                    *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Target(nPrim,nTrgt), A(3), B(3), AB(3)
*
*     Statement function for canonical indices
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      Ind1(ixyz,ix,iz) = ixyz*(ixyz+1)*(ixyz+2)/6
     &                 + (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
      iRout =25
      iPrint = nPrint(iRout)
*     Fast exit if HRR will not be applied.
*     Call GetMem('Enter_HRR','Check','Real',iDum,iDum)
      If (la.eq.0 .or. lb.eq.0) Then
         ipIn = 1
*        Call GetMem('Fxit_HRR','Check','Real',iDum,iDum)
         Return
      End If
*
      AB(1) = A(1) - B(1)
      AB(2) = A(2) - B(2)
      AB(3) = A(3) - B(3)
      If (lb.gt.la) Then
         AB(1) = -AB(1)
         AB(2) = -AB(2)
         AB(3) = -AB(3)
      End If
      ABSqrt=Sqrt(AB(1)**2+AB(2)**2+AB(3)**2)
      If (ABSqrt.eq.Zero) Then
         Call OCHRR(Target,nPrim,nTrgt,la,lb,ipIn)
      Else
*
         ib1Max=Min(la,lb)
         iaMin =Max(la,lb)
         ipRslt = 0
         Do 10 ib1 = 1, ib1Max
            ib=ib1-1
            iaMax = la+lb-ib1
            Do 15 ia = iaMax, iaMin, -1
               ia1=ia+1
               If (Mod(ib1,2).eq.0) Then
*                 Low loading of target integrals.
                  iab1 = nElem(ib1)*(Ind1(ia,ia,0)-Ind1(iaMin,iaMin,0))
*                 High access of source integrals.
                  iab  =  nTrgt - nElem(ib)*(Ind1(iaMax+1,0,iaMax+1) -
     &                    Ind1(ia,ia,0)+1)
                  ia1b =  nTrgt - nElem(ib)*(Ind1(iaMax+1,0,iaMax+1) -
     &                    Ind1(ia1,ia1,0)+1)
               Else
*                 High loading of target integrals.
                  iab1 = nTrgt - nElem(ib1)*(Ind1(iaMax,0,iaMax) -
     &                   Ind1(ia,ia,0)+1)
*                 Low access of source integrals.
                  iab  = nElem(ib)*(Ind1(ia,ia,0)-Ind1(iaMin,iaMin,0))
                  ia1b = nElem(ib)*(Ind1(ia1,ia1,0)-Ind1(iaMin,iaMin,0))
               End If
               ipRslt = iab1
*
*              Generate this block of integrals with the HRR
*
               Call HRR1(Target(1,iab1+1),nElem(ia)*nElem(ib1),
     &                   Target(1,ia1b+1),nElem(ia1)*nElem(ib),AB,
     &                   Target(1,iab+1),nElem(ia)*nElem(ib),
     &                   ia,ib,ia1,ib1,nPrim,la,lb)
*
 15         Continue
 10      Continue
         ipIn = nPrim*ipRslt + 1
      End If
*
      Return
      End
