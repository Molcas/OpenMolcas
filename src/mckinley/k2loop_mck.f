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
* Copyright (C) 1990-1992, Roland Lindh                                *
*               1990, IBM                                              *
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine k2Loop_mck(Coor,
     &                      iAnga,iCmpa,
     &                      iDCRR,nDCRR,Data, ijCmp,
     &                      Alpha,nAlpha,Beta,nBeta,
     &                      Coeff1,iBasn,Coeff2,jBasn,
     &                      nMemab,Con,
     &                      Wk002,m002,Wk003,m003,Wk004,m004,
     &                      iStb,jStb,
     &                      ipTmp1,ipTmp2,ipTmp3,
     &                      ipKnew,ipLnew,ipPnew,ipQnew)
************************************************************************
*                                                                      *
* Object: to compute zeta, kappa, and P.                               *
*         This is done for all unique pairs of centers                 *
*         generated from the symmetry unique centers A and B.          *
*                                                                      *
* Called from: Drvk2                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy   (ESSL)                                          *
*              DoZeta                                                  *
*              SchInt                                                  *
*              PckInt                                                  *
*              GetMem                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             June '91, modified to compute zeta, P, kappa and inte-   *
*             grals for Schwartz inequality in a k2 loop.              *
*             January '92 modified to gradient calculations.           *
*             April '92, modified to use the Cauchy-Schwarz inequality *
*              to estimate the integral derivatives.                   *
*                                                                      *
*             May   '95  modified (simplified) for hessian calculation *
*             By Anders Bernhardsson                                   *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "disp2.fh"
      Real*8 Coor(3,2), CoorM(3,4), Alpha(nAlpha), Beta(nBeta),
     &       Data(nAlpha*nBeta*nDArray+nDScalar,nDCRR),
     &       Coeff1(nAlpha,iBasn), Coeff2(nBeta,jBasn),
     &       Wk002(m002),Wk003(m003),Wk004(m004), Con(nAlpha*nBeta)
      Integer iDCRR(0:7), iAnga(4), iCmpa(4), mStb(2)
*
      Call k2Loop_mck_internal(Data)
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Con)
         Call Unused_real_array(Wk004)
         Call Unused_integer(ipTmp1)
         Call Unused_integer(ipTmp2)
         Call Unused_integer(ipTmp3)
         Call Unused_integer(ipKnew)
         Call Unused_integer(ipLnew)
         Call Unused_integer(ipPnew)
         Call Unused_integer(ipQnew)
      End If
*
*     This is to allow type punning without an explicit interface
      Contains
      SubRoutine k2Loop_mck_internal(Data)
      Use Iso_C_Binding
      Real*8, Target :: Data(nAlpha*nBeta*nDArray+nDScalar,nDCRR)
      Integer, Pointer :: iData(:)
      nZeta = nAlpha*nBeta
      mStb(1) = iStb
      mStb(2) = jStb
*
      CoorM(1,1) = Coor(1,1)
      CoorM(2,1) = Coor(2,1)
      CoorM(3,1) = Coor(3,1)
*
      Do 100 lDCRR = 0, nDCRR-1
*
         Call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
         call dcopy_(6,CoorM(1,1),1,CoorM(1,3),1)
*
*--------Compute Zeta, P and kappa.
*
         Call C_F_Pointer(C_Loc(Data(ip_IndZ(1,nZeta),lDCRR+1)),
     &                    iData,[nAlpha*nBeta+1])
         Call DoZeta(Alpha,nAlpha,Beta,nBeta,
     &               CoorM(1,1),CoorM(1,2),
     &               Data(ip_PCoor(1,nZeta),lDCRR+1),
     &               Data(ip_Z    (1,nZeta),lDCRR+1),
     &               Data(ip_Kappa(1,nZeta),lDCRR+1),
     &               Data(ip_ZInv (1,nZeta),lDCRR+1),
     &               Data(ip_Alpha(1,nZeta,1),lDCRR+1),
     &               Data(ip_Beta (1,nZeta,2),lDCRR+1),
     &               iData)
         Nullify(iData)
*
         Call SchInt_mck(CoorM,iAnga,iCmpa,nAlpha,nBeta,nMemab,
     &                   Data(ip_Z(1,nZeta),lDCRR+1),
     &                   Data(ip_ZInv(1,nZeta),lDCRR+1),
     &                   Data(ip_Kappa(1,nZeta),lDCRR+1),
     &                   Data(ip_PCoor(1,nZeta),lDCRR+1),
     &                   nZeta,Wk002,m002,Wk003,m003)
*
         Call PckInt_mck(Wk002,nZeta,ijCmp,
     &                   Data(ip_ab(1,nZeta),lDCRR+1),
     &                   Data(ip_Z(1,nZeta),lDCRR+1))
*                                                                      *
************************************************************************
*                                                                      *
*        Estimate the largest contracted integral.
*
         Call C_F_Pointer(C_Loc(Data(ip_IndZ(1,nZeta),lDCRR+1)),
     &                    iData,[nAlpha*nBeta+1])
         Data(ip_EstI(nZeta),lDCRR+1) =
     &                      EstI(Data(ip_Z(1,nZeta),lDCRR+1),
     &                           Data(ip_Kappa(1,nZeta),lDCRR+1),
     &                           nAlpha,nBeta,
     &                           Coeff1,iBasn,Coeff2,jBasn,
     &                           Data(ip_ab   (1,nZeta),lDCRR+1),
     &                           iCmpa(1)*iCmpa(2),
     &                           Wk002,m002,
     &                           iData)
*                                                                      *
************************************************************************
*                                                                      *
*------- Find the largest integral estimate (AO Basis).
*
         Tst  = -One
         Do  iZeta = 0, nZeta-1
             Tst=Max(Data(ip_Z(iZeta+1,nZeta),lDCRR+1),Tst)
         End Do
         Data(ip_ZetaM(nZeta),lDCRR+1) = tst
*
         Tst  = -One
         ZtMax=Zero
         abMax=Zero
         Do  iZeta = 1, nZeta
            tmp =          Data(ip_ab(iZeta,nZeta),lDCRR+1)
            If (Tst.lt.tmp) Then
               Tst = tmp
               ZtMax = Data(ip_Z (iZeta,nZeta),lDCRR+1)
               abMax = Data(ip_ab(iZeta,nZeta),lDCRR+1)
            End If
         End Do
         Data(ip_ZtMax(nZeta),lDCRR+1) = ZtMax
         Data(ip_abMax(nZeta),lDCRR+1) = abMax
 100  Continue
*
      Return
      End SubRoutine k2Loop_mck_internal
*
      End
