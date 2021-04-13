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
* Copyright (C) 2010, Jonas Bostrom                                    *
************************************************************************
      Subroutine Compute_B_1(irc,jDen,jSym,kSym,lSym,jVec1,nVec,
     &                       nBasFnc)
******************************************************************
*                                                                *
*     Author Jonas Bostrom, June 2010                            *
*                                                                *
*     Purpose: To pickup C_lk^J for a specific J from disk and   *
*              Rectangularize.                                   *
*                                                                *
******************************************************************
      use ExTerm, only: CijK
      Implicit Real*8 (a-h,o-z)
#include "exterm.fh"
#include "chomp2g_alaska.fh"
#include "WrkSpc.fh"

      If(imp2prpt .eq. 2) Then
         lBVec = nBasFnc*nBasFnc*nVec
         Do i = 1, 2
            iAdr = 1 + nBasFnc*nBasFnc*(jVec1-1)
            Call dDaFile(LuBVector(i),2,Work(ip_B_mp2(i)),lBVec,iAdr)
         End Do
      End If

      lCVec = nIJR(kSym,lSym,1)*nVec
      iAdr = nIJR(kSym,lSym,1)*(jVec1-1) + iAdrCVec(jSym,kSym,jDen)
      Call dDaFile(LuCVector(jSym,jDen),2,CijK,lCVec,iAdr)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(irc)
      End
*
******************************************************************
******************************************************************
******************************************************************
*
      Subroutine Compute_B_2(irc,lSOl,jVec1,kSym,lSym,nVec)
******************************************************************
*                                                                *
*     Author Jonas Bostrom, June 2010                            *
*                                                                *
*     Purpose: To do the first halftransformation                *
*              Ckl^J -> Cml^J .                                  *
*                                                                *
******************************************************************
      use ExTerm, only: CijK, VJ
      Implicit Real*8 (a-h,o-z)
#include "exterm.fh"
#include "WrkSpc.fh"
      Real*8, Pointer :: CijK_Rec(:,:)=>Null()
      Real*8, Pointer :: CmjK(:)=>Null()
*
      ijList(iSym,jSym,i,j) = ipijList + iChOrbR(iSym,jSym,1)
     &                 + i + (j-1)*(nChOrb(iSym-1,1)+1)
*
      jDen = 1
      ipMO = ip_CMOi(jDen) + iOff_CMOi(lSym,jDen) +
     &       (lSOl-1)*nChOrb(lSym-1,jDen)

      iE = jVec1*nChOrb(lSym-1,jDen)*nChOrb(kSym-1,jDen)
      n = nChOrb(lSym-1,1)
      m = nChOrb(kSym-1,1)

      iS = iE + 1
      iE = iE + n*m
      CijK_Rec(1:n,1:m) => CijK(iS:iE)

      CmjK(1:m) => VJ(1:m)

      Call FZero(CmjK,nChOrb(kSym-1,jDen))

      Do iI = 1, nChOrb(kSym-1,1)
         Do iJ = 1, iWork(ijList(lSym,kSym,0,iI))
            iJ1 = iWork(ijList(lSym,kSym,iJ,iI))
            CmjK(iI) = CmjK(iI) +
     &           CijK_Rec(iJ1,iI) * Work(ipMO + iJ1-1)
         End Do
      End Do

      CijK_Rec=>Null()
      CmjK=>Null()
      irc=0

c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nVec)
      End If
      End
*
******************************************************************
******************************************************************
******************************************************************
*
      Real*8 Function Compute_B_3(irc,kSym,lSym,kSOk)
******************************************************************
*                                                                *
*     Author Jonas Bostrom, June 2010                            *
*                                                                *
*     Purpose: To do the second halftransformation               *
*              Cmj^K -> Bmd^K .                                  *
*                                                                *
******************************************************************
      use ExTerm, only: VJ
      Implicit Real*8 (a-h,o-z)
#include "exterm.fh"
#include "WrkSpc.fh"
      Real*8, Pointer :: CmjK(:)=>Null()

      jDen = 1
*
      n = nChOrb(kSym-1,jDen)

      ipMO = ip_CMOi(jDen) + iOff_CMOi(kSym,jDen) + (kSOk-1)*n

      CmjK(1:n) => VJ(1:n)

      Compute_B_3=ddot_(n,Work(ipMO),1,CmjK,1)

      CmjK=>Null()

      irc=0
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(lSym)
      End If
      End
*
******************************************************************
******************************************************************
*
      Real*8 Function Compute_B_4(irc,kSOk,lSOl,jAOj,nBasFnc,iOpt)
******************************************************************
*                                                                *
*     Author Jonas Bostrom, June 2010                            *
*                                                                *
*     Purpose: To do part of MP2 gradient.                       *
*                                                                *
******************************************************************
      Implicit Real*8 (a-h,o-z)
#include "exterm.fh"
#include "chomp2g_alaska.fh"
#include "WrkSpc.fh"


      B_mp2 = 0.0d0
      iOff1 = (jAOj)*nBasFnc*nBasFnc + (kSOk-1)*nBasFnc + lSOl-1
      iOff2 = (jAOj)*nBasFnc*nBasFnc + (lSOl-1)*nBasFnc + kSOk-1
      B_mp2 = B_mp2 + (Work(ip_B_mp2(iOpt)+iOff1)+
     &                 Work(ip_B_mp2(iOpt)+iOff2))/2.0d0
      Compute_B_4 = B_mp2
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(irc)
      End
