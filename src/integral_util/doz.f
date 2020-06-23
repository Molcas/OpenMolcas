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
      SubRoutine DoZ   (Alpha,nAlpha,Beta,nBeta,A,B,P,Zeta,ZInv,
     &                  rKappa,IndZt,IncZet,SkipZt,Data,IndZ,
     &                  iphX,iphY,iphZ)
************************************************************************
* Object : to compute zeta, P and kappa.                               *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*             May '90 (modified for integral cutoff)                   *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden. Modified to pick up precomputed entities*
*             rather than to compute them. July '91.                   *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
      Real*8 Alpha(nAlpha), Beta(nBeta), Zeta(nAlpha*nBeta),
     &       A(3),          B(3),        P(nAlpha*nBeta,3),
     &       Data(nAlpha*nBeta*(nDarray-1)), rKappa(nAlpha*nBeta),
     &       ZInv(nAlpha*nBeta)
      Integer IndZt(nAlpha*nBeta), IndZ(nAlpha*nBeta)
      Logical SkipZt
*
      iRout = 68
      iPrint = nPrint(iRout)
*     Call qEnter('DoZ')
*
*     Transfer precomputer data
*
      nZeta=nAlpha*nBeta
      mZeta = IndZ(nAlpha*nBeta)
      SkipZt = mZeta.eq.0
      Call ICopy(mZeta,IndZ,1,IndZt,1)
      IndZt(nAlpha*nBeta)=mZeta
      call dcopy_(mZeta,Data(ip_Z    (1,nZeta)),1,Zeta,1)
      call dcopy_(mZeta,Data(ip_Kappa(1,nZeta)),1,rKappa,1)
      call dcopy_(mZeta,Data(ip_Pcoor(1,nZeta)),1,P(1,1),1)
      iOff=1+nZeta
      call dcopy_(mZeta,Data(ip_Pcoor(iOff,nZeta)),1,P(1,2),1)
      iOff=1+2*nZeta
      call dcopy_(mZeta,Data(ip_Pcoor(iOff,nZeta)),1,P(1,3),1)
      call dcopy_(mZeta,Data(ip_ZInv(1,nZeta)),1,ZInv,1)
*     Fix the phase due to operations on the whole charge density.
      If (iphY.ne.1) Call DScal_(mZeta,One*DBLE(iphY),P(1,2),1)
      If (iphX.ne.1) Call DScal_(mZeta,One*DBLE(iphX),P(1,1),1)
      If (iphZ.ne.1) Call DScal_(mZeta,One*DBLE(iphZ),P(1,3),1)
*
      If (iPrint.ge.99) Then
         Write (6,*) ' In DoZ'
         Call RecPrt(' Zeta',' ',Zeta,mZeta,1)
         Call RecPrt(' ZInv',' ',ZInv,mZeta,1)
         Call RecPrt(' Kappa',' ',rKappa,mZeta,1)
         Call RecPrt(' Px',' ',P(1,1),mZeta,1)
         Call RecPrt(' Py',' ',P(1,2),mZeta,1)
         Call RecPrt(' Pz',' ',P(1,3),mZeta,1)
         Write (6,*) ' phase factors=',iphX, iphY, iphZ
         Write (6,*) ' IndZt=',IndZt
         Call RecPrt(' Data',' ',Data,nAlpha*nBeta,nDArray)
         Write (6,*) ' Exit DoZ'
      End If
*
*     Call qExit('DoZ')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(A)
         Call Unused_real_array(B)
         Call Unused_integer(IncZet)
      End If
      End
