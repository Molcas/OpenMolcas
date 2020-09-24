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
      SubRoutine Rysg1(iAnga,nRys,nT,
     &                 Alpha,Beta,Gamma,Delta,
     &                 Zeta,ZInv,nZeta,Eta,EInv,nEta,
     &                 P,lP,Q,lQ,Coori,Coora,CoorAC,
     &                 Array,nArray,
     &                 Tvalue,ModU2,Cff2D,
     &                 PAO,nPAO,Grad,nGrad,IfGrad,IndGrd,kOp,iuvwx)
************************************************************************
*                                                                      *
* Object: to compute the gradient of the two-electron integrals.       *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              Tvalue                                                  *
*              RtsWgh                                                  *
*              vRysRW                                                  *
*              ModU2                                                   *
*              Cff2D                                                   *
*              Rys2Dm                                                  *
*              HrrCtl                                                  *
*              Rys2Dg                                                  *
*              Assg1                                                   *
*              Distg1                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Modified to 1st order derivatives October '91            *
************************************************************************
      use vRys_RW
      use Symmetry_Info, only: iOper
      use Real_Info, only: ChiI2
      use Temporary_Parameters, only: IsChi
      Implicit Real*8 (A-H,O-Z)
      External Tvalue, ModU2, Cff2D
      External Exp_1, Exp_2
#include "real.fh"
#include "notab.fh"
#include "print.fh"
      Real*8 Zeta(nZeta), ZInv(nZeta), P(lP,3),
     &       Eta(nEta),   EInv(nEta),  Q(lQ,3),
     &       Alpha(nZeta), Beta(nZeta), Gamma(nEta), Delta(nEta),
     &       CoorAC(3,2), Coora(3,4), Coori(3,4), Array(nArray),
     &       PAO(nT,nPAO), Grad(nGrad), Temp(9)
      Integer iAnga(4), IndGrd(3,4), Index(3,4),
     &          kOp(4), iuvwx(4),   JndGrd(3,4), lOp(4)
      Logical AeqB, CeqD, AeqC, EQ, IfGrad(3,4),
     &          JfGrad(3,4)
*
      lOp(1) = iOper(kOp(1))
      lOp(2) = iOper(kOp(2))
      lOp(3) = iOper(kOp(3))
      lOp(4) = iOper(kOp(4))
#ifdef _DEBUGPRINT_
      call dcopy_(lP-nZeta,[Zero],0,P(nZeta+1,1),1)
      call dcopy_(lP-nZeta,[Zero],0,P(nZeta+1,2),1)
      call dcopy_(lP-nZeta,[Zero],0,P(nZeta+1,3),1)
      Call RecPrt(' In Rysg1:P',' ',P,lP,3)
      call dcopy_(lQ-nEta,[Zero],0,Q(nEta+1,1),1)
      call dcopy_(lQ-nEta,[Zero],0,Q(nEta+1,2),1)
      call dcopy_(lQ-nEta,[Zero],0,Q(nEta+1,3),1)
      Call RecPrt(' In Rysg1:Q',' ',Q,lQ,3)
      Call RecPrt(' In Rysg1:Zeta',' ',Zeta,nZeta,1)
      Call RecPrt(' In Rysg1:ZInv',' ',ZInv,nZeta,1)
      Call RecPrt(' In Rysg1:Eta',' ',Eta,nEta,1)
      Call RecPrt(' In Rysg1:EInv',' ',EInv,nEta,1)
      Call RecPrt(' In Rysg1:Alpha',' ',Alpha,nZeta,1)
      Call RecPrt(' In Rysg1:Beta ',' ',Beta ,nZeta,1)
      Call RecPrt(' In Rysg1:Gamma',' ',Gamma,nEta,1)
      Call RecPrt(' In Rysg1:Delta',' ',Delta,nEta,1)
      Call RecPrt(' In Rysg1:Coora',' ',Coora,3,4)
      Call RecPrt(' In Rysg1:Coori',' ',Coori,3,4)
      Call RecPrt(' In Rysg1:CoorAC',' ',CoorAC,3,2)
      Write (6,*) ' In Rysg1: iAnga=',iAnga
#endif
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      AeqB = EQ(Coora(1,1),Coora(1,2))
      CeqD = EQ(Coora(1,3),Coora(1,4))
      AeqC = EQ(Coora(1,1),Coora(1,3))
      lla = 0
      llb = 0
      llc = 0
      lld = 0
      Do 10 i = 1, 3
         If (IfGrad(i,1)) Then
            lla = 1
         End If
         If (IfGrad(i,2)) Then
            llb = 1
         End If
         If (IfGrad(i,3)) Then
            llc = 1
         End If
         If (IfGrad(i,4)) Then
            lld = 1
         End If
 10   Continue
      lab = Max(lla,llb)
      lcd = Max(llc,lld)
      nabMax = la + lb + lab
      ncdMax = lc + ld + lcd
*
*     Allocate memory for the integral gradients.
*
      ip = 1
*     ipAC = ip
*     ip = ip + nT*nPAO * 9
*
*     Allocate memory for the 2D-integrals.
*
      ip2D0 = ip
      n2D0 = Max( (nabMax+1)*(ncdMax+1),
     &            (la+2)*(lb+2)*(ncdMax+1),
     &            (la+2)*(lb+2)*(lc+2)*(ld+2))
      ip = ip + n2D0*3*nT*nRys
*
*     Allocate memory for the 1st order derivatives of the 2D-integrals.
*
      ip2D1 = ip
      n2D1 = Max( (nabMax+1)*(ncdMax+1),
     &            (la+2)*(lb+2)*(ncdMax+1),
     &            (la+1)*(lb+1)*(lc+1)*(ld+1)*3)
      ip = ip + n2D1*3*nT*nRys
*
*     Allocate memory for the coefficients in the recurrence relation
*     of the 2D-integrals.
*
      nTR=nT*nRys
      ipPAQP = ip
      ip = ip + nTR*3
      ipQCPQ = ip
      ip = ip + nTR*3
      ipB10 = ip
      lB10=Max(Min(nabMax-1,1),0)
      ip = ip + nTR*3*lB10
      labMax = Min(nabMax,ncdMax)
      ipB00 = ip
      lB00=Max(Min(labMax,1),0)
      ip = ip + nTR*3*lB00
      ipB01 = ip
      lB01=Max(Min(ncdMax-1,1),0)
      ip = ip + nTR*3*lB01
*     Allocate memory for the roots.
      ipU2 = ip
      ip = ip + nT*nRys
*     Allocate memory for Zeta, ZInv, Eta, EInv
      ipZeta = ip
      ip = ip + nT
      ipEta  = ip
      ip = ip + nT
      ipZInv = ip
      ip = ip + nT
      ipEInv = ip
      ip = ip + nT
*     Allocate memory for P and Q
      ipP = ip
      ip = ip + 3*nT
      ipQ = ip
      ip = ip + 3*nT
*     Allocate memory for the inverse.
      ipDiv = ip
      ip = ip + nT
*     Allocate memory for the arguments.
      ipTv = ip
      ip = ip + nT
*define _CHECK_
#ifdef _CHECK_
      If (ip-1.gt.nArray) Then
         Call WarningMessage(2,'Rysg1: ip-1 =/= nArray (pos.1)')
         Write (6,*) ' nArray=',nArray
         Write (6,*) ' ip-1  =',ip-1
         Write (6,*) ' nRys  =',nRys
         Write (6,*) ' nZeta =',nZeta
         Write (6,*) ' nEta  =',nEta
         Call Abend()
      End If
#endif
*
*     Expand Zeta, ZInv, Eta ,EInv, P, and Q
*
      Do iEta = 1, nEta
         iOff = (iEta-1)*nZeta
         call dcopy_(nZeta,Zeta,  1,Array(iOff+ipZeta  ),1)
         call dcopy_(nZeta,ZInv,  1,Array(iOff+ipZInv  ),1)
         call dcopy_(nZeta,P(1,1),1,Array(iOff+ipP     ),1)
         iOff = iOff + nZeta*nEta
         call dcopy_(nZeta,P(1,2),1,Array(iOff+ipP     ),1)
         iOff = iOff + nZeta*nEta
         call dcopy_(nZeta,P(1,3),1,Array(iOff+ipP     ),1)
      End Do
      Do iZeta = 1, nZeta
         iOff = iZeta-1
         call dcopy_(nEta, Eta,  1,Array(iOff+ipEta) ,nZeta)
         call dcopy_(nEta,EInv,  1,Array(iOff+ipEInv),nZeta)
         call dcopy_(nEta,Q(1,1),1,Array(iOff+ipQ     ),nZeta)
         iOff = iOff + nZeta*nEta
         call dcopy_(nEta,Q(1,2),1,Array(iOff+ipQ     ),nZeta)
         iOff = iOff + nZeta*nEta
         call dcopy_(nEta,Q(1,3),1,Array(iOff+ipQ     ),nZeta)
      End Do
*
*     Compute tha arguments for which we will compute the roots and
*     the weights.
*
      Call Tvalue(Array(ipZeta),Array(ipEta),Array(ipP),Array(ipQ),nT,
     &            Array(ipTv),Array(ipDiv),IsChi,ChiI2)
*
*     Compute roots and weights. Make sure that the weights ends up in
*     the array where the z component of the 2D integrals will be.
*     Call vRysRW if roots and weights are tabulated in various Taylor
*     expansions. If not tabulated call RtsWgh.
*
*     Pointer to z-component of 2D-integrals where the weights will be
*     put directly. This corresponds to xyz2D(1,1,3,0,0).
      ipWgh = ip2D0 + 2*nT*nRys
      If (nRys.gt.nMxRys .or. NoTab) Then
#ifdef _CHECK_
         If (ip-1.gt.nArray) Then
            Call WarningMessage(2,'Rysg1: ip-1 =/= nArray (pos.2)')
            Write (6,*) ' nArray=',nArray
            Write (6,*) ' ip-1  =',ip-1
            Call Abend()
         End If
#endif
         Call RtsWgh(Array(ipTv),nT,Array(ipU2),Array(ipWgh),nRys)
      Else
#ifdef _CHECK_
         If (ip-1.gt.nArray) Then
            Call WarningMessage(2,'Rysg1: ip-1 =/= nArray (pos.3)')
            Write (6,*) ' nArray=',nArray
            Write (6,*) ' ip-1  =',ip-1
            Call Abend()
         End If
#endif
         Call vRysRW(la+1,lb,lc,ld,Array(ipTv),Array(ipU2),Array(ipWgh),
     &               nT,nRys)
      End If
*-----Drop ipTv
      ip = ip - nT
*
*     Modify the roots.
*
      Call ModU2(Array(ipU2),nT,nRys,Array(ipDiv))
*-----Drop ipDiv
      ip = ip - nT
*
*     Compute coefficients for the recurrence relations of the
*     2D-integrals
*
      Call Cff2D(Max(nabMax-1,0),Max(ncdMax-1,0),nRys,
     &           Array(ipZeta),Array(ipZInv),Array(ipEta),Array(ipEInv),
     &           nT,Coori,CoorAC,Array(ipP),Array(ipQ),
     &           la+lab,lb,lc+lcd,ld,
     &           Array(ipU2),Array(ipPAQP),Array(ipQCPQ),
     &           Array(ipB10),Array(ipB00),labMax,Array(ipB01))
*-----Drop ipU2
      ip = ip - nT*nRys
*     Let go of Zeta, ZInv, Eta, and EInv
      ip = ip - nT*4
*     Let go of P and Q
      ip = ip - 6*nT
*
*     Compute the intermediate 2D-integrals from the roots and weights.
*
      Call vRys2Dm(Array(ip2D0),nT,nRys,nabMax,
     &            ncdMax,Array(ipPAQP),Array(ipQCPQ),
     &            Array(ipB10),Max(nabMax-1,0),
     &            Array(ipB00),labMax,
     &            Array(ipB01),Max(ncdMax-1,0),
     &            la,lb,lc,ld,IfGrad)
*-----Drop ipB01
      ip = ip - nTR*3*lB01
*-----Drop ipB00
      ip = ip - nTR*3*lB00
*-----Drop ipB10
      ip = ip - nTR*3*lB10
*-----Drop ipQCPQ
      ip = ip - nTR*3
*-----Drop ipPAQP
      ip = ip - nTR*3
*
*-----Apply the transfer equation to the intermediate 2D-integrals.
*
      Call HrrCtl(Array(ip2D0),n2D0,Array(ip2D1),n2D1,
     &            la,lb,lc,ld,nabmax,ncdmax,
     &            nTR,Coora(1,1),Coora(1,2),Coora(1,3),Coora(1,4),
     &            IfGrad)
*
*     Compute the gradients of the 2D-integrals. Copy some information
*     which will be modified. This has to be done in order to facilitate
*     partitioning.
*
      ipScr = ip
      ip = ip + nT*nRys
      ipTmp = ip
      ip = ip + nT
      Call ICopy(12,IndGrd,1,JndGrd,1)
      Do 8877 i = 1, 3
         Do 7788 j = 1, 4
            JfGrad(i,j) = IfGrad(i,j)
 7788    Continue
 8877 Continue
      Call Rys2Dg(Array(ip2D0),nT,nRys,la,lb,lc,ld,
     &            Array(ip2D1),JfGrad,
     &            JndGrd,Coora,Alpha,Beta,Gamma,Delta,nZeta,nEta,
     &            Array(ipScr),Array(ipTmp),Index,
     &            Exp_1,Exp_2,nZeta,nEta)
*-----Drop ipScr
      ip = ip - nTR
*-----Drop ipTmp
      ip = ip - nT
*
*     Assemble the gradients of the ERI's
*
      Call Assg1(Temp,PAO,nT,nRys,la,lb,lc,ld,Array(ip2D0),
     &           Array(ip2D1),JfGrad,Index,mVec)
*-----Drop ip2D1
      ip = ip - nTR*3*n2D1
*-----Drop ip2D0
      ip = ip - nTR*3*n2D0
*
*     Distribute the contributions to the molecular gradient
*
      Call Distg1(Temp,mVec,Grad,nGrad,JfGrad,JndGrd,
     &            iuvwx,lOp)
*-----Drop ipAC
*     ip = ip - nT*nPAO * 9
#ifdef _CHECK_
      If (ip.ne.1) Then
         Call WarningMessage(2,'Rysg1: ip=/=1')
         Call Abend()
      End If
#endif
*
      Return
      End
