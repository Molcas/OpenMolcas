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
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Rysg2(iAnga,nRys,nT,
     &                 Alpha,Beta,Gamma,Delta,
     &                 Zeta,ZInv,nZeta,Eta,EInv,nEta,
     &                 P,lP,Q,lQ,Coori,Coora,CoorAC,
     &                 Array,nArray,
     &                 Tvalue,ModU2,Cff2D,
     &                 PAO,nPAO,Hess,nHess,IfGrd,IndGrd,
     &                 IfHss,IndHss,nOp,iuvwx,IfG,
     &                 mvec,Index_Out,lGrad,lHess,Tr)
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
*              Assg2                                                   *
*              Distg2                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Modified to 1st order derivatives October '91            *
*             Modified to 2nd order derviatives Mars '95 by            *
*             Anders Bernhardsson Theoretical Chemistry,               *
*             University of Lund                                       *
************************************************************************
*   @parameter iAnga Angular momenta for each center
*   @parameter nRys  Order of rys polynomia
*   @parameter nT    Number of alpha-beta-gamma-delta multiplies
*   @parameter Alpha Explonents on 1st center
*   @parameter Beta Explonents on 2nd center
*   @parameter Gamma Explonents on 3rd center
*   @parameter Delta Explonents on 4ht center
*   @parameter Zeta  Alpha*Beta
*   @parameter Zeta  Zeta invers
*   @parameter nZeta Alpha Beta multiplies
*   @parameter Eta  Gamma*Delta
*   @parameter Eta  Eta invers
*   @parameter nEta gamma delta multiplies
*   @parameter P
*   @parameter lP  Length of p
*   @parameter Q
*   @parameter lQ  Length of Q
*   @parameter Coori Coordinates of center just used to check AeqB CeqD etc
*   @parameter Coora Coordinates of center used in hor. recursion
*   @parameter CoorAC Coordinates of center <max(la,lb),max(lc,ld)>
*   @parameter Array Scratch and output area for 1st derivatives
*   @parameter nArray Size of scratch area
*   @parameter PAO Density
*   @parameter nPAO lenth of density
*   @parameter Hess Output area for hessian (added)
*   @parameter nHess Size of hessian
*   @parameter ifgrad true for all 1st derivatives that are needed
*   @parameter indgrad  index in gradient on which integrals should be added
*   @parameter ifhss true for all 2nd derivatives that are needed
*   @parameter indhss  index in hess on which contracted integrals should be added
*   @parameter nop oerator number for the operator that generates center
*   @parameter iuvwx number of stabilazors
*   @parameter ifg true for all centers on which derivatives should be calculated
*   @parameter index_out  index where first derivatives are stored (out)
*   @parameter lhess true if 2nd derivatives should be calculated
*   @parameter lgrad true if 1st derivatives should be calculated
*   @parameter tr  true for all centers on which should be calculated via translation invarians
      use vRys_RW
      use Symmetry_Info, only: iOper
      use Real_Info, only: ChiI2
      use Temporary_Parameters, only: IsChi
      Implicit Real*8 (A-H,O-Z)
      External Tvalue, ModU2, Cff2D
      External Exp_1, Exp_2
#include "itmax.fh"
#include "info.fh"
#include "notab.fh"
#include "real.fh"
      Real*8 Zeta(nZeta), ZInv(nZeta), P(lP,3),
     &       Eta(nEta),   EInv(nEta),  Q(lQ,3),
     &       Alpha(nZeta), Beta(nZeta), Gamma(nEta), Delta(nEta),
     &       CoorAC(3,2), Coora(3,4), Coori(3,4), Array(nArray),
     &       PAO(nT,nPAO), Hess(nHess)
*    &       rkab(nZeta),rdcd(neta)
      Integer iAnga(4), IndGrd(3,4,0:nirrep-1), Index1(3,4),
     &          nOp(4), iuvwx(4),   JndGrd(3,4,0:7), lOp(4),
     &          IndHss(4,3,4,3,0:7),Index2(3,4,4),
     &          Index3(3,3),Index4(2,6,3),ng(3),nh(3),Index_Out(3,4)
      Logical AeqB, CeqD, AeqC, EQ, IfGrd(3,4),
     &          JfGrd(3,4),IfG(4),KfGrd(3,4),
     &          ifhss(4,3,4,3),lgrad,lhess,Tr(4)
      nElem(i) = (i+1)*(i+2)/2
*

      Call LCopy(12,[.false.],0,JfGrd,1)
      Call LCopy(12,[.false.],0,KfGrd,1)
      lOp(1) = iOper(nOp(1))
      lOp(2) = iOper(nOp(2))
      lOp(3) = iOper(nOp(3))
      lOp(4) = iOper(nOp(4))
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
         If (IfHss(1,i,1,i)) lla=2
         If (IfGrd(i,1)) lla=max(lla,1)
         If (IfHss(2,i,2,i)) llb=2
         If (IfGrd(i,2)) llb=max(llb,1)
         If (IfHss(3,i,3,i)) llc=2
         If (IfGrd(i,3)) llc=max(llc,1)
         If (IfHss(4,i,4,i)) lld=2
         If (IfGrd(i,4)) lld=max(lld,1)
 10   Continue
      lab = Max(lla,llb)
      lcd = Max(llc,lld)
      nabMax = la + lb + lab
      ncdMax = lc + ld + lcd
*
*     Allocate memory for the integral gradients.
*
*
      ip = 1

      ipFinal=ip
      MemFinal=9*nt*nElem(la)*nElem(lb)*nElem(lc)*nElem(ld)
      ip=ip+MemFinal
*
*     Allocate memory for the 2D-integrals.
*
      ip2D0 = ip
      n2D0 = Max( (nabMax+1)*(ncdMax+1),
     &            (la+3)*(lb+3)*(ncdMax+1),
     &            (la+3)*(lb+3)*(lc+3)*(ld+3))
      ip = ip + n2D0*3*nT*nRys
*
*     Allocate memory for the 2nd order derivatives of the 2D-integrals.
*
      ip2D1 = ip
      n2D1 = Max( (nabMax+1)*(ncdMax+1),
     &            (la+3)*(lb+3)*(ncdMax+1),
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
      If (ip-1.gt.nArray) Then
         Call WarningMessage(2,'Rysg2: ip-1.gt.nArray (pos. 1)')
         Write (6,*) 'ip,nArray=',ip,nArray
         Call Abend()
      End If
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
         If (ip-1.gt.nArray) Then
            Call WarningMessage(2,'Rysg2: ip-1.gt.nArray (pos. 2)')
            Write (6,*) 'ip,nArray=',ip,nArray
            Call Abend()
         End If
*
         Call RtsWgh(Array(ipTv),nT,Array(ipU2),Array(ipWgh),nRys)
      Else
         If (ip-1.gt.nArray) Then
            Call WarningMessage(2,'Rysg2: ip-1.gt.nArray (pos. 3)')
            Write (6,*) 'ip,nArray=',ip,nArray
            Call Abend()
         End If

*        Make sure rys11/her11 is called  (la+1)
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
     &        Array(ipU2),Array(ipPAQP),Array(ipQCPQ),
     &        Array(ipB10),Array(ipB00),labMax,Array(ipB01))
*-----Drop ipU2
      ip = ip - nT*nRys
*     Let go of Zeta, ZInv, Eta, and EInv
      ip = ip - nT*4
*     Let go of P and Q
      ip = ip - 6*nT
*
*     -------------------------------
*     Compute the intermediate 2D-integrals from the roots and weights.
*
      Call Rs2Dmm(Array(ip2D0),nT,nRys,nabMax,
     &            ncdMax,Array(ipPAQP),Array(ipQCPQ),
     &            Array(ipB10),Max(nabMax-1,0),
     &            Array(ipB00),labMax,
     &            Array(ipB01),Max(ncdMax-1,0),
     &            la,lb,lc,ld,IfHss,IfGrd)
*
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
      Call HrrCtl_mck(Array(ip2D0),n2D0,Array(ip2D1),n2D1,
     &            la,lb,lc,ld,nabmax,ncdmax,
     &            nTR,Coora(1,1),Coora(1,2),Coora(1,3),Coora(1,4),
     &            IfHss,IfGrd,nt,nrys)
*
*     Compute the gradients of the 2D-integrals. Copy some information
*     which will be modified. This has to be done in order to facilitate
*     partioning.
*
      ip2D2= ip
      n2D2 =  (la+1)*(lb+1)*(lc+1)*(ld+1)*18
      ip = ip + n2D2*nT*nRys
      ipScr = ip
      ip = ip + nT*nRys
      ipTmp = ip
      ip = ip + nrys*nT
      ipScr2=ip
      ip=ip+nT*nRys
      If (ip-1.gt.nArray) Then
         Call WarningMessage(2,'Rysg2: ip-1.gt.nArray (pos. 4)')
         Write (6,*) 'ip,nArray=',ip,nArray
         Call Abend()
      End If
      Call ICopy(12*nirrep,IndGrd,1,JndGrd,1)
      Do 8877 i = 1, 3
         Do 7788 j = 1, 4
            JfGrd(i,j) = IfGrd(i,j)
 7788    Continue
 8877 Continue
      Do iCent=1,4
       Do jCar=1,3
        Do kCent=1,i
         If (iCent.eq.kCent) then
           iStop=3
         Else
          iStop=jCar-1
         End If
         Do lCar=1,istop
          if (jCar.ne.lCar) Then
             if (ifhss(iCent,jCar,kCent,lCar)) Then
                kfgrd(jCar,iCent)=.true.
                kfgrd(lCar,kCent)=.true.
             End If
          End If
         End Do
        End Do
       End Do
      End Do
      Do jCent=1,4
       Do iCar=1,3
        if (kfgrd(iCar,jCent).or.Ifgrd(iCar,jCent)) Then
         kfgrd(iCar,jCent)=.true.
        End If
       End Do
      End Do
*
      Call Rs2Dgh(Array(ip2D0),nT,nRys,la,lb,lc,ld,
     &            Array(ip2D1),Array(ip2D2),
     &            IfHss,IndHss,KfGrd,JndGrd,IfG,
     &            Coora,Alpha,Beta,Gamma,Delta,nZeta,nEta,
     &            Array(ipScr),Array(ipScr2),Array(ipTmp),
     &            Index1,Index2,Index3,Index4,ng,nh,
     &            Exp_1,Exp_2,nZeta,nEta,nIrrep,Tr)
*-----Drop ipScr
      ip = ip - nTR
*-----Drop ipScr2
      ip = ip - nTR
*-----Drop ipTmp
      ip = ip - nTR
*
*     G2
*
      ipg2 = ip
      ip = ip + 78
      If (ip-1.gt.nArray) Then
         Call WarningMessage(2,'Rysg2: ip-1.gt.nArray (pos. 5)')
         Write (6,*) 'ip,nArray=',ip,nArray
         Call Abend()
      End If
*
      If (lGrad) Then
       Do iCe=1,4
        Do iCa=1,3
         kfGrd(iCa,iCe)=(kfgrd(iCa,iCe).and.ifgrd(iCa,iCe))
        End Do
       End Do
       Call Assg1_mck(Array,nT,nRys,la,lb,lc,ld,Array(ip2D0),
     &                Array(ip2D1),kfGrd,Index1,mVec_,Index_out)
       mVec = mVec_

      End If
*
      If (lHess) Then
         Call Assg2(Array(ipg2),nT,nRys,la,lb,lc,ld,Array(ip2D0),
     &              Array(ip2D1),Array(ip2D2),IfHss,
     &              Index3,Index4,ng,nh,PAO)
*
*        Distribute the contributions to the molecular static hessian
*
         Call Distg2(Array(ipg2),Hess,nHess,JndGrd,
     &            IfHss,IndHss,iuvwx,lOp,nop,Tr,IfG)
      End If
*
*-----Drop ipg2
      ip = ip - 78
*-----Drop ip2D1
      ip = ip - nTR*3*n2D1
*-----Drop ip2D0
      ip = ip - nTR*3*n2D0
*-----Drop ip2D2
      ip=ip-n2D2*nT*nRys
*
      If (ip.ne.1+MemFinal) Then
         Call WarningMessage(2,'Rysg2: ip.ne.1+MemFinal (pos. 5)')
         Write (6,*) 'ip,MemFinal=',ip,MemFinal
         Call Abend()
      End If
      Call lCopy(12,kfgrd,1,ifgrd,1)
*
      Return
      End
