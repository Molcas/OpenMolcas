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
* Copyright (C) 1993, Roland Lindh                                     *
*               1993, Per Boussard                                     *
************************************************************************
      SubRoutine M1Int(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of the M1 integrals used  *
*         ECP calculations. The operator is the nuclear attraction     *
*         operator times a s-type gaussian function.                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              mHrr                                                    *
*              DCR                                                     *
*              Rys                                                     *
*              DaXpY   (ESSL)                                          *
*              Hrr                                                     *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Alpha : exponents of bra gaussians                              *
*      nAlpha: number of primitives (exponents) of bra gaussians       *
*      Beta  : as Alpha but for ket gaussians                          *
*      nBeta : as nAlpha but for the ket gaussians                     *
*      Zeta  : sum of exponents (nAlpha x nBeta)                       *
*      ZInv  : inverse of Zeta                                         *
*      rKappa: gaussian prefactor for the products of bra and ket      *
*              gaussians.                                              *
*      P     : center of new gaussian from the products of bra and ket *
*              gaussians.                                              *
*      Final : array for computed integrals                            *
*      nZeta : nAlpha x nBeta                                          *
*      nComp : number of components in the operator (e.g. dipolmoment  *
*              operator has three components)                          *
*      la    : total angular momentum of bra gaussian                  *
*      lb    : total angular momentum of ket gaussian                  *
*      A     : center of bra gaussian                                  *
*      B     : center of ket gaussian                                  *
*      nRys  : order of Rys- or Hermite-Gauss polynomial               *
*      Array : Auxiliary memory as requested by ECPMem                 *
*      nArr  : length of Array                                         *
*      Ccoor : coordinates of the operator, zero for symmetric oper.   *
*      NOrdOp: Order of the operator                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
*             Physics, University of Stockholm, Sweden, October '93.   *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External TNAI, Fake, Cff2D, XRys2D
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), C(3),
     &       Array(nZeta*nArr), Ccoor(3), TC(3), CoorAC(3,2),
     &       Coori(3,4), Coora(3,4)
      Character*80 Label
      Integer iStabM(0:nStabM-1), iDCRT(0:7), lOper(nComp),
     &          iChO(nComp),  iAnga(4)
      Logical EQ, NoSpecial
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      Call qEnter('M1Int')
      iRout = 193
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In M1Int: A',' ',A,1,3)
         Call RecPrt(' In M1Int: RB',' ',RB,1,3)
         Call RecPrt(' In M1Int: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In M1Int: P',' ',P,nZeta,3)
         Write (6,*) ' In M1Int: la,lb=',' ',la,lb
      End If
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = 0
      iAnga(4) = 0
      mabMin = nabSz(Max(la,lb)-1)+1
      mabMax = nabSz(la+lb)
      If (EQ(A,RB)) mabMin=nabSz(la+lb-1)+1
      mAInt  = (mabMax-mabMin+1)
*
*     Find center to accumulate angular momentum on. (HRR)
*
      If (la.ge.lb) Then
       call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
       call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
*
*-----Compute FLOP's and size of work array which HRR will use.
*
      Call mHrr(la,lb,nFlop,nMem)
*
*-----Allocate Scratch for primitives and work area for HRR
*
      ip = 1
      ipAInt = ip
      k = nabSz(la+lb) - nabSz(Max(la,lb)-1)
      ip = ip + nZeta*Max(k,nMem)
      ipK = ip
      ip = ip + nZeta
      ipZ = ip
      ip = ip + nZeta
      ipZI = ip
      ip = ip + nZeta
      ipPx = ip
      ip = ip + nZeta
      ipPy = ip
      ip = ip + nZeta
      ipPz = ip
      ip = ip + nZeta
      If (ip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'M1Int: ip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr,nZeta=',nArr,nZeta
         Write (6,*) ' nMem=',nMem
         Call Abend()
      End If
      ipTmp = ip
      mArray = nArr*nZeta - ip + 1
*
      call dcopy_(nZeta*Max(k,nMem),[Zero],0,Array(ipAInt),1)
*
*-----Loop over nuclear centers.
*
      kdc = 0
      Do 100 kCnttp = 1, nCnttp
         If (.Not.ECP(kCnttp)) Go To 111
         If (nM1(kCnttp).eq.0) Go To 111
         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3)= dbsc(kCnttp)%Coor(1:3,kCnt)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+kCnt),nStab(kdc+kCnt),iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            Do 102 lDCRT = 0, nDCRT-1
               TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
               TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
               TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
               call dcopy_(3,A,1,Coora(1,1),1)
               call dcopy_(3,RB,1,Coora(1,2),1)
               call dcopy_(6,Coora(1,1),1,Coori(1,1),1)
               If (.Not.EQ(A,RB) .or. .Not.EQ(A,TC)) Then
                  Coori(1,1) = Coori(1,1)+One
*                 Coora(1,1) = Coora(1,1)+One
               End If
               call dcopy_(3,TC,1,CoorAC(1,2),1)
               call dcopy_(3,TC,1,Coori(1,3),1)
               call dcopy_(3,TC,1,Coori(1,4),1)
               call dcopy_(3,TC,1,Coora(1,3),1)
               call dcopy_(3,TC,1,Coora(1,4),1)
*
               Do 1011 iM1xp=0, nM1(kCnttp)-1
                  Gamma = Work(ipM1xp(kCnttp)+iM1xp)
*
*-----------------Modify the original basis. Observe that
*                 simplification due to A=B are not valid for the
*                 exponent index, eq. P-A=/=0.
*
                  Do 1012 iZeta = 1, nZeta
                     PTC2 = (P(iZeta,1)-TC(1))**2
     &                    + (P(iZeta,2)-TC(2))**2
     &                    + (P(iZeta,3)-TC(3))**2
                     Tmp0 = Zeta(iZeta)+Gamma
                     Tmp1 = Exp(-Zeta(iZeta)*Gamma*PTC2/Tmp0)
                     Array(ipK+iZeta-1)  = rKappa(iZeta) * Tmp1
                     Array(ipZ+iZeta-1)  = Tmp0
                     Array(ipZI+iZeta-1) = One/Tmp0
                     Array(ipPx+iZeta-1) =
     &                  (Zeta(iZeta)*P(iZeta,1)+Gamma*TC(1))/Tmp0
                     Array(ipPy+iZeta-1) =
     &                  (Zeta(iZeta)*P(iZeta,2)+Gamma*TC(2))/Tmp0
                     Array(ipPz+iZeta-1) =
     &                  (Zeta(iZeta)*P(iZeta,3)+Gamma*TC(3))/Tmp0
 1012             Continue
*
*-----------------Compute integrals with the Rys quadrature.
*
                  nT = nZeta
                  NoSpecial=.True.
                  Call Rys(iAnga,nT,Array(ipZ),Array(ipZI),nZeta,
     &                     [One],[One],1,Array(ipPx),nZeta,TC,1,
     &                     Array(ipK),[One],Coori,Coora,CoorAC,
     &                     mabmin,mabmax,0,0,Array(ipTmp),mArray,
     &                     TNAI,Fake,Cff2D,XRys2D,NoSpecial)
*
*-----------------Accumulate result for all nuclei. Take the charge on
*                 the center into account.
*
                  Factor = -Charge(kCnttp)*Work(ipM1cf(kCnttp)+iM1xp)
     &                   * Fact
                  Call DaXpY_(nZeta*mAInt,Factor,Array(ipTmp),1,
     &                       Array(ipAInt),1)
                  If (iPrint.ge.99) Then
                     Call Recprt(' [a+b,0|A|0] in Array',' ',
     &                           Array(ipTmp),nZeta,mAInt)
                     Call RecPrt(' [a+b,0|A|0] in AInt',' ',
     &                           Array(ipAInt),nZeta,mAInt)
                  End If
*
 1011          Continue
 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
 100  Continue
*
*-----Use the HRR to compute the required primitive integrals.
*
      Call HRR(la,lb,A,RB,Array(ipAInt),nZeta,nMem,ipIn)
      ii = ipAInt + ipIn - 1
*-----Move result
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,Array(ii),1,Final,1)
*
      If (iPrint.ge.99) Then
         Write (6,*) ' Result in M1Int'
         Do 150 ia = 1, nElem(la)
            Do 250 ib = 1, nElem(lb)
               Write (Label,'(A,I2,A,I2,A)')
     &               ' Final(',ia,',',ib,')'
               Call RecPrt(Label,' ',Final(1,ia,ib,1),nAlpha,nBeta)
 250        Continue
 150     Continue
      End If
*
*     Call GetMem(' Exit M1Int','LIST','REAL',iDum,iDum)
      Call QExit('M1Int')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
         Call Unused_integer(nRys)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iChO)
         Call Unused_integer(nOrdOp)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
