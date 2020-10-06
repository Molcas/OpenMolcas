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
* Copyright (C) 1990,1991,1993,1999, Roland Lindh                      *
*               1990, IBM                                              *
************************************************************************
      SubRoutine k2Loop(Coor,
     &                  iAnga,iCmpa,iShll,
     &                  iDCRR,nDCRR,Data,
     &                  Alpha,nAlpha,Beta, nBeta,
     &                  Alpha_,Beta_,
     &                  Coeff1,iBasn,Coeff2,jBasn,
     &                  Zeta,ZInv,Kappab,P,IndP,nZeta,IncZZ,Con,
     &                  Wrk,nWork2,
     &                  Cmpct,nScree,mScree,iStb,jStb,
     &                  Dij,nDij,nDCR,nHm,ijCmp,DoFock,
     &                  Scr,nScr,
     &                  Knew,Lnew,Pnew,Qnew,nNew,DoGrad,HMtrx,nHrrMtrx)
************************************************************************
*                                                                      *
* Object: to compute zeta, kappa, P, and the integrals [nm|nm] for     *
*         prescreening. This is done for all unique pairs of centers   *
*         generated from the symmetry unique centers A and B.          *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN.                              *
*             June '91, modified to compute zeta, P, kappa and inte-   *
*             grals for Schwartz inequality in a k2 loop.              *
*             Modified for direct SCF, January '93                     *
************************************************************************
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep, iOper
      use Real_Info, only: CutInt, RadMax, cdMax, EtMax
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
      External TERIS, ModU2, Cmpct, Cff2DS, Rys2D
#include "real.fh"
#include "Molcas.fh"
#include "disp.fh"
#include "print.fh"
      Real*8 Coor(3,4), CoorM(3,4), Coori(3,4), Coora(3,4), CoorAC(3,2),
     &       Alpha(nAlpha), Beta(nBeta), Dij(nDij,nDCR),
     &       Data((nZeta*(nDArray+2*ijCmp)+nDScalar+nHm),nDCRR),
     &       Zeta(nZeta), ZInv(nZeta), Kappab(nZeta), P(nZeta,3),
     &       Wrk(nWork2), Q(3), TA(3), TB(3), Scr(nScr,3),
     &       Con(nZeta), Coeff1(nAlpha,iBasn), Coeff2(nBeta,jBasn),
     &       Alpha_(nZeta),Beta_(nZeta), HMtrx(nHrrMtrx,2)
      Real*8  Knew(nNew), Lnew(nNew), Pnew(nNew*3), Qnew(nNew*3)
      Integer   iDCRR(0:7), iAnga(4), iCmpa(4), mStb(2),
     &          iShll(2), IndP(nZeta)
      Integer isave(1024)
      Logical AeqB, EQ, NoSpecial,
     &        DoFock, DoGrad
      External EQ
      Dimension Dummy(1)
*                                                                      *
************************************************************************
*                                                                      *
      Call k2loop_internal(Data)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine k2loop_internal(Data)
      Use Iso_C_Binding
      Real*8, Target :: Data((nZeta*(nDArray+2*ijCmp)+nDScalar+nHm),
     &                       nDCRR)
      Integer, Pointer :: iData(:)
      Logical TF, TstFnc
      External TstFnc
      Interface
         SubRoutine Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &                  Eta,EInv,nEta,
     &                  P,lP,Q,lQ,rKapab,rKapcd,Coori,Coora,CoorAC,
     &                  mabMin,mabMax,mcdMin,mcdMax,Array,nArray,
     &                  Tvalue,ModU2,Cff2D,Rys2D,NoSpecial)
         Integer iAnga(4), nT, nZeta, nEta, lP, lQ, mabMin, mabMax,
     &           mcdMin, mcdMax, nArray
         External Tvalue, ModU2, Cff2D, Rys2D
         Real*8 Zeta(nZeta), ZInv(nZeta), P(lP,3), rKapab(nZeta),
     &          Eta(nEta),   EInv(nEta),  Q(lQ,3), rKapcd(nEta),
     &          CoorAC(3,2), Coora(3,4), Coori(3,4), Array(nArray)
         Logical NoSpecial
         End Subroutine Rys
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function to compute canonical index
*
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
      TF(mdc,iIrrep,iComp) = TstFnc(dc(mdc)%iCoSet,
     &                              iIrrep,iComp,dc(mdc)%nStab)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 241
      iPrint = nPrint(iRout)
*     iPrint = 99
      call dcopy_(3,[One],0,Q,1)
      nData=nZeta*(nDArray+2*ijCmp)+nDScalar+nHm
      call dcopy_(nData*nDCRR,[Zero],0,Data,1)
      mStb(1) = iStb
      mStb(2) = jStb
      la = iAnga(1)
      lb = iAnga(2)
      iSmAng=la+lb+la+lb
      iCmpa_= iCmpa(1)
      jCmpb_= iCmpa(2)
      iShlla = iShll(1)
      jShllb = iShll(2)
*                                                                      *
************************************************************************
*                                                                      *
      call dcopy_(3,Coor(1,1),1,CoorM(1,1),1)
*                                                                      *
************************************************************************
*                                                                      *
      Do 100 lDCRR = 0, nDCRR-1
*
         Call ICopy(1024,nPrint,1,iSave,1)
         Call ICopy(1024,[5],0,nPrint,1)
         iR = iDCRR(lDCRR)
*
         Call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
         AeqB = EQ(CoorM(1,1),CoorM(1,2))
*        Branch out if integrals are zero by symmetry.
         If (AeqB .and. Mod(iSmAng,2).eq.1) Go To 100
         call dcopy_(6,CoorM(1,1),1,CoorM(1,3),1)
         If (iPrint.ge.99) Call RecPrt(' Actual centers',
     &                                  ' ',CoorM,3,4)
*                                                                      *
************************************************************************
*                                                                      *
*        Compute zeta, P and kappa.
*
*        No triangulatization applied at this level
         Call DoZeta(Alpha,nAlpha,Beta,nBeta,
     &               CoorM(1,1),CoorM(1,2),
     &               P,
     &               Zeta,
     &               Kappab,
     &               ZInv,
     &               Alpha_,
     &               Beta_,
     &               IndP)
*                                                                      *
************************************************************************
*                                                                      *
*        Generate transformation matrix from intermediate integrals
*        to final angular composition.
*
         mabMin=nabSz(Max(la,lb)-1)+1
         If (EQ(CoorM(1,1),CoorM(1,2))) mabMin = nabSz(la+lb-1)+1
         mabMax=nabSz(la+lb)
         ne=(mabMax-mabMin+1)
         Do iIrrep = 0, nIrrep-1
            i13_=ip_HrrMtrx(nZeta)+(iIrrep*nHm)/nIrrep
            Call OA(iOper(iIrrep),CoorM(1:3,1),TA)
            Call OA(iOper(iIrrep),CoorM(1:3,2),TB)
            Call HrrMtrx(Data(i13_,lDCRR+1),
     &                   ne,la,lb,TA,TB,
     &                   Shells(iShlla)%Transf,RSph(ipSph(la)),iCmpa_,
     &                   Shells(jShllb)%Transf,RSph(ipSph(lb)),jCmpb_)
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*        Compute primitive integrals to be used in the prescreening
*        by the Schwartz inequality.
*
         call dcopy_(12,CoorM(1,1),1,Coora(1,1),1)
         call dcopy_(12,CoorM(1,1),1,Coori(1,1),1)
*
*        Compute actual size of [a0|c0] block
*
         mcdMin=mabMin
         mcdMax=mabMax
         mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
*
*        Find the proper centers to start of with the angular
*        momentum on. If la.eq.lb there will excist an
*        ambiguity to which center that angular momentum should
*        be accumulated on. In that case we will use A and C of
*        the order as defined by the basis functions types.
*
         If (iAnga(1).ge.iAnga(2)) Then
            call dcopy_(3,Coora(1,1),1,CoorAC(1,1),1)
         Else
            call dcopy_(3,Coora(1,2),1,CoorAC(1,1),1)
         End If
         call dcopy_(3,CoorAC(1,1),1,CoorAC(1,2),1)
*
*        Compute [a0|c0], ijkl,a,c
*
         Jnd = 0
         nScree = nScree + nZeta
         Do iZeta = 1, nZeta, IncZZ
            mZeta = Min(nZeta-iZeta+1,IncZZ)
*
            nT = mZeta*1
            NoSpecial=.True.
            Call Rys(iAnga,nT,
     &               Zeta(iZeta),ZInv(iZeta),mZeta,[One],[One],1,
     &               P(iZeta,1),nZeta,Q,1,Kappab(iZeta),[One],
     &               Coori,Coora,CoorAC,
     &               mabMin,mabMax,mcdMin,mcdMax,
     &               Wrk,nWork2,TERIS,ModU2,Cff2DS,
     &               Rys2D,NoSpecial)
            If (iPrint.ge.59)
     &         Call RecPrt(' In k2Loop: ijkl,[a0|c0]',' ',Wrk,
     &                               mZeta,mabcd)
*
*---------- Apply a transpose prior to Tnsctl to fake the action
*           of Cntrct.
*
            iW3=1+mZeta*mabcd
            Call DGeTMO(Wrk,mZeta,mZeta,mabcd,Wrk(iW3),mabcd)
            call dcopy_(mabcd*mZeta,Wrk(iW3),1,Wrk,1)
            Call TnsCtl(Wrk,nWork2,Coora,
     &                  mabcd,mZeta,mabMax,mabMin,mabMax,mabMin,
     &                  Data(ip_HrrMtrx(nZeta),lDCRR+1),
     &                  Data(ip_HrrMtrx(nZeta),lDCRR+1),
     &                  la,lb,la,lb,
     &                  iCmpa_,jCmpb_,iCmpa_,jCmpb_,
     &                  iShlla,jShllb,iShlla,jShllb,i_Int)
            If (i_Int.eq.1) Then
               iW2=1
               iW3=1+mZeta*(iCmpa_*jCmpb_)**2
            Else
               iW2=i_Int
               iW3=1
            End If
*                                                                      *
************************************************************************
*                                                                      *
         Call C_F_Pointer(C_Loc(Data(ip_IndZ(1,nZeta),lDCRR+1)),iData,
     &                    [nAlpha*nBeta+1])
*                                                                      *
************************************************************************
*                                                                      *
*-----------Store data in core
*
            Call Cmpct(Wrk(iW2),iCmpa_,jCmpb_,nZeta,mZeta,
     &                 Zeta(iZeta),Kappab(iZeta),
     &                 P(iZeta,1),IndP(iZeta),Con,
     &                 Data(ip_Z    (1,nZeta),lDCRR+1),
     &                 Data(ip_Kappa(1,nZeta),lDCRR+1),
     &                 Data(ip_Pcoor(1,nZeta),lDCRR+1),
     &                 iData,iZeta-1,Jnd,
     &                 Data(ip_ZInv (1,nZeta),lDCRR+1),CutInt,RadMax,
     &                 cdMax,EtMax,AeqB,
     &                 Data(ip_ab   (1,nZeta),lDCRR+1),
     &                 Data(ip_abCon(1,nZeta),lDCRR+1),
     &                 Alpha_(iZeta),
     &                 Data(ip_Alpha(1,nZeta,1),lDCRR+1),
     &                 Beta_(iZeta),
     &                 Data(ip_Beta (1,nZeta,2),lDCRR+1))
*
         End Do ! iZeta
         mScree = mScree + Jnd
*                                                                      *
************************************************************************
*                                                                      *
*        Estimate the largest contracted integral.
*
         Data(ip_EstI(nZeta),lDCRR+1) =
     &                      EstI(Data(ip_Z(1,nZeta),lDCRR+1),
     &                           Data(ip_Kappa(1,nZeta),lDCRR+1),
     &                           nAlpha,nBeta,
     &                           Coeff1,iBasn,Coeff2,jBasn,
     &                           Data(ip_ab   (1,nZeta),lDCRR+1),
     &                           iCmpa_*jCmpb_,
     &                           Wrk,nWork2,
     &                           iData)
         Nullify(iData)
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
         iOffZ = nDij-nZeta-1
         ZtMax=One
         abMax=Zero
         ZtMaxD=One
         abMaxD=Zero
         Do iZeta = 0, Jnd-1
            tmp = Data(ip_abCon(iZeta+1,nZeta),lDCRR+1)
            If (abMax.lt.tmp) Then
               abMax = tmp
               ZtMax = Data(ip_Z    (iZeta+1,nZeta),lDCRR+1)
            End If
            If (DoFock) Then
               tmp = Data(ip_ab(iZeta+1,nZeta),lDCRR+1)
     &             * Dij(iOffZ+iZeta,lDCRR+1)
               If (abMaxD.lt.tmp) Then
                 abMaxD = tmp
                 ZtMaxD = Data(ip_Z(iZeta+1,nZeta),lDCRR+1)
               End If
            Else
                 ZtMaxD=-One
                 abMaxD=Zero
            End If
         End Do
         Data(ip_ZtMax(nZeta),lDCRR+1) = ZtMax
         Data(ip_abMax (nZeta),lDCRR+1) = abMax
         Data(ip_ZtMaxD(nZeta),lDCRR+1) = ZtMaxD
         Data(ip_abMaxD(nZeta),lDCRR+1) = abMaxD
         Call ICopy(1024,iSave,1,nPrint,1)
*                                                                      *
************************************************************************
*                                                                      *
*------- Compute data for gradient evaluation
*
         mZeta=Jnd
         If (DoGrad.and.mZeta.gt.0) Then
*
*---------- Compute primitive integrals to be used in the prescreening
*           by the Cauchy-Schwarz inequality.
*
            Do iZeta = 1, mZeta, IncZZ
               lZeta = Min(mZeta-iZeta+1,IncZZ)
               Call SchInt(CoorM,iAnga,iCmpa,lZeta,
     &                     Data(ip_Z    (iZeta,nZeta),lDCRR+1),
     &                     Data(ip_ZInv (iZeta,nZeta),lDCRR+1),
     &                     Data(ip_Kappa(iZeta,nZeta),lDCRR+1),
     &                     Data(ip_PCoor(iZeta,nZeta),lDCRR+1),
     &                     Data(ip_Kappa(iZeta,nZeta),lDCRR+1),
     &                     Data(ip_PCoor(iZeta,nZeta),lDCRR+1),
     &                     nZeta,Wrk,nWork2,HMtrx,
     &                     nHrrMtrx,iShlla,jShllb,i_Int)
               Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                     Data(ip_abG(nZeta,nHm)+iZeta-1,lDCRR+1),
     &                     Data(ip_Kappa(iZeta,nZeta),lDCRR+1),.True.,
     &                     Data(ip_Z   (iZeta,nZeta),lDCRR+1),nZeta,
     &                     Dummy)
            End Do
*
*---------- Second order numerical differentiation. The gradients are
*           restricted to only those with respect to symmetrical
*           displacements. The symmetric three point formula is used in
*           the numerical procedure.
*
            iIrrep=0
            Delta = 1.0D-03
            iOff_g=ip_abG(nZeta,nHm)+ijCmp*nZeta
            Call FZero(Data(iOff_g,lDCRR+1),nZeta*ijCmp)
            Scr(1:nZeta*ijCmp,:)=Zero
*
*---------- Loop over center A and B.
*
            Do iCnt = 1, 2
*
               nDisp=IndDsp(mStb(iCnt),iIrrep)
               Do iComp = 1, 3
                  iCmp=2**(iComp-1)
                  If (TF(mStb(iCnt),iIrrep,iCmp) .and.
     &                Direct(nDisp+1)) Then
                     nDisp = nDisp + 1
                     temp = CoorM(iComp,iCnt)
*
                     CoorM(iComp,iCnt  ) = temp + Delta
                     CoorM(iComp,iCnt+2) = temp + Delta
                     Call NewPK(CoorM(1,1),CoorM(1,2),
     &                          Pnew,mZeta,nZeta,
     &                          Knew,
     &                          Data(ip_Alpha(1,nZeta,1),lDCRR+1),
     &                          Data(ip_Beta (1,nZeta,2),lDCRR+1))
                     Do iZeta = 1, mZeta, IncZZ
                        lZeta = Min(mZeta-iZeta+1,IncZZ)
                        Call SchInt(CoorM,
     &                              iAnga,iCmpa,lZeta,
     &                              Data(ip_Z   (iZeta,nZeta),lDCRR+1),
     &                              Data(ip_ZInv(iZeta,nZeta),lDCRR+1),
     &                              Knew(iZeta),
     &                              Pnew(iZeta),
     &                              Knew(iZeta),
     &                              Pnew(iZeta),
     &                              nZeta,Wrk,nWork2,
     &                              HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
                        Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                              Scr(iZeta,1),
     &                              Knew(iZeta),
     &                              .False.,
     &                              Data(ip_Z    (iZeta,nZeta),lDCRR+1),
     &                              nZeta,
     &                              Knew(iZeta))
                     End Do
*
                     CoorM(iComp,iCnt  ) = temp - Delta
                     CoorM(iComp,iCnt+2) = temp - Delta
                     Call NewPK(CoorM(1,1),CoorM(1,2),
     &                          Qnew,mZeta,nZeta,
     &                          Lnew,
     &                          Data(ip_Alpha(1,nZeta,1),lDCRR+1),
     &                          Data(ip_Beta (1,nZeta,2),lDCRR+1))
                     Do iZeta = 1, mZeta, IncZZ
                        lZeta = Min(mZeta-iZeta+1,IncZZ)
                        Call SchInt(CoorM,
     &                              iAnga,iCmpa,lZeta,
     &                              Data(ip_Z   (iZeta,nZeta),lDCRR+1),
     &                              Data(ip_ZInv(iZeta,nZeta),lDCRR+1),
     &                              Lnew(iZeta),
     &                              Qnew(iZeta),
     &                              Lnew(iZeta),
     &                              Qnew(iZeta),
     &                              nZeta,Wrk,nWork2,
     &                              HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
                        Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                              Scr(iZeta,2),
     &                              Lnew(iZeta),
     &                              .False.,
     &                              Data(ip_Z   (iZeta,nZeta),lDCRR+1),
     &                              nZeta,
     &                              Lnew(iZeta))
                     End Do
*
                     Call DaXpY_(nZeta*ijCmp, One,Scr(1,2),1,
     &                                           Scr(1,1),1)
*
                     CoorM(iComp,iCnt  ) = temp + Delta
                     CoorM(iComp,iCnt+2) = temp - Delta
                     Do iZeta = 1, mZeta, IncZZ
                        lZeta = Min(mZeta-iZeta+1,IncZZ)
                        Call SchInt(CoorM,
     &                              iAnga,iCmpa,lZeta,
     &                              Data(ip_Z   (iZeta,nZeta),lDCRR+1),
     &                              Data(ip_ZInv(iZeta,nZeta),lDCRR+1),
     &                              Knew(iZeta),
     &                              Pnew(iZeta),
     &                              Lnew(iZeta),
     &                              Qnew(iZeta),
     &                              nZeta,Wrk,nWork2,
     &                              HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
                        Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                              Scr(iZeta,3),
     &                              Knew(iZeta),
     &                              .False.,
     &                              Data(ip_Z    (iZeta,nZeta),lDCRR+1),
     &                              nZeta,
     &                              Lnew(iZeta))
                     End Do
*
                     Call DaXpY_(nZeta*ijCmp,-One,Scr(1,3),1,
     &                                           Scr(1,1),1)
*
                     CoorM(iComp,iCnt  ) = temp - Delta
                     CoorM(iComp,iCnt+2) = temp + Delta
                     Do iZeta = 1, mZeta, IncZZ
                        lZeta = Min(mZeta-iZeta+1,IncZZ)
                        Call SchInt(CoorM,
     &                              iAnga,iCmpa,lZeta,
     &                              Data(ip_Z   (iZeta,nZeta),lDCRR+1),
     &                              Data(ip_ZInv(iZeta,nZeta),lDCRR+1),
     &                              Lnew(iZeta),
     &                              Qnew(iZeta),
     &                              Knew(iZeta),
     &                              Pnew(iZeta),
     &                              nZeta,Wrk,nWork2,
     &                              HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
                        Call PckInt(Wrk(i_Int),lZeta,ijCmp,
     &                              Scr(iZeta,3),
     &                              Lnew(iZeta),
     &                              .False.,
     &                              Data(ip_Z    (iZeta,nZeta),lDCRR+1),
     &                              nZeta,
     &                              Knew(iZeta))
                     End Do
*
                     Call DaXpY_(nZeta*ijCmp,-One,Scr(1,3),1,
     &                                           Scr(1,1),1)
*
                     Call DScal_(nZeta*ijCmp,One/(Four*Delta**2),
     &                                           Scr(1,1),1)
                     Call AbsAdd(nZeta*ijCmp,    Scr(1,1),1,
     &                                   Data(iOff_g,lDCRR+1),1)
*
                     CoorM(iComp,iCnt  ) = temp
                     CoorM(iComp,iCnt+2) = temp
                  End If
               End Do
*
            End Do
         End If       ! DoGrad
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUG_
#ifdef _DEBUG_
         Write (6,*)
         Write (6,*) 'lDCRR=',lDCRR
         Call WrCheck('Zeta ',Data(ip_Z    (1,nZeta),  lDCRR+1),nZeta)
         Call WrCheck('Kappa',Data(ip_Kappa(1,nZeta),  lDCRR+1),nZeta)
         Call WrCheck('P    ',Data(ip_PCoor(1,nZeta),  lDCRR+1),nZeta*3)
         Call WrCheck('xA   ',Data(ip_Alpha(1,nZeta,1),lDCRR+1),nZeta)
         Call WrCheck('xB   ',Data(ip_Beta (1,nZeta,2),lDCRR+1),nZeta)
         Call WrCheck('ZInv ',Data(ip_ZInv (1,nZeta),  lDCRR+1),nZeta)
         If (DoGrad) Then
            Call WrCheck('ab   ',
     &         Data(ip_abG  (nZeta,nHm),  lDCRR+1),nZeta*ijCmp)
            iOff_g=ip_abG(nZeta,nHm)+ijCmp*nZeta
            Call WrCheck('abG  ',
     &         Data(iOff_g,               lDCRR+1),nZeta*ijCmp)
         End If
         Write (6,*)
         Write (6,*) ' ERI(Max)=',Data(  ip_EstI(nZeta),lDCRR+1)
         Write (6,*) ' ZtMax   =',Data( ip_ZtMax(nZeta),lDCRR+1)
         Write (6,*) ' abMax   =',Data(ip_abMax (nZeta),lDCRR+1)
         Write (6,*) ' ZtMaxD  =',Data(ip_ZtMaxD(nZeta),lDCRR+1)
         Write (6,*) ' abMaxD  =',Data(ip_abMaxD(nZeta),lDCRR+1)
         Call WrCheck(' HrrMtrx',
     &        Data(ip_HrrMtrx(nZeta),lDCRR+1),
     &        ne*iCmpa_*jCmpb_)
#endif
 100  Continue ! lDCRR
*
      Return
      End Subroutine k2loop_internal
*
      End
