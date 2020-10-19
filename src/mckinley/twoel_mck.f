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
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine TwoEl_mck(Coor,
     &     iAngV,iCmp,iShell,iShll,iAO,iAOst,
     &     iStb,jStb,kStb,lStb,nRys,
     &     Data1,nab,nData1,Data2,ncd,nData2,Pren,Prem,
     &     Alpha,nAlpha,iPrInc, Beta, nBeta,jPrInc,
     &     Gamma,nGamma,kPrInc,Delta,nDelta,lPrInc,
     &     Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl,
     &     Zeta,ZInv,P,rKab,nZeta,Eta,EInv,Q,rKcd,nEta,
     &     xA,xB,xG,xD,xPre,Hess,nhess,
     &     IfGrd,IndGrd,ifHss,IndHss,IfG,
     &     PSO,nPSO,Work2,nWork2,Work3,nWork3,Work4,nWork4,
     &     Aux,nAux,WorkX,nWorkX,
     &     Shijij,
     &     Dij1,Dij2,mDij,nDij,Dkl1,Dkl2,mDkl,nDkl,
     &     Dik1,Dik2,mDik,nDik,Dil1,Dil2,mDil,nDil,
     &     Djk1,Djk2,mDjk,nDjk,Djl1,Djl2,mDjl,nDjl,
     &     icmpi,Fin,nfin,Temp,nTemp,nTwo2,nFt,
     &     IndZet,IndEta,TwoHam,ipdens,Buffer,nBuffer,
     &     lgrad,ldot,n8,ltri,Dan,Din,
     &     moip,naco,rMOIN,nMOIN,new_fock)
************************************************************************
*                                                                      *
*     Input:                                                           *
*     Data1     :                                                      *
*     Data2                                                            *
*     PSO                                                              *
*     Work2                                                            *
*     Work3                                                            *
*     Work4                                                            *
*     AUX                                                              *
*     Fin    : Area for cntrctd int in sph hmn                         *
*     Temp   : Working place for F gen and n8                          *
*     TwoHam : Final results fock matrix and MO's                      *
*                                                                      *
*     Object:      To construct the first order derivatives of the AO- *
*     integrals and add them up to the MO derivatives and              *
*     the Fock matrix derivatives and contract the second              *
*     order derivatives of the AO's with the second order              *
*     density matrix.                                                  *
*                                                                      *
*     Authors: Roland Lindh, IBM Almaden Research Center, San Jose, CA *
*     March '90                                                        *
*     Anders Bernhardsson Theoretical Chemistry 95                     *
************************************************************************
*                                                                      *
*     When we are calculating the second order derivatives we need     *
*     the derivatives of the two electron integrals in three ways:     *
*                                                                      *
*     (2)                                                              *
*     1)  To calculate the static term H                               *
*     -                                                                *
*     2)  To calculate the non-zero part of <0|[E  ,H]|0>              *
*     pq                                                               *
*                                                                      *
*                                                                      *
*     3)  To calculate the derivatives of the MO orbitals with all     *
*     four indexes in the active space.                                *
*                                                                      *
*     In this implementation all contributions are calculated at       *
*     the same time.                                                   *
*                                                                      *
*     (2)                                                              *
*     The H    is calculated by contracting the second order           *
*     derivatives on the flight with the second order density matrix   *
*                                                                      *
*     (1)  (1)       (1)                                               *
*     the    F  - F    and MO     are calculated by first contracting  *
*     pq   qp                                                          *
*                                                                      *
*     the primitives and transform the integrals to  spherical         *
*     harmonics and then construct the Fock matrix as a direct SCF     *
*     The Fock matrixes is transformed to MO base and then added up to *
*     total Fock matrix                                                *
*                                                                      *
************************************************************************
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Phase_Info
      use Real_Info, only: CutInt
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External TERI1, ModU2, Cff2D
#include "Molcas.fh"
#include "ndarray.fh"
#include "real.fh"
#include "disp.fh"
#include "disp2.fh"
#include "buffer.fh"
#include "cputime.fh"
#include "print.fh"
*
      Real*8 Coor(3,4), CoorM(3,4), CoorAC(3,2),
     &     Alpha(nAlpha), Beta(nBeta), Gamma(nGamma), Delta(nDelta),
     &     xA(nZeta),xB(nZeta), xG(nEta), xD(nEta),
     &     Data1(nZeta*nDArray+nDScalar,nData1),Hess(*),
     &     Data2( nEta*nDArray+nDScalar,nData2),rKab(nZeta),rKcd(nEta),
     &     Zeta(nZeta), ZInv(nZeta),  P(nZeta,3),
     &     Eta(nEta),   EInv(nEta),    Q(nEta,3),
     &     Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj),
     &     Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl),
     &     PSO(iBasi*jBasj*kBask*lBasl,nPSO), Work2(nWork2),
     &     Work3(nWork3), Work4(nWork4),Aux(nAux),
     &     xpre(nGamma*nDelta*nAlpha*nBeta),Fin(nfin),
     &     Dij1(mDij,nDij),Dkl1(mDkl,nDkl),Dik1(mDik,nDik),
     &     Dil1(mDil,nDil),Djk1(mDjk,nDjk),Djl1(mDjl,nDjl),
     &     Dij2(mDij,nDij),Dkl2(mDkl,nDkl),Dik2(mDik,nDik),
     &     Dil2(mDil,nDil),Djk2(mDjk,nDjk),Djl2(mDjl,nDjl),
     &     WorkX(nWorkX),Temp(nTemp),TwoHam(nTwo2),
     &     Buffer(nBuffer),rMOIN(nMOIN),Din(*),Dan(* )
*
      Integer iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), iStabN(0:7),
     &     iStabM(0:7),  IndGrd(3,4,0:7), iAO(4),
     &     iCmp(4), iShell(4), iShll(4),
     &     nOp(4), iAngV(4), iAOst(4),JndGrd(3,4,0:7),icmpi(4),
     &     IndZet(nAlpha*nBeta),Indeta(nGamma*nDelta), iuvwx(4),
     &     IndHss(4,3,4,3,0:7), JndHss(4,3,4,3,0:7),
     &     Index(3,4), moip(0:7)
*
      Logical Shijij, AeqB, CeqD, AeqC, ABeqCD, ABeq, CDeq, IfGrd(3,4),
     &        JfGrd(3,4), first,IfHss(4,3,4,3),JfHss(4,3,4,3),IfG(4),
     &        ltri,Tr(4),ldot,ldot2,lgrad,n8,log,no_integrals,new_fock
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function to compute canonical index
*
      nElem(i) = (i+1)*(i+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      Call TwoEl_mck_Internal(Data1,Data2)

      Contains
      Subroutine TwoEl_mck_Internal(Data1,Data2)
      Use Iso_C_Binding
      Real*8, Target :: Data1(nZeta*nDArray+nDScalar,nData1),
     &                  Data2( nEta*nDArray+nDScalar,nData2)
      Integer, Pointer :: iData1(:), iData2(:)
      Integer :: lZeta=0, lEta=0
      Logical EQ, lEmpty
      External EQ, lEmpty
*                                                                      *
************************************************************************
*                                                                      *
*     P R O L O G
*                                                                      *
************************************************************************
*                                                                      *
      nGr=0
      ABeq = EQ(Coor(1,1),Coor(1,2))
      CDeq = EQ(Coor(1,3),Coor(1,4))
      la = iAngV(1)
      lb = iAngV(2)
      lc = iAngV(3)
      ld = iAngV(4)
      ldot2=ldot
      iSmAng=la+lb+lc+ld
      iCmpa = iCmp(1)
      jCmpb = iCmp(2)
      kCmpc = iCmp(3)
      lCmpd = iCmp(4)
      iShlla = iShll(1)
      jShllb = iShll(2)
      kShllc = iShll(3)
      lShlld = iShll(4)
      IncZet=nAlpha*jPrInc
      IncEta=nGamma*lPrInc
      LmbdT=0
      nijkl = iBasi*jBasj*kBask*lBasl
      nabcd=iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)
      mab = nElem(la)*nElem(lb)
      mcd = nElem(lc)*nElem(ld)
*
*     Scratch space for Fock Matrix construction
*
      ip=1
      ipS1=ip
      nS1=nijkl*nabcd
      ip=ip+nS1
      ipS2=ip
      nS2=max(nS1,nijkl+max(iBasi*lBasl,jBasj*lBasl,
     &     iBasi*kBask,jBasj*kBask))
      ip=ip+nS2
      ipFT=ip
      ip=ip+nFT
      ipTemp=ip
      nTe=nijkl*nabcd
      ip=ip+nTe
      If (ip-1.gt.nTemp) Then
         Write (6,*) 'TwoEl_McK: ip-1.gt.nTemp'
         Write (6,*) 'ip,nTemp=',ip,nTemp
         Call Abend()
      End If
*
      iuvwx(1) = dc(iStb)%nStab
      iuvwx(2) = dc(jStb)%nStab
      iuvwx(3) = dc(kStb)%nStab
      iuvwx(4) = dc(lStb)%nStab
*
      iffab = 9
      iffcd = 9
*                                                                      *
************************************************************************
*                                                                      *
*             - - - - - - E N D   P R O L O G - - - - - -
*                                                                      *
************************************************************************
*                                                                      *
*-----Find the Double Coset Representatives for center A and B
*
      If (nIrrep.eq.1) Then
         nDCRR=1
         iDCRR(0)=0
         LmbdR=1
      Else
         Call DCR(LmbdR,dc(iStb)%iStab,dc(iStb)%nStab,
     &                  dc(jStb)%iStab,dc(jStb)%nStab,
     &                  iDCRR,nDCRR)
      End If
      u = DBLE(dc(iStb)%nStab)
      v = DBLE(dc(jStb)%nStab)
*
*--------Find stabilizer for center A and B
*
      If (nIrrep.eq.1) Then
         lStabM=1
         iStabM(0)=0
      Else
         Call Inter(dc(iStb)%iStab,dc(iStb)%nStab,
     &              dc(jStb)%iStab,dc(jStb)%nStab,iStabM,lStabM)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Find the Double Coset Representatives for center C and D.
*
      If (nIrrep.eq.1) Then
         nDCRS=1
         iDCRS(0)=0
         LmbdS=1
      Else
         Call DCR(LmbdS,dc(kStb)%iStab,dc(kStb)%nStab,
     &                  dc(lStb)%iStab,dc(lStb)%nStab,
     &                  iDCRS,nDCRS)
      End If
      w = DBLE(dc(kStb)%nStab)
      x = DBLE(dc(lStb)%nStab)
*
*-----------Find stabilizer for center C and D
*
      If (nIrrep.eq.1) Then
         lStabN=1
         iStabN(0)=0
      Else
         Call Inter(dc(kStb)%iStab,dc(kStb)%nStab,
     &              dc(lStb)%iStab,dc(lStb)%nStab,iStabN,lStabN)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Find the Double Coset Representatives for the two charge
*     distributions.
*
      If (nIrrep.eq.1) Then
         nDCRT=1
         iDCRT(0)=0
         LmbdT=1
      Else
         Call DCR(LmbdT,iStabM,lStabM,iStabN,lStabN,iDCRT,nDCRT)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Factor due to summation over DCR
*
      If (MolWgh.eq.1) Then
         Fact = DBLE(nIrrep) / DBLE(LmbdT)
      Else If (MolWgh.eq.0) Then
         Fact = u*v*w*x / DBLE(nIrrep**3 * LmbdT)
      Else
         Fact = Sqrt(u*v*w*x)/DBLE(nIrrep*LmbdT)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      nOp(1)=NrOpr(0)
      call dcopy_(3,Coor(1,1),1,CoorM(1,1),1)
*                                                                      *
************************************************************************
*                                                                      *
*     - - - - Loop over first set
*                                                                      *
************************************************************************
*                                                                      *
      Do 100 lDCRR = 0, nDCRR-1
         nOp(2)=NrOpr(iDCRR(lDCRR))
         Call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
         AeqB = EQ(CoorM(1,1),CoorM(1,2))
*                                                                      *
************************************************************************
*                                                                      *
*     - - - - Loop over second set
*                                                                      *
************************************************************************
*                                                                      *
         Do 200 lDCRS = 0, nDCRS-1
            call dcopy_(3,Coor(1,3),1,CoorM(1,3),1)
            Call OA(iDCRS(lDCRS),Coor(1:3,4),CoorM(1:3,4))
            CeqD = EQ(Coor(1,3),CoorM(1,4))
*                                                                      *
************************************************************************
*                                                                      *
*     - - - - Loop over third set
*                                                                      *
************************************************************************
*                                                                      *
            Do 300 lDCRT = nDCRT-1, 0, -1

               nOp(3) = NrOpr(iDCRT(lDCRT))
               nOp(4) = NrOpr(iEor(iDCRT(lDCRT),iDCRS(lDCRS)))
*
               iDCRTS=iEor(iDCRT(lDCRT),iDCRS(lDCRS))
               Call OA(iDCRTS,Coor(1:3,4),CoorM(1:3,4))
               Call OA(iDCRT(lDCRT),Coor(1:3,3),CoorM(1:3,3))
*
               AeqC = EQ(CoorM(1,1),CoorM(1,3))
               ABeqCD = AeqB .and. CeqD .and. AeqC
*--------------No contribution to geometric derivatives from one-center
*              integrals
               If (ABeqCD) Go To 302
*
*--------------Find the proper centers to start of with the angular
*              momentum on. If la.eq.lb there will excist an
*              ambiguity to which center that angular momentum should
*              be accumulated on. In that case we will use A and C of
*              the order as defined by the basis functions types.
*
               If (iAngV(1).ge.iAngV(2)) Then
                  call dcopy_(3,CoorM(1,1),1,CoorAC(1,1),1)
               Else
                  call dcopy_(3,CoorM(1,2),1,CoorAC(1,1),1)
               End If
               If (iAngV(3).ge.iAngV(4)) Then
                  call dcopy_(3,CoorM(1,3),1,CoorAC(1,2),1)
               Else
                  call dcopy_(3,CoorM(1,4),1,CoorAC(1,2),1)
               End If
*
*     Calculate the desymmetrized two-electron density matrix in
*     cartisian AO base.
*
               Call Timing(dum1,Time,dum2,dum3)
               If (ldot2)
     &              Call TwoDns(iAngV,iCmp,shijij,ishll,ishell,
     &                   iAO,nOp,iBasi,jBasj,kBask,lBasl,
     &                   Aux,nAux,Work2,nWork2,Work3,nWork3,work4,
     &                   nWork4,PSO,nPSO,Fact)
*
               Call Timing(dum1,Time,dum2,dum3)
               CpuStat(nTwoDens)=CpuStat(nTwoDens)+Time
*
*----------------------------------------------------------------*
*
*     Loops to partion the primitives
*
*----------------------------------------------------------------*
               lDCR1=NrOpr(iDCRR(lDCRR))+1
               lDCR2=NrOpr(iDCRS(lDCRS))+1
               ix2 = iPhase(1,iDCRT(lDCRT))
               iy2 = iPhase(2,iDCRT(lDCRT))
               iz2 = iPhase(3,iDCRT(lDCRT))
*
               Call C_F_Pointer(C_Loc(Data1(ip_IndZ(1,nZeta),lDCR1)),
     &                         iData1,[nZeta+1])
               Call C_F_Pointer(C_Loc(Data2(ip_IndZ(1,nEta ),lDCR2)),
     &                         iData2,[nEta +1])
               nZeta_Tot=iData1(nZeta+1)
               nEta_Tot =iData2(nEta +1)
*
               no_integrals=.true.
               first=.true.
               nGr=0
               Do 400 iZeta = 1, nZeta_Tot, IncZet
                  mZeta=Min(IncZet,nZeta_Tot-iZeta+1)
*-----------------Check that subblock of contraction matrix has non-zero
*                 elements.
                  If (lEmpty(Coeff2,nBeta,nBeta,jBasj))
     &                 Go To 401
                  Do 410 iEta  = 1, nEta_Tot, IncEta
                     mEta=Min(IncEta,nEta_Tot-iEta+1)
*-----------------Check that subblock of contraction matrix has non-zero
*                 elements.
                     If (lEmpty(Coeff4,nDelta,nDelta,lBasl))
     &                    Go To 411
                     Pren = Pren + DBLE(mab*mcd*mZeta*mEta)
*-------------------------------------------------------------------------
*
*     Fix the control matrixes for derivatives
*     and try to use translation invariance as
*     efficient as possible.
*
*     OBS DETTA SKALL FLYTTAS UT UR INRE LOOPEN
*
*-------------------------------------------------------------------------
                     Call LCopy(144,IfHss,1,JfHss,1)
                     Call LCopy(12,IfGrd,1,JfGrd,1)
                     Call LCopy(4,[.true.],0,ifg,1)
                     Call LCopy(4,[.false.],0,Tr,1)
                     Call ICopy(144*nIrrep,IndHss,1,JndHss,1)
                     Call ICopy(12*nIrrep,IndGrd,1,JndGrd,1)
*
*     Delete one center that should be calculated with
*     translation invariance
*
                     call Translation(ifg,jfgrd,jfhss,tr,jndgrd,jndhss,
     &                                coorm, nirrep,indgrd,indhss)

                     if (.not.ldot)  Call LCopy(144,[.false.],0,JfHss,1)
                     if (.not.ldot)  Call iCopy(144*8,[0],0,JndHss,1)
*-------------------------------------------------------------*
*     PRE PRESCREENING                                        *
*-------------------------------------------------------------*
*
                     lZeta=mZeta
                     lEta =mEta
*
*-----------------Decontract the 2nd order density matrix
*
*     Work4->Work2  Work3:scratch
                     Call Timing(dum1,Time,dum2,dum3)
                     If (ldot2)
     &                    Call Tcrtnc_h(
     &                    Coeff1,nAlpha,iBasi,
     &                    Coeff2,nBeta,jBasj,
     &                    Coeff3,nGamma,kBask,
     &                    Coeff4,nDelta,lBasl,
     &                    Work4,mab*mcd,Work3,nWork3/2,Work2,
     &                    iData1(iZeta:iZeta+mZeta-1),mZeta,
     &                    iData2(iEta :iEta +mEta -1),mEta)
                     Call Timing(dum1,Time,dum2,dum3)
                     CPUStat(nTwoDens)=CPUStat(nTwoDens)+Time
*
*-----------------Transfer k2 data and prescreen
*
*     Work2:PAO-> Work2
*     Work3 Scratch
                     Call Timing(dum1,Time,dum2,dum3)
                     Call Screen_mck(Work2,Work3,mab*mcd,nZeta,nEta,
     &                    mZeta,mEta,lZeta,lEta,
     &                    Zeta,ZInv,P,xA,xB,rKab,
     &                    Data1(ip_Z(iZeta,nZeta),lDCR1),
     &                    iData1(iZeta:iZeta+mZeta-1),
     &                    Data1(ip_ZtMax(nZeta),ldcr1),
     &                    Data1(ip_abMax(nZeta),ldcr1),
     &                    Data1(ip_ZetaM(nZeta),ldcr1),
     &                    nAlpha,nBeta,
     &                    Eta, EInv,Q,xG,xD,rKcd,
     &                    Data2(ip_Z(iEta,nEta),lDCR2),
     &                    iData2(iEta :iEta +mEta -1),
     &                    Data2(ip_ZtMax(nEta),ldcr2),
     &                    Data2(ip_abMax(nEta),ldcr2),
     &                    Data2(ip_ZetaM(nEta),ldcr2),
     &                    nGamma,nDelta,
     &                    xpre,
     &                    1,1,1,ix2,iy2,iz2,
     &                    CutInt,
     &                    PreScr,
     &                    IndZet,IndEta,ldot2)
                     Call Timing(dum1,Time,dum2,dum3)
                     CPUStat(nScreen)=CPUStat(nScreen)+Time
*
                     Prem = Prem + DBLE(mab*mcd*lZeta*lEta)
                     If (lzeta*leta.ne.0) no_integrals=.false.
                     If (lZeta*lEta.eq.0) Go To 411
*
*-----------------Compute integral derivative and accumulate
*     contribution to the molecular gradient.
*
*     Work2:PAO
*     Work3:Work area  The PO integrals are stored in the begining
*     of Work3
*
                     Call Timing(dum1,Time,dum2,dum3)
*
                     Call Rysg2(iAngV,nRys,lZeta*lEta,
     &                    xA,xB,xG,xD,
     &                    Zeta,ZInv,lZeta,
     &                    Eta,EInv,lEta,
     &                    P,nZeta,Q,nEta,
     &                    CoorM,CoorM,CoorAC,Work3,nWork3,
     &                    TERI1,ModU2,Cff2D,
     &                    Work2,mab*mcd,
     &                    Hess,nHess,JfGrd,JndGrd,
     &                    JfHss,JndHss,nOp,iuvwx,IfG,
     &                    nGr,Index,lgrad,ldot,Tr)
                     Call Timing(dum1,Time,dum2,dum3)
                     CPUStat(nIntegrals)=CPUStat(nIntegrals)+Time

*     Work3 AO
*     Work3_3  Scratch
*     ->    Work3_2
*--------------------------------------------------------------*
*
*--------------Transform integrals ta AO base
*
*--------------------------------------------------------------*
                     ip2=nGr*mab*mcd*lZeta*lEta+1
                     Call Timing(dum1,Time,dum2,dum3)
                     Call Cntrct_mck(First,
     &                    Coeff1,nAlpha,iBasi,
     &                    Coeff2,nBeta ,jBasj,
     &                    Coeff3,nGamma,kBask,
     &                    Coeff4,nDelta,lBasl,
     &                    Work3,nGr*mab*mcd,
     &                    Work3(ip2),nwork3-ip2,
     &                    xpre,WorkX,nWorkX,
     &                    lZeta*lEta,
     &                    IndZet,nZeta,lZeta,IndEta,nEta,lEta)
 411                 Continue
 410              Continue
 401              Continue
 400           Continue
*
*
*     Mark which derivatives that should be calculated with translation
*     invarians.
*
               If (nGr.eq.0) goto 911
               Do iCNT=1,4
                  If (Tr(iCnt)) Then
                     Do iCar=1,3
                        log=.false.
                        Do iIrr=0,nIrrep-1
                           log=(log.or.indgrd(iCar,iCnt,iIrr).ne.0)
                        End Do
                        If (log) Index(iCar,iCnt)=-1
                     End Do
                  End If
               End Do
*
               If (MolWgh.eq.1) Then
                  FactNd = DBLE(nIrrep) / DBLE(LmbdT)
               Else If (MolWgh.eq.0) Then
                  FactNd = u*v*w*x / DBLE(nIrrep**3 * LmbdT)
               Else
                  factNd = sqrt(u*v*w*x)/DBLE(nirrep*lmbdt)
               End If
*
               If (FactNd.ne.One)
     &              Call DScal_(nGr*mab*mcd*nijkl,FactNd,WorkX,1)
*
*-----------------------------------------------------------------*
*
*     Transpose abcd,g,IJKL -> bcd,g,IJKL,A Work3 -> Work3_2
*
*-----------------------------------------------------------------*
*
               niag=nijkl*nElem(lb)*mcd*nGr
               Call CrSph_mck(WorkX,niag,(la+1)*(la+2)/2,
     &              RSph(ipSph(la)),la,
     &              Shells(iShlla)%Transf,
     &              Shells(iShlla)%Prjct,
     &              Work3,iCmpa)
               nw3=niag*iCmpa
               ip2=1+nw3
*-----------------------------------------------------------------*
*
*     Transpose   bcd,g,IJKL,A -> cd,g,IJKL,AB Work3_2->Work3
*
*-----------------------------------------------------------------*
               niag=nijkl*mcd*nGr*iCmpa
               nw3_2=niag*jCmpb
               If (nw3+nw3_2.gt.nWork3) Then
                  Write (6,*) '1: nw3+nw3_2.gt.nWork3'
                  Call Abend()
               End If
               Call CrSph_mck(Work3,niag,(lb+1)*(lb+2)/2,
     &              RSph(ipSph(lb)),lb,
     &              Shells(jShllb)%Transf,
     &              Shells(jShllb)%Prjct,
     &              Work3(ip2),jCmpb)
*-----------------------------------------------------------------*
*
*     Transpose  cd,g,IJKL,AB -> d,g,IJKL,ABC  Work3->Work3_2
*
*-----------------------------------------------------------------*
*
               niag=nijkl*nGr*nElem(ld)*iCmpa*jCmpb
               Call CrSph_mck(Work3(ip2),niag,(lc+1)*(lc+2)/2,
     &              RSph(ipSph(lc)),lc,
     &              Shells(kShllc)%Transf,
     &              Shells(kShllc)%Prjct,
     &              Work3,kCmpc)
               If (niag*kCmpc.gt.nw3) Then
                  Write (6,*) 'niag*kCmpc.gt.nw3'
                  Call Abend()
               End If
               nw3=niag*kCmpc
               ip2=nw3+1
*-----------------------------------------------------------------*
*
*     Transpose   d,g,IJKL,ABC -> g,IJKL,ABCD Work3_2->Work3
*
*-----------------------------------------------------------------*
               niag=nijkl*nGr*iCmpa*jCmpb*kCmpc
               nw3_2=niag*lCmpd
               If (nw3+nw3_2.gt.nWork3) Then
                  Write (6,*) '2: nw3+nw3_2.gt.nWork3'
                  Call Abend()
               End If
               Call CrSph_mck(Work3,niag,(ld+1)*(ld+2)/2,
     &              RSph(ipSph(ld)),ld,
     &              Shells(lShlld)%Transf,
     &              Shells(lShlld)%Prjct,
     &              Work3(ip2),lCmpd)
*-----------------------------------------------------------------*
*
*     Transpose g,IJKL,ABCD -> IJKL,ABCD,g Work3->Buffer
*
*-----------------------------------------------------------------*
               niag=nijkl*iCmpa*jCmpb*kCmpc*lCmpd
               Call DGetMO(Work3(ip2),nGr,nGr,niag,Fin,niag)
*
*     D E B U G  (calculates gradient from transformed integrals)
*
*
               Call Timing(dum1,Time,dum2,dum3)
               CPUStat(nTrans)=CPUStat(nTrans)+Time
*
*-----------------------------------------------------------------*
*
*     Send the integrals to clrbuffer for construction of
*
*-----------------------------------------------------------------*
*
*
               Call ClrBuf(idcrr(ldcrr),idcrs(ldcrs),
     &              idcrt(ldcrt),nGr,
     &              iStb,
     &              jStb,
     &              kStb,
     &              lStb,
     &              Shijij,iAngV,iCmpi,iCmp,
     &              iShll,iShell,iShell,
     &              iBasi,jBasj,kBask,lBasl,
     &              Dij1,Dij2,mDij,nDij,
     &              Dkl1,Dkl2,mDkl,nDkl,
     &              Dik1,Dik2,mDik,nDik,
     &              Dil1,Dil2,mDil,nDil,
     &              Djk1,Djk2,mDjk,nDjk,
     &              Djl1,Djl2,mDjl,nDjl,
     &              fin,nfin,
     &              Temp(ipFT),nFT,
     &              Temp(ipS1),nS1,Temp(ipS2),nS2,
     &              Temp(ipTemp),nTe,TwoHam,nTwo2,
     &              JndGrd,Index,iao,iaost,iuvwx,ifG,n8,ltri,
     &              moip,nAcO,rMoin,nmoin,ntemp,Buffer,
     &              coor,nOp,Din,Dan,new_fock)
 911           Continue
*
 302           Continue
*
 300        Continue
*
 200     Continue
*
 100  Continue
      Return
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nab)
         Call Unused_integer(ncd)
         Call Unused_real_array(Alpha)
         Call Unused_integer(iPrInc)
         Call Unused_real_array(Beta)
         Call Unused_real_array(Gamma)
         Call Unused_integer(kPrInc)
         Call Unused_real_array(Delta)
         Call Unused_integer(ipdens)
      End If
      End Subroutine Twoel_Mck_Internal
      End Subroutine Twoel_Mck
