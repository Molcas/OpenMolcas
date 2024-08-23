!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,1991,1993, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine TwoEl_NoSym(iS_,jS_,kS_,lS_,                           &
     &           Coor,                                                  &
     &           iAnga,iCmp,iShell,iShll,iAO,iAOst,                     &
     &           NoInts,iStabs,                                         &
     &           nAlpha,iPrInc, nBeta,jPrInc,                           &
     &           nGamma,kPrInc,nDelta,lPrInc,                           &
     &           nData1,nData2,                                         &
     &           k2Data1,k2Data2,                                       &
     &           IJeqKL,kOp,                                            &
     &           Dij,mDij,mDCRij,Dkl,mDkl,mDCRkl,Dik,mDik,mDCRik,       &
     &           Dil,mDil,mDCRil,Djk,mDjk,mDCRjk,Djl,mDjl,mDCRjl,       &
     &           Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl,   &
     &           FckTmp,nFT,nZeta,nEta,                                 &
     &           SOInt,nSOInt,Wrk,nWork2,                               &
     &           Shijij,Aux,nAux)
!***********************************************************************
!                                                                      *
! Object: to generate the SO integrals for four fixed centers and      *
!         fixed basis set types.                                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!          Roland Lindh, Dept. of Theoretical Chemistry, University of *
!          Lund, SWEDEN. Modified to use Schwartz inequality for pre-  *
!          screening, July 1991.                                       *
!          Modified for direct SCF, January '93                        *
!***********************************************************************
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Gateway_Info, only: ThrInt, CutInt
      use Symmetry_Info, only: nIrrep
      use Int_Options, only: DoIntegrals, DoFock, FckNoClmb, FckNoExch
      use Int_Options, only: ExFac, Thize, W2Disc, IntOnly=>PreSch
      use Int_Options, only: Disc_Mx, Disc, Quad_ijkl
      use k2_arrays, only: TwoHam=>pFq, Dens=>pDq
      use Breit, only: nComp
      use Constants
#ifdef _DEBUGPRINT_
#endif
      use k2_structure, only: k2_type
      Implicit None
#include "twoswi.fh"
      Integer iS_,jS_,kS_,lS_,                                          &
     &        nAlpha,iPrInc, nBeta,jPrInc,                              &
     &        nGamma,kPrInc,nDelta,lPrInc,                              &
     &        nData1,nData2,                                            &
     &        mDij,mDCRij,mDkl,mDCRkl,mDik,mDCRik,                      &
     &        mDil,mDCRil,mDjk,mDCRjk,mDjl,mDCRjl,                      &
     &        iBasi,jBasj,kBask,lBasl,                                  &
     &        nFT,nZeta,nEta,                                           &
     &        nSOInt,nWork2,                                            &
     &        nAux
      Real*8 Coor(3,4)
      Integer iAnga(4), iCmp(4), iShell(4), iShll(4), iAO(4), iAOst(4)
      Logical NoInts
      Integer iStabs(4)
      Type(k2_type) k2data1(nData1), k2Data2(nData2)
      Logical IJeqKL
      Integer kOp(4)
      Real*8 Dij(mDij,mDCRij),Dkl(mDkl,mDCRkl),Dik(mDik,mDCRik),        &
     &       Dil(mDil,mDCRil),Djk(mDjk,mDCRjk),Djl(mDjl,mDCRjl),        &
     &       Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj),                 &
     &       Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl),                &
     &       FckTmp(nFT)
      Real*8 SOInt(iBasi*jBasj*kBask*lBasl,nSOInt),                     &
     &       Wrk(nWork2)
      Logical Shijij
      Real*8 Aux(nAux)

!Local variables
      Real*8 CoorAC(3,2), QInd(2)
      Integer iWR(2)
      Logical NoPInts, AeqB, CeqD, AeqC, ABeqCD,                        &
     &        EQ, Do_TnsCtl, IeqK,JeqL,                                 &
     &        Pij, Pkl, Pijkl, Pik, Pjl,                                &
     &        lEmpty, Prescreen_On_Int_Only, DoCoul, DoExch,            &
     &        Scrij, Scrkl, Scrik, Scril, Scrjk, Scrjl,                 &
     &        Batch_On_Disk, DoAOBatch, All_Spherical
      Logical ::  Copy=.True., NoCopy=.False.
#include "SysDef.fh"
      External EQ, lEmpty
      Integer iStb, jStb, kStb, lStb

      Integer la, lb, lc, ld, ISMAng, nab, ncd, nijkl, nInts, ipAOInt,  &
     &        iW3, iW4, kInts, mInts, mabMin, mabMax, mcdMin, mcdMax,   &
     &        mabcd, IncZet, IncEta, mWork2, nZeta_Tot, nEta_Tot, kabcd,&
     &        ipAOInt_, mZeta, mEta, iOpt, i_Int, nByte,                &
     &        iPer, nabcd, iW4_, iZeta, iEta
      Real*8 RST_Triplet, vijkl, vij, vkl, vik, vil, vjk, vjl, q4

!                                                                      *
!***********************************************************************
!                                                                      *
!     Declaration of statement functions to compute canonical index
!
      Integer ixyz, nabSz
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
      iStb = iStabs(1)
      jStb = iStabs(2)
      kStb = iStabs(3)
      lStb = iStabs(4)
!
      All_Spherical=Shells(iShll(1))%Prjct.and.                         &
     &              Shells(iShll(2))%Prjct.and.                         &
     &              Shells(iShll(3))%Prjct.and.                         &
     &              Shells(iShll(4))%Prjct
!
#ifdef _DEBUGPRINT_
      Call RecPrt('Coeff1',' ',Coeff1,nAlpha,iBasi)
      Call RecPrt('Coeff2',' ',Coeff2,nBeta,jBasj)
      Call RecPrt('Coeff3',' ',Coeff3,nGamma,kBask)
      Call RecPrt('Coeff4',' ',Coeff4,nDelta,lBasl)
#endif
!
      RST_triplet=One
      QInd(2)=RST_triplet
      kOp(1)=0
      kOp(2)=0
      kOp(3)=0
      kOp(4)=0
!
      la = iAnga(1)
      lb = iAnga(2)
      lc = iAnga(3)
      ld = iAnga(4)
      iSmAng=la+lb+lc+ld
!
!switch (to generate better start orbitals...)
      AeqB = EQ(Coor(1,1),Coor(1,2))
      If (NDDO .AND. .NOT.AeqB) Go To 99
      CeqD = EQ(Coor(1,3),Coor(1,4))
      If (NDDO .AND. .NOT.CeqD) Go To 99
!switch
!
      AeqC = EQ(Coor(1,1),Coor(1,3))
      ABeqCD = AeqB .and. CeqD .and. AeqC
      If (ABeqCD .and. Mod(iSmAng,2).eq.1) Go To 99
!     For Spherical Gaussians, batches like
!     (DS|SS), (FP|SS) and (FS|PS) vanish as well
      If (ABeqCD .and. All_Spherical .and.                              &
     &    2*Max(la,lb,lc,ld).gt.iSmAng) Go To 99
!
      nab = iCmp(1)*iCmp(2)
      ncd = iCmp(3)*iCmp(4)
      nijkl = iBasi*jBasj*kBask*lBasl*nComp
      nabcd = nab*ncd
      nInts =nijkl*nabcd
      ipAOInt=1
      iW3=1+nInts
      iW4=1
!
      vijkl = k2Data1(1)%abMax * k2Data2(1)%abMax
!
      Batch_On_Disk = (vijkl.gt.Thize) .and.                            &
     &       (Disc+DBLE(nInts+2+2/RtoI).le.Disc_Mx)
!
      Prescreen_On_Int_Only = IntOnly
      If (DoIntegrals) Prescreen_On_Int_Only=.True.
      If (Batch_On_Disk) Prescreen_On_Int_Only = .True.
!
      If (DoFock) Then
         vij = Dij(mDij,1)
         vkl = Dkl(mDkl,1)
!        Coulomb contributions
         Scrkl = vij*vijkl.ge.ThrInt
         Scrij = vkl*vijkl.ge.ThrInt
         DoCoul  = Scrij .or. Scrkl
!
!        Exchange contributions
         vik = Dik(mDik,1)/Four
         vil = Dil(mDil,1)/Four
         vjk = Djk(mDjk,1)/Four
         vjl = Djl(mDjl,1)/Four
         Scrjl = vik*vijkl.ge.ThrInt
         Scrjk = vil*vijkl.ge.ThrInt
         Scril = vjk*vijkl.ge.ThrInt
         Scrik = vjl*vijkl.ge.ThrInt
         DoExch = Scrjl .or. Scrjk .or. Scril .or. Scrik
         If (FckNoClmb) DoCoul=.False.
         If (FckNoExch) DoExch=.False.
      Else
         DoCoul=.False.
         DoExch=.False.
      End If
      DoAOBatch=(DoIntegrals.and.vijkl.gt.CutInt).or.                   &
     &          (DoFock.and.(DoCoul.or.DoExch)) .or.                    &
     &          (Batch_On_Disk.and.W2Disc)
!
!-----Branch out if crude estimate indicates no contributions!
!
      If (.Not.DoAOBatch) Then
         If (.Not.Batch_On_Disk) Then
            Go To 99
         Else If (Batch_On_Disk.and..Not.W2Disc) Then
 1111       Continue
            Call iRBuf(iWR,2,Copy)
            Call dRBuf(QInd,2,Copy)
            Call Store_QLast(QInd)
            kInts=iWR(1)
            mInts=iWR(2)
            If (QInd(1).eq.Quad_ijkl) Then
               If (kInts.ne.nInts) Then
                  Call WarningMessage(2,                                &
     &                        'Twoel: kInts.ne.nInts!')
                  Write (6,*) 'Twoel: kInts,mInts,nInts=',              &
     &                                kInts,mInts,nInts
                  Write (6,*) 'Index,1:',QInd(1),Quad_ijkl
                  Call Abend()
               End If
               If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,NoCopy)
               Disc = Disc + DBLE(2/RtoI + 2 + mInts)
               Go To 99
            Else If (QInd(1).lt.Quad_ijkl) Then
               If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,NoCopy)
               Disc = Disc + DBLE(2/RtoI + 2 + mInts)
               Go To 1111
            Else
               Call WarningMessage(2,'Twoel: batch is lost!')
               Write (6,*) 'Index,1:',QInd(1),QInd(2),                  &
     &                      Quad_ijkl,RST_triplet
               Call Abend()
            End If
         End If
      End If
!
!
!-----Branch point for partial integral storage
!
      If (Batch_On_Disk.and..Not.W2Disc) Go To 6767
!                                                                      *
!***********************************************************************
!                                                                      *
!     Here if the AO batch will be computed!
!
!-----Compute actual size of the {a0|c0} block
!
      mabMin=nabSz(Max(la,lb)-1)+1
      If (EQ(Coor(1,1),Coor(1,2))) mabMin = nabSz(la+lb-1)+1
      mabMax=nabSz(la+lb)
      mcdMin=nabSz(Max(lc,ld)-1)+1
      If (EQ(Coor(1,3),Coor(1,4))) mcdMin = nabSz(lc+ld-1)+1
      mcdMax=nabSz(lc+ld)
      mabcd=(mabMax-mabMin+1)*(mcdMax-mcdMin+1)
!
!-----Find the proper centers to start of with the angular
!     momentum on. If la.eq.lb there will exist an
!     ambiguity to which center that angular momentum should
!     be accumulated on. In that case we will use A and C of
!     the order as defined by the basis functions types.
!
      If (iAnga(1).ge.iAnga(2)) Then
         call dcopy_(3,Coor(1,1),1,CoorAC(1,1),1)
      Else
         call dcopy_(3,Coor(1,2),1,CoorAC(1,1),1)
      End If
      If (iAnga(3).ge.iAnga(4)) Then
         call dcopy_(3,Coor(1,3),1,CoorAC(1,2),1)
      Else
         call dcopy_(3,Coor(1,4),1,CoorAC(1,2),1)
      End If
!
!-----Set flags if triangularization will be used
!
      IeqK = EQ(Coor(1,1),Coor(1,3))
      JeqL = EQ(Coor(1,2),Coor(1,4))
      IJeqKL = IeqK .and. JeqL
!
!-----Loops to partion the primitives
!
      IncZet=nAlpha*jPrInc
      IncEta=nGamma*lPrInc
      If (nZeta.ne.IncZet.or.nEta.ne.IncEta) Then
         mWork2 = nWork2 - nijkl*mabcd
         ipAOInt=1+nijkl*mabcd
      Else
         mWork2 = nWork2
         ipAOInt=1
      End If
!
      nZeta_Tot=k2Data1(1)%IndZ(nZeta+1)
      nEta_Tot =k2Data2(1)%IndZ(nEta +1)
#ifdef _DEBUGPRINT_
      Write (6,*) 'nZeta_Tot, IncZet=',nZeta_Tot, IncZet
      Write (6,*) 'nEta_Tot,  IncEta=',nEta_Tot,  IncEta
#endif
!
      kabcd=0
      Do_TnsCtl=.False.
      NoInts=.True.
      NoPInts = .True.
      ipAOInt_=ipAOInt
      iW4_=iW4
!
      Do iZeta = 1, nZeta_Tot, IncZet
         mZeta=Min(IncZet,nZeta_Tot-iZeta+1)
         If (lEmpty(Coeff2,nBeta,nBeta,jBasj)) Cycle
!
         Do iEta  = 1, nEta_Tot,  IncEta
            mEta=Min(IncEta,nEta_Tot-iEta+1)
            If (lEmpty(Coeff4,nDelta,nDelta,lBasl)) Cycle
!
            Call DrvRys(iZeta,iEta,nZeta,nEta,mZeta,mEta,               &
     &                  nZeta_Tot,nEta_Tot,                             &
     &                  k2data1(1),                                     &
     &                  k2data2(1),                                     &
     &                  nAlpha,nBeta,nGamma,nDelta,                     &
     &                  1,1,1,1,1,1,ThrInt,CutInt,                      &
     &                  vij,vkl,vik,vil,vjk,vjl,                        &
     &                  Prescreen_On_Int_Only,NoInts,iAnga,             &
     &                  Coor,CoorAC,                                    &
     &                  mabMin,mabMax,mcdMin,mcdMax,nijkl/nComp,        &
     &                  nabcd,mabcd,Wrk,ipAOInt_,iW4_,                  &
     &                  nWork2,mWork2,                                  &
     &                  k2Data1(1)%HrrMtrx(:,1),                        &
     &                  k2Data2(1)%HrrMtrx(:,1),                        &
     &                  la,lb,lc,ld,                                    &
     &                  iCmp,iShll,NoPInts,                             &
     &                  Dij(1,1),mDij,Dkl(1,1),mDkl,Do_TnsCtl,kabcd,    &
     &                  Coeff1,iBasi,Coeff2,jBasj,                      &
     &                  Coeff3,kBask,Coeff4,lBasl)
!
         End Do
      End Do
!
      If (NoPInts) Then
         If (W2Disc) Then
            If (Batch_On_Disk) Then
               iOpt=0
               mInts=0
!
               iWR(1)=nInts
               iWR(2)=mInts
               Call iWBuf(iWR,2)
               QInd(1)=Quad_ijkl
               Call dWBuf(QInd,2)
               Call Store_QLast(QInd)
!
               Disc = Disc + DBLE(2/RtoI + 2 + mInts)
            End If
         End If
         Go To 99
      End If
!
!-----Apply the transfer equation and transform the spherical
!     harmonic gaussian.
!
      If (Do_TnsCtl) Then
         Call TnsCtl(Wrk(iW4),nWork2,                                   &
     &               nijkl,mabMax,mabMin,mcdMax,mcdMin,                 &
     &               k2Data1(1)%HrrMtrx(:,1),                           &
     &               k2Data2(1)%HrrMtrx(:,1),                           &
     &               la,lb,lc,ld,                                       &
     &               iCmp(1),iCmp(2),iCmp(3),iCmp(4),                   &
     &               iShll(1),iShll(2),iShll(3),iShll(4),i_Int)
         ipAOInt=i_Int
         If (i_Int.eq.1) Then
            iW3=1+nijkl*nabcd
         Else
            iW3=1
         End If
      Else
!
!--------Undo the late Cntrct
!
         call dcopy_(nijkl*nabcd,Wrk(ipAOInt),1,Wrk(iW3),1)
         Call DGeTMO(Wrk(iW3),nabcd,nabcd,nijkl,Wrk(ipAOInt),nijkl)
!
      End If
#ifdef _DEBUGPRINT_
      Call RecPrt('(AB|CD)',' ',Wrk(ipAOInt),nijkl/nComp,               &
     &            nComp*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4))
#endif
!
!-----Branch point for partial integral storage
!
      If (Batch_On_Disk.and.W2Disc) Then
!
!--------Write integrals to current position on disc.
!
         iOpt=0
         Call PkR8(iOpt,nInts,nByte,Wrk(ipAOInt),Wrk(iW3))
         mInts=(nByte+RtoB-1)/RtoB
!
         iWR(1)=nInts
         iWR(2)=mInts
!        Write (*,*) 'nInts,mInts=',nInts,mInts
         Call iWBuf(iWR,2)
         QInd(1)=Quad_ijkl
         Call dWBuf(QInd,2)
         Call Store_QLast(QInd)
         Call dWBuf(Wrk(iW3),mInts)
!
         Disc = Disc + DBLE(2/RtoI + 2 + mInts)
!
      End If
!
 6767 Continue
      If (Batch_On_Disk.and..Not.W2Disc) Then
 1112    Continue
         Call iRBuf(iWR,2,Copy)
         Call dRBuf(QInd,2,Copy)
         Call Store_QLast(QInd)
         kInts=iWR(1)
         mInts=iWR(2)
         If (QInd(1).eq.Quad_ijkl) Then
            If (kInts.ne.nInts) Then
               Call WarningMessage(2,                                   &
     &                     'Twoel: kInts.ne.nInts!')
               Write (6,*) 'Twoel: kInts,mInts,nInts=',                 &
     &                             kInts,mInts,nInts
               Write (6,*) 'Index,1:',QInd(1),Quad_ijkl
               Call Abend()
            End If
            If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,Copy)
            Disc = Disc + DBLE(2/RtoI + 2 + mInts)
            If (mInts.eq.0) Go To 99
         Else If (QInd(1).lt.Quad_ijkl) Then
            If (mInts.gt.0) Call dRBuf(Wrk(iW3),mInts,NoCopy)
            Disc = Disc + DBLE(2/RtoI + 2 + mInts)
            Go To 1112
         Else
            Call WarningMessage(2,'Twoel: batch is lost!')
            Write (6,*) 'Index,1:',QInd(1),QInd(2),                     &
     &                   Quad_ijkl,RST_triplet
            Call Abend()
         End If
!
         iOpt=0
         Call UpkR8(iOpt,nInts,nByte,Wrk(iW3),Wrk(ipAOInt))
      End If
!
!-----Accumulate contributions directly to the Fock matrix.
!
      If (DoFock)                                                       &
     &Call FckAcc_NoSymq(iCmp(1),iCmp(2),iCmp(3),iCmp(4),               &
     &                  Shijij, iShell, nijkl,                          &
     &                  Wrk(ipAOInt),TwoHam,Dens,Size(TwoHam),          &
     &                  iAO,iAOst,                                      &
     &                  iBasi,jBasj,kBask,lBasl,DoCoul,DoExch,          &
     &                  vij,vkl,vik,vil,vjk,vjl,ExFac)
!
      If (DoIntegrals) Then
         If (ipAOInt.ne.1) Then
            call dcopy_(nijkl*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4),          &
     &                  Wrk(ipAOInt),1,Wrk(1),1)
            ipAOInt=1
         ENd If
         iPer = 1
         Pij= iS_.eq.jS_
         Pkl= kS_.eq.lS_
         Pik= iS_.eq.kS_
         Pjl= jS_.eq.lS_
         Pijkl= Pij .and. Pkl .and. Pik .and. Pjl
         If (Pij)   iPer = iPer*2
         If (Pkl)   iPer = iPer*2
         If (Pijkl) iPer = iPer*2
         q4 = DBLE(8)/DBLE(iPer)
         If (nIrrep.eq.1) q4 = One
         If (q4.ne.One) Call DScal_(nijkl*iCmp(1)*iCmp(2)               &
     &                            *iCmp(3)*iCmp(4),q4,Wrk(ipAOInt),1)
      End If
  99  Continue
!
      Return
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iStb)
         Call Unused_integer(jStb)
         Call Unused_integer(kStb)
         Call Unused_integer(lStb)
         Call Unused_integer(iPrInc)
         Call Unused_integer(kPrInc)
         Call Unused_integer(nData1)
         Call Unused_integer(nData2)
         Call Unused_real_array(FckTmp)
         Call Unused_real_array(SoInt)
         Call Unused_real_array(Aux)
      End If
!
      End SubRoutine TwoEl_NoSym
