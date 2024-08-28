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
subroutine TwoEl_NoSym(iS_,jS_,kS_,lS_,Coor,iAnga,iCmp,iShell,iShll,iAO,iAOst,NoInts,iStabs,nAlpha,iPrInc,nBeta,jPrInc,nGamma, &
                       kPrInc,nDelta,lPrInc,nData1,nData2,k2Data1,k2Data2,IJeqKL,kOp,Dij,mDij,mDCRij,Dkl,mDkl,mDCRkl,Dik,mDik, &
                       mDCRik,Dil,mDil,mDCRil,Djk,mDjk,mDCRjk,Djl,mDjl,mDCRjl,Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl, &
                       FckTmp,nFT,nZeta,nEta,SOInt,nSOInt,Wrk,nWork2,Shijij,Aux,nAux)
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

use Basis_Info, only: Shells
use Gateway_Info, only: ThrInt, CutInt
use Symmetry_Info, only: nIrrep
use Int_Options, only: DoIntegrals, DoFock, FckNoClmb, FckNoExch
use Int_Options, only: ExFac, Thize, W2Disc, IntOnly => PreSch
use Int_Options, only: Disc_Mx, Disc, Quad_ijkl
use k2_arrays, only: TwoHam => pFq, Dens => pDq
use Breit, only: nComp
use Constants, only: One, Four, Eight
use Definitions, only: wp, u6
use k2_structure, only: k2_type

implicit none
#include "twoswi.fh"
integer iS_, jS_, kS_, lS_, nAlpha, iPrInc, nBeta, jPrInc, nGamma, kPrInc, nDelta, lPrInc, nData1, nData2, mDij, mDCRij, mDkl, &
        mDCRkl, mDik, mDCRik, mDil, mDCRil, mDjk, mDCRjk, mDjl, mDCRjl, iBasi, jBasj, kBask, lBasl, nFT, nZeta, nEta, nSOInt, &
        nWork2, nAux
real*8 Coor(3,4)
integer iAnga(4), iCmp(4), iShell(4), iShll(4), iAO(4), iAOst(4)
logical NoInts
integer iStabs(4)
type(k2_type) k2data1(nData1), k2Data2(nData2)
logical IJeqKL
integer kOp(4)
real*8 Dij(mDij,mDCRij), Dkl(mDkl,mDCRkl), Dik(mDik,mDCRik), Dil(mDil,mDCRil), Djk(mDjk,mDCRjk), Djl(mDjl,mDCRjl), &
       Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj), Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl), FckTmp(nFT)
real*8 SOInt(iBasi*jBasj*kBask*lBasl,nSOInt), Wrk(nWork2)
logical Shijij
real*8 Aux(nAux)
real*8 CoorAC(3,2), QInd(2)
integer iWR(2)
logical NoPInts, AeqB, CeqD, AeqC, ABeqCD, EQ, Do_TnsCtl, IeqK, JeqL, Pij, Pkl, Pijkl, Pik, Pjl, lEmpty, Prescreen_On_Int_Only, &
        DoCoul, DoExch, Scrij, Scrkl, Scrik, Scril, Scrjk, Scrjl, Batch_On_Disk, DoAOBatch, All_Spherical
logical :: Copy = .true., NoCopy = .false.
#include "SysDef.fh"
external EQ, lEmpty
integer iStb, jStb, kStb, lStb
integer la, lb, lc, ld, ISMAng, nab, ncd, nijkl, nInts, ipAOInt, iW3, iW4, kInts, mInts, mabMin, mabMax, mcdMin, mcdMax, mabcd, &
        IncZet, IncEta, mWork2, nZeta_Tot, nEta_Tot, kabcd, ipAOInt_, mZeta, mEta, iOpt, i_Int, nByte, iPer, nabcd, iW4_, iZeta, &
        iEta
real*8 RST_Triplet, vijkl, vij, vkl, vik, vil, vjk, vjl, q4
! Statement function
integer ixyz, nabSz
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

iStb = iStabs(1)
jStb = iStabs(2)
kStb = iStabs(3)
lStb = iStabs(4)

All_Spherical = Shells(iShll(1))%Prjct .and. Shells(iShll(2))%Prjct .and. Shells(iShll(3))%Prjct .and. Shells(iShll(4))%Prjct

#ifdef _DEBUGPRINT_
call RecPrt('Coeff1',' ',Coeff1,nAlpha,iBasi)
call RecPrt('Coeff2',' ',Coeff2,nBeta,jBasj)
call RecPrt('Coeff3',' ',Coeff3,nGamma,kBask)
call RecPrt('Coeff4',' ',Coeff4,nDelta,lBasl)
#endif

RST_triplet = One
QInd(2) = RST_triplet
kOp(1) = 0
kOp(2) = 0
kOp(3) = 0
kOp(4) = 0

la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
iSmAng = la+lb+lc+ld

! switch (to generate better start orbitals...)
AeqB = EQ(Coor(1,1),Coor(1,2))
if (NDDO .and. (.not. AeqB)) Go To 99
CeqD = EQ(Coor(1,3),Coor(1,4))
if (NDDO .and. (.not. CeqD)) Go To 99
! switch

AeqC = EQ(Coor(1,1),Coor(1,3))
ABeqCD = AeqB .and. CeqD .and. AeqC
if (ABeqCD .and. (mod(iSmAng,2) == 1)) Go To 99
! For Spherical Gaussians, batches like
!(DS|SS), (FP|SS) and (FS|PS) vanish as well
if (ABeqCD .and. All_Spherical .and. (2*max(la,lb,lc,ld) > iSmAng)) Go To 99

nab = iCmp(1)*iCmp(2)
ncd = iCmp(3)*iCmp(4)
nijkl = iBasi*jBasj*kBask*lBasl*nComp
nabcd = nab*ncd
nInts = nijkl*nabcd
ipAOInt = 1
iW3 = 1+nInts
iW4 = 1

vijkl = k2Data1(1)%abMax*k2Data2(1)%abMax

Batch_On_Disk = (vijkl > Thize) .and. (Disc+real(nInts+2+2/RtoI,kind=wp) <= Disc_Mx)

Prescreen_On_Int_Only = IntOnly
if (DoIntegrals) Prescreen_On_Int_Only = .true.
if (Batch_On_Disk) Prescreen_On_Int_Only = .true.

if (DoFock) then
  vij = Dij(mDij,1)
  vkl = Dkl(mDkl,1)
  ! Coulomb contributions
  Scrkl = vij*vijkl >= ThrInt
  Scrij = vkl*vijkl >= ThrInt
  DoCoul = Scrij .or. Scrkl

  ! Exchange contributions
  vik = Dik(mDik,1)/Four
  vil = Dil(mDil,1)/Four
  vjk = Djk(mDjk,1)/Four
  vjl = Djl(mDjl,1)/Four
  Scrjl = vik*vijkl >= ThrInt
  Scrjk = vil*vijkl >= ThrInt
  Scril = vjk*vijkl >= ThrInt
  Scrik = vjl*vijkl >= ThrInt
  DoExch = Scrjl .or. Scrjk .or. Scril .or. Scrik
  if (FckNoClmb) DoCoul = .false.
  if (FckNoExch) DoExch = .false.
else
  DoCoul = .false.
  DoExch = .false.
end if
DoAOBatch = (DoIntegrals .and. (vijkl > CutInt)) .or. (DoFock .and. (DoCoul .or. DoExch)) .or. (Batch_On_Disk .and. W2Disc)

! Branch out if crude estimate indicates no contributions!

if (.not. DoAOBatch) then
  if (.not. Batch_On_Disk) Go To 99
  if (.not. W2Disc) then
1111 continue
    call iRBuf(iWR,2,Copy)
    call dRBuf(QInd,2,Copy)
    call Store_QLast(QInd)
    kInts = iWR(1)
    mInts = iWR(2)
    if (QInd(1) == Quad_ijkl) then
      if (kInts /= nInts) then
        call WarningMessage(2,'Twoel: kInts /= nInts!')
        write(u6,*) 'Twoel: kInts,mInts,nInts=',kInts,mInts,nInts
        write(u6,*) 'Index,1:',QInd(1),Quad_ijkl
        call Abend()
      end if
      if (mInts > 0) call dRBuf(Wrk(iW3),mInts,NoCopy)
      Disc = Disc+real(2/RtoI+2+mInts,kind=wp)
      Go To 99
    else if (QInd(1) < Quad_ijkl) then
      if (mInts > 0) call dRBuf(Wrk(iW3),mInts,NoCopy)
      Disc = Disc+real(2/RtoI+2+mInts,kind=wp)
      Go To 1111
    else
      call WarningMessage(2,'Twoel: batch is lost!')
      write(u6,*) 'Index,1:',QInd(1),QInd(2),Quad_ijkl,RST_triplet
      call Abend()
    end if
  end if
end if

! Branch point for partial integral storage

if (Batch_On_Disk .and. (.not. W2Disc)) Go To 6767
!                                                                      *
!***********************************************************************
!                                                                      *
! Here if the AO batch will be computed!

! Compute actual size of the {a0|c0} block

mabMin = nabSz(max(la,lb)-1)+1
if (EQ(Coor(1,1),Coor(1,2))) mabMin = nabSz(la+lb-1)+1
mabMax = nabSz(la+lb)
mcdMin = nabSz(max(lc,ld)-1)+1
if (EQ(Coor(1,3),Coor(1,4))) mcdMin = nabSz(lc+ld-1)+1
mcdMax = nabSz(lc+ld)
mabcd = (mabMax-mabMin+1)*(mcdMax-mcdMin+1)

! Find the proper centers to start of with the angular
! momentum on. If la == lb there will exist an
! ambiguity to which center that angular momentum should
! be accumulated on. In that case we will use A and C of
! the order as defined by the basis functions types.

if (iAnga(1) >= iAnga(2)) then
  call dcopy_(3,Coor(1,1),1,CoorAC(1,1),1)
else
  call dcopy_(3,Coor(1,2),1,CoorAC(1,1),1)
end if
if (iAnga(3) >= iAnga(4)) then
  call dcopy_(3,Coor(1,3),1,CoorAC(1,2),1)
else
  call dcopy_(3,Coor(1,4),1,CoorAC(1,2),1)
end if

! Set flags if triangularization will be used

IeqK = EQ(Coor(1,1),Coor(1,3))
JeqL = EQ(Coor(1,2),Coor(1,4))
IJeqKL = IeqK .and. JeqL

! Loops to partion the primitives

IncZet = nAlpha*jPrInc
IncEta = nGamma*lPrInc
if ((nZeta /= IncZet) .or. (nEta /= IncEta)) then
  mWork2 = nWork2-nijkl*mabcd
  ipAOInt = 1+nijkl*mabcd
else
  mWork2 = nWork2
  ipAOInt = 1
end if

nZeta_Tot = k2Data1(1)%IndZ(nZeta+1)
nEta_Tot = k2Data2(1)%IndZ(nEta+1)
#ifdef _DEBUGPRINT_
write(u6,*) 'nZeta_Tot, IncZet=',nZeta_Tot,IncZet
write(u6,*) 'nEta_Tot,  IncEta=',nEta_Tot,IncEta
#endif

kabcd = 0
Do_TnsCtl = .false.
NoInts = .true.
NoPInts = .true.
ipAOInt_ = ipAOInt
iW4_ = iW4

do iZeta=1,nZeta_Tot,IncZet
  mZeta = min(IncZet,nZeta_Tot-iZeta+1)
  if (lEmpty(Coeff2,nBeta,nBeta,jBasj)) cycle

  do iEta=1,nEta_Tot,IncEta
    mEta = min(IncEta,nEta_Tot-iEta+1)
    if (lEmpty(Coeff4,nDelta,nDelta,lBasl)) cycle

    call DrvRys(iZeta,iEta,nZeta,nEta,mZeta,mEta,nZeta_Tot,nEta_Tot,k2data1(1),k2data2(1),nAlpha,nBeta,nGamma,nDelta,1,1,1,1,1,1, &
                ThrInt,CutInt,vij,vkl,vik,vil,vjk,vjl,Prescreen_On_Int_Only,NoInts,iAnga,Coor,CoorAC,mabMin,mabMax,mcdMin,mcdMax, &
                nijkl/nComp,nabcd,mabcd,Wrk,ipAOInt_,iW4_,nWork2,mWork2,k2Data1(1)%HrrMtrx(:,1),k2Data2(1)%HrrMtrx(:,1),la,lb,lc, &
                ld,iCmp,iShll,NoPInts,Dij(1,1),mDij,Dkl(1,1),mDkl,Do_TnsCtl,kabcd,Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4, &
                lBasl)

  end do
end do

if (NoPInts) then
  if (W2Disc) then
    if (Batch_On_Disk) then
      iOpt = 0
      mInts = 0

      iWR(1) = nInts
      iWR(2) = mInts
      call iWBuf(iWR,2)
      QInd(1) = Quad_ijkl
      call dWBuf(QInd,2)
      call Store_QLast(QInd)

      Disc = Disc+real(2/RtoI+2+mInts,kind=wp)
    end if
  end if
  Go To 99
end if

! Apply the transfer equation and transform the spherical
! harmonic gaussian.

if (Do_TnsCtl) then
  call TnsCtl(Wrk(iW4),nWork2,nijkl,mabMax,mabMin,mcdMax,mcdMin,k2Data1(1)%HrrMtrx(:,1),k2Data2(1)%HrrMtrx(:,1),la,lb,lc,ld, &
              iCmp(1),iCmp(2),iCmp(3),iCmp(4),iShll(1),iShll(2),iShll(3),iShll(4),i_Int)
  ipAOInt = i_Int
  if (i_Int == 1) then
    iW3 = 1+nijkl*nabcd
  else
    iW3 = 1
  end if
else

  ! Undo the late Cntrct

  call dcopy_(nijkl*nabcd,Wrk(ipAOInt),1,Wrk(iW3),1)
  call DGeTMO(Wrk(iW3),nabcd,nabcd,nijkl,Wrk(ipAOInt),nijkl)

end if
#ifdef _DEBUGPRINT_
call RecPrt('(AB|CD)',' ',Wrk(ipAOInt),nijkl/nComp,nComp*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4))
#endif

! Branch point for partial integral storage

if (Batch_On_Disk .and. W2Disc) then

  ! Write integrals to current position on disc.

  iOpt = 0
  call PkR8(iOpt,nInts,nByte,Wrk(ipAOInt),Wrk(iW3))
  mInts = (nByte+RtoB-1)/RtoB

  iWR(1) = nInts
  iWR(2) = mInts
  !write(u6,*) 'nInts,mInts=',nInts,mInts
  call iWBuf(iWR,2)
  QInd(1) = Quad_ijkl
  call dWBuf(QInd,2)
  call Store_QLast(QInd)
  call dWBuf(Wrk(iW3),mInts)

  Disc = Disc+real(2/RtoI+2+mInts,kind=wp)

end if

6767 continue
if (Batch_On_Disk .and. (.not. W2Disc)) then
1112 continue
  call iRBuf(iWR,2,Copy)
  call dRBuf(QInd,2,Copy)
  call Store_QLast(QInd)
  kInts = iWR(1)
  mInts = iWR(2)
  if (QInd(1) == Quad_ijkl) then
    if (kInts /= nInts) then
      call WarningMessage(2,'Twoel: kInts /= nInts!')
      write(u6,*) 'Twoel: kInts,mInts,nInts=',kInts,mInts,nInts
      write(u6,*) 'Index,1:',QInd(1),Quad_ijkl
      call Abend()
    end if
    if (mInts > 0) call dRBuf(Wrk(iW3),mInts,Copy)
    Disc = Disc+real(2/RtoI+2+mInts,kind=wp)
    if (mInts == 0) Go To 99
  else if (QInd(1) < Quad_ijkl) then
    if (mInts > 0) call dRBuf(Wrk(iW3),mInts,NoCopy)
    Disc = Disc+real(2/RtoI+2+mInts,kind=wp)
    Go To 1112
  else
    call WarningMessage(2,'Twoel: batch is lost!')
    write(u6,*) 'Index,1:',QInd(1),QInd(2),Quad_ijkl,RST_triplet
    call Abend()
  end if

  iOpt = 0
  call UpkR8(iOpt,nInts,nByte,Wrk(iW3),Wrk(ipAOInt))
end if

! Accumulate contributions directly to the Fock matrix.

if (DoFock) call FckAcc_NoSymq(iCmp(1),iCmp(2),iCmp(3),iCmp(4),Shijij,iShell,nijkl,Wrk(ipAOInt),TwoHam,Dens,size(TwoHam),iAO, &
                               iAOst,iBasi,jBasj,kBask,lBasl,DoCoul,DoExch,vij,vkl,vik,vil,vjk,vjl,ExFac)

if (DoIntegrals) then
  if (ipAOInt /= 1) then
    call dcopy_(nijkl*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4),Wrk(ipAOInt),1,Wrk(1),1)
    ipAOInt = 1
  end if
  iPer = 1
  Pij = iS_ == jS_
  Pkl = kS_ == lS_
  Pik = iS_ == kS_
  Pjl = jS_ == lS_
  Pijkl = Pij .and. Pkl .and. Pik .and. Pjl
  if (Pij) iPer = iPer*2
  if (Pkl) iPer = iPer*2
  if (Pijkl) iPer = iPer*2
  q4 = Eight/real(iPer,kind=wp)
  if (nIrrep == 1) q4 = One
  if (q4 /= One) call DScal_(nijkl*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4),q4,Wrk(ipAOInt),1)
end if
99 continue

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(iStb)
  call Unused_integer(jStb)
  call Unused_integer(kStb)
  call Unused_integer(lStb)
  call Unused_integer(iPrInc)
  call Unused_integer(kPrInc)
  call Unused_integer(nData1)
  call Unused_integer(nData2)
  call Unused_real_array(FckTmp)
  call Unused_real_array(SoInt)
  call Unused_real_array(Aux)
end if

end subroutine TwoEl_NoSym
