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
subroutine TwoEl_NoSym( &
#                      define _CALLING_
#                      include "twoel_interface.fh"
                      )
!***********************************************************************
!                                                                      *
! Object: to generate the SO integrals for four fixed centers and      *
!         fixed basis set types.                                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!          Roland Lindh, Dept. of Theoretical Chemistry, University of *
!          Lund, SWEDEN. Modified to use Schwarz inequality for pre-   *
!          screening, July 1991.                                       *
!          Modified for direct SCF, January '93                        *
!***********************************************************************

use Index_Functions, only: iTri, nTri3_Elem1
use Basis_Info, only: Shells
use iSD_data, only: nSD
use Gateway_Info, only: CutInt, ThrInt
use Symmetry_Info, only: nIrrep
use Int_Options, only: Disc, Disc_Mx, DoFock, DoIntegrals, ExFac, FckNoClmb, FckNoExch, PreSch, Quad_ijkl, Thize, W2Disc
use k2_arrays, only: DeDe, pDq, pFq
use k2_structure, only: IndK2, k2_type, k2data
use Breit, only: nComp
use NDDO, only: twoel_NDDO
use Dens_stuff, only: ipDDij, ipDDik, ipDDil, ipDDjk, ipDDjl, ipDDkl, mDCRij, mDCRik, mDCRil, mDCRjk, mDCRjl, mDCRkl, mDij, mDik, &
                      mDil, mDjk, mDjl, mDkl
use Constants, only: Zero, One, Four, Eight
use Definitions, only: wp, iwp, u6, RtoB, RtoI

implicit none
#include "twoel_interface.fh"
integer(kind=iwp) :: i_Int, iAng(4), iAO(4), iAOst(4), iBasi, iCmp(4), iEta, ijS, ik2, IncEta, IncZet, iOpt, ipAOInt, ipAOInt_, &
                     iPer, iS, iShell(4), iShll(4), ISMAng, iW3, iW4, iW4_, iWR(2), iZeta, jBasj, jk2, jPrInc, jS, kabcd, kBask, &
                     kInts, klS, kS, la, lb, lBasl, lc, ld, lPrInc, lS, mabcd, mabMax, mabMin, mcdMax, mcdMin, mEta, mInts, &
                     mWork2, mZeta, nab, nabcd, nAlpha, nBeta, nByte, ncd, nDCRR, nDCRS, nDelta, nEta, nEta_Tot, nGamma, nInts, &
                     nZeta, nZeta_Tot
real(kind=wp) :: CoorAC(3,2), q4, QInd(2), RST_Triplet, vij, vijkl, vik, vil, vjk, vjl, vkl
logical(kind=iwp) :: ABeqCD, AeqB, AeqC, All_Spherical, Batch_On_Disk, CeqD, Do_TnsCtl, DoAOBatch, DoCoul, DoExch, NoPInts, Pij, &
                     Pijkl, Pik, Pjl, Pkl, Prescreen_On_Int_Only, Scrij, Scrik, Scril, Scrjk, Scrjl, Scrkl, Shijij
real(kind=wp), pointer :: Coeff1(:,:), Coeff2(:,:), Coeff3(:,:), Coeff4(:,:), Dij(:,:), Dik(:,:), Dil(:,:), Djk(:,:), Djl(:,:), &
                          Dkl(:,:)
type(k2_type), pointer :: k2data1(:), k2data2(:)
logical(kind=iwp), parameter :: Copy = .true., NoCopy = .false.
logical(kind=iwp), external :: EQ

#include "macros.fh"
unused_var(SoInt)

nZeta = iSD4(5,1)*iSD4(5,2)
nEta = iSD4(5,3)*iSD4(5,4)

Shijij = ((iSD4(11,1) == iSD4(11,3)) .and. (iSD4(11,2) == iSD4(11,4)))

iShell(:) = iSD4(11,:)
iShll(:) = iSD4(0,:)
iAng(:) = iSD4(1,:)
iAO(:) = iSD4(7,:)
iAOst(:) = iSD4(8,:)

jPrInc = iSD4(6,2)
lPrInc = iSD4(6,4)

!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up pointers to k2 entities.

iS = iShell(1)
jS = iShell(2)
kS = iShell(3)
lS = iShell(4)
ijS = iTri(iS,jS)
klS = iTri(kS,lS)
nDCRR = IndK2(2,ijS)
ik2 = IndK2(3,ijS)
nDCRS = IndK2(2,klS)
jk2 = IndK2(3,klS)
k2data1(1:nDCRR) => k2Data(1:nDCRR,ik2)
k2data2(1:nDCRS) => k2Data(1:nDCRS,jk2)

if (DoFock) then
  Dij(1:mDij,1:mDCRij) => DeDe(ipDDij:ipDDij+mDij*mDCRij-1)
  Dkl(1:mDkl,1:mDCRkl) => DeDe(ipDDkl:ipDDkl+mDkl*mDCRkl-1)
  Dik(1:mDik,1:mDCRik) => DeDe(ipDDik:ipDDik+mDik*mDCRik-1)
  Dil(1:mDil,1:mDCRil) => DeDe(ipDDil:ipDDil+mDil*mDCRil-1)
  Djk(1:mDjk,1:mDCRjk) => DeDe(ipDDjk:ipDDjk+mDjk*mDCRjk-1)
  Djl(1:mDjl,1:mDCRjl) => DeDe(ipDDjl:ipDDjl+mDjl*mDCRjl-1)
else
  ! Dummy association
  Dij(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Dkl(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Dik(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Dil(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Djk(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Djl(1:1,1:1) => DeDe(lbound(DeDe,1):)
end if

All_Spherical = (Shells(iShll(1))%Prjct .and. Shells(iShll(2))%Prjct .and. Shells(iShll(3))%Prjct .and. Shells(iShll(4))%Prjct)

nAlpha = iSD4(5,1)
nBeta = iSD4(5,2)
nGamma = iSD4(5,3)
nDelta = iSD4(5,4)
iBasi = iSD4(19,1)
jBasj = iSD4(19,2)
kBask = iSD4(19,3)
lBasl = iSD4(19,4)

Coeff1(1:nAlpha,1:iBasi) => Shells(iShll(1))%pCff(:,iAOst(1)+1:)
Coeff2(1:nBeta,1:jBasj) => Shells(iShll(2))%pCff(:,iAOst(2)+1:)
Coeff3(1:nGamma,1:kBask) => Shells(iShll(3))%pCff(:,iAOst(3)+1:)
Coeff4(1:nDelta,1:lBasl) => Shells(iShll(4))%pCff(:,iAOst(4)+1:)

#ifdef _DEBUGPRINT_
call RecPrt('Coeff1',' ',Coeff1,nAlpha,iBasi)
call RecPrt('Coeff2',' ',Coeff2,nBeta,jBasj)
call RecPrt('Coeff3',' ',Coeff3,nGamma,kBask)
call RecPrt('Coeff4',' ',Coeff4,nDelta,lBasl)
#endif

RST_triplet = One
QInd(2) = RST_triplet

la = iAng(1)
lb = iAng(2)
lc = iAng(3)
ld = iAng(4)
iSmAng = la+lb+lc+ld

! switch (to generate better start orbitals...)
AeqB = EQ(Coor(:,1),Coor(:,2))
if (twoel_NDDO .and. (.not. AeqB)) return
CeqD = EQ(Coor(:,3),Coor(:,4))
if (twoel_NDDO .and. (.not. CeqD)) return
! switch

AeqC = EQ(Coor(:,1),Coor(:,3))
ABeqCD = AeqB .and. CeqD .and. AeqC
if (ABeqCD .and. (mod(iSmAng,2) == 1)) return
! For Spherical Gaussians, batches like
!(DS|SS), (FP|SS) and (FS|PS) vanish as well
if (ABeqCD .and. All_Spherical .and. (2*max(la,lb,lc,ld) > iSmAng)) return

iCmp(:) = iSD4(2,:)
nab = iCmp(1)*iCmp(2)
ncd = iCmp(3)*iCmp(4)
nabcd = nab*ncd
nInts = nijkl*nabcd
ipAOInt = 1
iW3 = 1+nInts
iW4 = 1

vijkl = k2Data1(1)%abConMax*k2Data2(1)%abConMax

Batch_On_Disk = (vijkl > Thize) .and. (Disc+real(nInts+2+2/RtoI,kind=wp) <= Disc_Mx)

Prescreen_On_Int_Only = PreSch
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
  if (.not. Batch_On_Disk) return
  if (.not. W2Disc) then
    do
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
        return
      else if (QInd(1) < Quad_ijkl) then
        if (mInts > 0) call dRBuf(Wrk(iW3),mInts,NoCopy)
        Disc = Disc+real(2/RtoI+2+mInts,kind=wp)
      else
        call WarningMessage(2,'Twoel: batch is lost!')
        write(u6,*) 'Index,1:',QInd(1),QInd(2),Quad_ijkl,RST_triplet
        call Abend()
      end if
    end do
  end if
end if

! Branch point for partial integral storage

if ((.not. Batch_On_Disk) .or. W2Disc) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Here if the AO batch will be computed!

  ! Compute actual size of the {a0|c0} block

  mabMin = nTri3_Elem1(max(la,lb)-1)
  if (EQ(Coor(:,1),Coor(:,2))) mabMin = nTri3_Elem1(la+lb-1)
  mabMax = nTri3_Elem1(la+lb)-1
  mcdMin = nTri3_Elem1(max(lc,ld)-1)
  if (EQ(Coor(:,3),Coor(:,4))) mcdMin = nTri3_Elem1(lc+ld-1)
  mcdMax = nTri3_Elem1(lc+ld)-1
  mabcd = (mabMax-mabMin+1)*(mcdMax-mcdMin+1)

  ! Find the proper centers to start of with the angular
  ! momentum on. If la == lb there will exist an
  ! ambiguity to which center that angular momentum should
  ! be accumulated on. In that case we will use A and C of
  ! the order as defined by the basis functions types.

  if (la >= lb) then
    CoorAC(:,1) = Coor(:,1)
  else
    CoorAC(:,1) = Coor(:,2)
  end if
  if (lc >= ld) then
    CoorAC(:,2) = Coor(:,3)
  else
    CoorAC(:,2) = Coor(:,4)
  end if

  ! Loops to partition the primitives

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
# ifdef _DEBUGPRINT_
  write(u6,*) 'nZeta=',nZeta
  write(u6,*) 'nEta=',nEta
  write(u6,*) 'nZeta_Tot, IncZet=',nZeta_Tot,IncZet
  write(u6,*) 'nEta_Tot,  IncEta=',nEta_Tot,IncEta
# endif

  kabcd = 0
  Do_TnsCtl = .false.
  NoPInts = .true.
  ipAOInt_ = ipAOInt
  iW4_ = iW4

  do iZeta=1,nZeta_Tot,IncZet
    mZeta = min(IncZet,nZeta_Tot-iZeta+1)
    if (all(Coeff2(:,:) == Zero)) cycle

    do iEta=1,nEta_Tot,IncEta
      mEta = min(IncEta,nEta_Tot-iEta+1)
      if (all(Coeff4(:,:) == Zero)) cycle

      call DrvRys(iZeta,iEta,nZeta,nEta,mZeta,mEta,nZeta_Tot,nEta_Tot,k2data1(1),k2data2(1),nAlpha,nBeta,nGamma,nDelta,1,1,1,1,1, &
                  1,ThrInt,CutInt,vij,vkl,vik,vil,vjk,vjl,Prescreen_On_Int_Only,NoInts,iAng,Coor,CoorAC,mabMin,mabMax,mcdMin, &
                  mcdMax,nijkl/nComp,nabcd,mabcd,Wrk,ipAOInt_,iW4_,nWork2,mWork2,k2Data1(1)%HrrMtrx(:,1),k2Data2(1)%HrrMtrx(:,1), &
                  la,lb,lc,ld,iCmp,iShll,NoPInts,Dij(:,1),mDij,Dkl(:,1),mDkl,Do_TnsCtl,kabcd,Coeff1,iBasi,Coeff2,jBasj,Coeff3, &
                  kBask,Coeff4,lBasl)

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
    return
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

    Wrk(iW3:iW3+nijkl*nabcd-1) = Wrk(ipAOInt:ipAOInt+nijkl*nabcd-1)
    call DGeTMO(Wrk(iW3),nabcd,nabcd,nijkl,Wrk(ipAOInt),nijkl)

  end if
# ifdef _DEBUGPRINT_
  call RecPrt('(AB|CD)',' ',Wrk(ipAOInt),nijkl/nComp,nComp*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4))
# endif

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
end if

if (Batch_On_Disk .and. (.not. W2Disc)) then
  do
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
      if (mInts == 0) return
      exit
    else if (QInd(1) < Quad_ijkl) then
      if (mInts > 0) call dRBuf(Wrk(iW3),mInts,NoCopy)
      Disc = Disc+real(2/RtoI+2+mInts,kind=wp)
    else
      call WarningMessage(2,'Twoel: batch is lost!')
      write(u6,*) 'Index,1:',QInd(1),QInd(2),Quad_ijkl,RST_triplet
      call Abend()
    end if
  end do

  iOpt = 0
  call UpkR8(iOpt,nInts,nByte,Wrk(iW3),Wrk(ipAOInt))
end if

! Accumulate contributions directly to the Fock matrix.

if (DoFock) call FckAcc_NoSymq(iCmp(1),iCmp(2),iCmp(3),iCmp(4),Shijij,iShell,nijkl,Wrk(ipAOInt),pFq,pDq,size(pFq),iAO,iAOst, &
                               iBasi,jBasj,kBask,lBasl,DoCoul,DoExch,vij,vkl,vik,vil,vjk,vjl,ExFac)

if (DoIntegrals) then
  if (ipAOInt /= 1) then
    Wrk(1:nijkl*nabcd) = Wrk(ipAOInt:ipAOInt+nijkl*nabcd-1)
    ipAOInt = 1
  end if
  iPer = 1
  Pij = (iSD4(20,1) == iSD4(20,2))
  Pkl = (iSD4(20,3) == iSD4(20,4))
  Pik = (iSD4(20,1) == iSD4(20,3))
  Pjl = (iSD4(20,2) == iSD4(20,4))
  Pijkl = (Pij .and. Pkl .and. Pik .and. Pjl)
  if (Pij) iPer = iPer*2
  if (Pkl) iPer = iPer*2
  if (Pijkl) iPer = iPer*2
  q4 = Eight/real(iPer,kind=wp)
  if (nIrrep == 1) q4 = One
  if (q4 /= One) Wrk(ipAOInt:ipAOInt+nijkl*nabcd-1) = q4*Wrk(ipAOInt:ipAOInt+nijkl*nabcd-1)
end if

end subroutine TwoEl_NoSym
