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
! Copyright (C) 1990, Roland Lindh                                     *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine TwoEl_mck(Coor,nRys,Hess,nHess,IfGrd,IndGrd,IfHss,IndHss,IfG,PSO,nijkl,nPSO,Work2,nWork2,Work3,nWork3,Work4,nWork4,Aux, &
                     nAux,WorkX,nWorkX,Fin,nfin,Temp,nTemp,nTwo2,nFt,TwoHam,Buffer,nBuffer,lgrad,ldot,n8,ltri,Dan,Din,moip,naco, &
                     rMOIN,nMOIN,iSD4)
!***********************************************************************
!                                                                      *
!     Input:                                                           *
!     PSO                                                              *
!     Work2                                                            *
!     Work3                                                            *
!     Work4                                                            *
!     AUX                                                              *
!     Fin    : Area for cntrctd int in sph hmn                         *
!     Temp   : Working place for F gen and n8                          *
!     TwoHam : Final results fock matrix and MO's                      *
!                                                                      *
!     Object:      To construct the first order derivatives of the AO- *
!     integrals and add them up to the MO derivatives and              *
!     the Fock matrix derivatives and contract the second              *
!     order derivatives of the AO's with the second order              *
!     density matrix.                                                  *
!                                                                      *
!     Authors: Roland Lindh, IBM Almaden Research Center, San Jose, CA *
!     March '90                                                        *
!     Anders Bernhardsson Theoretical Chemistry 95                     *
!***********************************************************************
!                                                                      *
!     When we are calculating the second order derivatives we need     *
!     the derivatives of the two electron integrals in three ways:     *
!                                                                      *
!     (2)                                                              *
!     1)  To calculate the static term H                               *
!     -                                                                *
!     2)  To calculate the non-zero part of <0|[E  ,H]|0>              *
!     pq                                                               *
!                                                                      *
!                                                                      *
!     3)  To calculate the derivatives of the MO orbitals with all     *
!     four indexes in the active space.                                *
!                                                                      *
!     In this implementation all contributions are calculated at       *
!     the same time.                                                   *
!                                                                      *
!     (2)                                                              *
!     The H    is calculated by contracting the second order           *
!     derivatives on the flight with the second order density matrix   *
!                                                                      *
!     (1)  (1)       (1)                                               *
!     the    F  - F    and MO     are calculated by first contracting  *
!     pq   qp                                                          *
!                                                                      *
!     the primitives and transform the integrals to  spherical         *
!     harmonics and then construct the Fock matrix as a direct SCF     *
!     The Fock matrixes is transformed to MO base and then added up to *
!     total Fock matrix                                                *
!                                                                      *
!***********************************************************************

use McKinley_global, only: CPUStat, nIntegrals, nScreen, nTrans, nTwoDens, PreScr
use iSD_data, only: nSD
use Index_Functions, only: iTri, nTri_Elem1
use k2_structure, only: Indk2, k2_type, k2Data
use k2_arrays, only: BraKet, DeDe, DeDe2
use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: Shells
use Phase_Info, only: iPhase
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use Rys_interfaces, only: cff2d_kernel, modu2_kernel, tval1_kernel
use Dens_stuff, only: ipDDij, ipDDij2, ipDDik, ipDDik2, ipDDil, ipDDil2, ipDDjk, ipDDjk2, ipDDjl, ipDDjl2, ipDDkl, ipDDkl2, &
                      mDCRij, mDCRik, mDCRil, mDCRjk, mDCRjl, mDCRkl, mDij, mDik, mDil, mDjk, mDjl, mDkl
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nRys, nHess, IndGrd(3,4,0:7), IndHss(4,3,4,3,0:7), nijkl, nPSO, nWork2, nWork3, nWork4, nAux, &
                                 nWorkX, nfin, nTemp, nTwo2, nFt, nBuffer, moip(0:7), naco, nMOIN, iSD4(0:nSD,4)
real(kind=wp), intent(in) :: Coor(3,4), PSO(nijkl,nPSO), Dan(*), Din(*)
real(kind=wp), intent(inout) :: Hess(nHess), WorkX(nWorkX), TwoHam(nTwo2), Buffer(nBuffer), rMOIN(nMOIN)
logical(kind=iwp), intent(in) :: IfGrd(3,4), IfHss(4,3,4,3), lgrad, ldot, n8, ltri
logical(kind=iwp), intent(out) :: IfG(4)
real(kind=wp), intent(out) :: Work2(nWork2), Work3(nWork3), Work4(nWork4), Aux(nAux), Fin(nfin), Temp(nTemp)
integer(kind=iwp) :: iAngV(4), iAO(4), iAOst(4), iBasi, iCar, iCmp(4), iCmpa, iCNT, iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), iDCRTS, &
                     iEta, iIrr, ijS, ik2, IncEta, IncZet, Indx(3,4), ip, ip2, ipFT, ipS1, ipS2, ipTemp, iS, iShell(4), iShll(4), &
                     iShlla, iuvwx(4), ix2, iy2, iz2, iZeta, jBasj, jCmpb, jk2, JndGrd(3,4,0:7), JndHss(4,3,4,3,0:7), jPrInc, jS, &
                     jShllb, kBask, kCmpc, klS, kS, kShllc, la, lb, lBasl, lc, lCmpd, ld, lDCR1, lDCR2, lDCRR, lDCRS, lDCRT, lEta, &
                     lPrInc, lS, lShlld, lZeta, mab, mcd, mEta, mZeta, n, nabcd, nAlpha, nBeta, nDCR1, nDCR2, nDCRR, nDCRS, nDCRT, &
                     nDelta, nEta, nEta_Tot, nGamma, nGr, niag, nOp(4), nS1, nS2, nTe, nw3, nw3_2, nZeta, nZeta_Tot
real(kind=wp) :: CoorAC(3,2), CoorM(3,4), dum1, dum2, dum3, Fact, Time
logical(kind=iwp) :: ABeqCD, AeqB, AeqC, CeqD, first, JfGrd(3,4), JfHss(4,3,4,3), l_og, ldot2, Tr(4), Shijij
procedure(cff2d_kernel) :: Cff2D
procedure(modu2_kernel) :: ModU2
procedure(tval1_kernel) :: TERI1
type(k2_type), pointer :: k2data1(:), k2data2(:)
real(kind=wp), pointer :: Coeff1(:,:), Coeff2(:,:), Coeff3(:,:), Coeff4(:,:), Dij1(:,:), Dij2(:,:), Dik1(:,:), Dik2(:,:), &
                          Dil1(:,:), Dil2(:,:), Djk1(:,:), Djk2(:,:), Djl1(:,:), Djl2(:,:), Dkl1(:,:), Dkl2(:,:)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ

!                                                                      *
!***********************************************************************
!                                                                      *
! PROLOGUE
!                                                                      *
!***********************************************************************
!                                                                      *
jPrInc = iSD4(6,2)
lPrInc = iSD4(6,4)

iAO(:) = iSD4(7,:)
iCmp(:) = iSD4(2,:)
iShll(:) = iSD4(0,:)
iShell(:) = iSD4(11,:)
iAngV(:) = iSD4(1,:)
iAOst(:) = iSD4(8,:)

nAlpha = iSD4(5,1)
nBeta = iSD4(5,2)
nGamma = iSD4(5,3)
nDelta = iSD4(5,4)

iBasi = iSD4(19,1)
jBasj = iSD4(19,2)
kBask = iSD4(19,3)
lBasl = iSD4(19,4)

nZeta = iSD4(5,1)*iSD4(5,2)
nEta = iSD4(5,3)*iSD4(5,4)

Shijij = ((iSD4(11,1) == iSD4(11,3)) .and. (iSD4(11,2) == iSD4(11,4)))

nGr = 0
ldot2 = ldot
la = iAngV(1)
lb = iAngV(2)
lc = iAngV(3)
ld = iAngV(4)
iCmpa = iCmp(1)
jCmpb = iCmp(2)
kCmpc = iCmp(3)
lCmpd = iCmp(4)
iShlla = iShll(1)
jShllb = iShll(2)
kShllc = iShll(3)
lShlld = iShll(4)
IncZet = nAlpha*jPrInc
IncEta = nGamma*lPrInc
nabcd = iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)
mab = nTri_Elem1(la)*nTri_Elem1(lb)
mcd = nTri_Elem1(lc)*nTri_Elem1(ld)

! Scratch space for Fock Matrix construction

ip = 1
ipS1 = ip
nS1 = nijkl*nabcd
ip = ip+nS1
ipS2 = ip
nS2 = max(nS1,nijkl+max(iBasi*lBasl,jBasj*lBasl,iBasi*kBask,jBasj*kBask))
ip = ip+nS2
ipFT = ip
ip = ip+nFT
ipTemp = ip
nTe = nijkl*nabcd
ip = ip+nTe
if (ip-1 > nTemp) then
  write(u6,*) 'TwoEl_McK: ip-1 > nTemp'
  write(u6,*) 'ip,nTemp=',ip,nTemp
  call Abend()
end if
#ifdef _WARNING_WORKAROUND_
! Avoid some warnings about unset output arguments
Temp(nTemp) = One
#endif
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
nDCR1 = IndK2(2,ijS)
ik2 = IndK2(3,ijS)
nDCR2 = IndK2(2,klS)
jk2 = IndK2(3,klS)
k2data1(1:nDCR1) => k2Data(1:nDCR1,ik2)
k2data2(1:nDCR2) => k2Data(1:nDCR2,jk2)

Coeff1(1:nAlpha,1:iBasi) => Shells(iShll(1))%pCff(:,iAOst(1)+1:)
Coeff2(1:nBeta,1:jBasj) => Shells(iShll(2))%pCff(:,iAOst(2)+1:)
Coeff3(1:nGamma,1:kBask) => Shells(iShll(3))%pCff(:,iAOst(3)+1:)
Coeff4(1:nDelta,1:lBasl) => Shells(iShll(4))%pCff(:,iAOst(4)+1:)

if (lbound(DeDe,1) >= 1) then
  Dij1(1:mDij,1:mDCRij) => DeDe(ipDDij:ipDDij+mDij*mDCRij-1)
  Dkl1(1:mDkl,1:mDCRkl) => DeDe(ipDDkl:ipDDkl+mDkl*mDCRkl-1)
  Dik1(1:mDik,1:mDCRik) => DeDe(ipDDik:ipDDik+mDik*mDCRik-1)
  Dil1(1:mDil,1:mDCRil) => DeDe(ipDDil:ipDDil+mDil*mDCRil-1)
  Djk1(1:mDjk,1:mDCRjk) => DeDe(ipDDjk:ipDDjk+mDjk*mDCRjk-1)
  Djl1(1:mDjl,1:mDCRjl) => DeDe(ipDDjl:ipDDjl+mDjl*mDCRjl-1)
else
  ! dummy associations
  Dij1(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Dkl1(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Dik1(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Dil1(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Djk1(1:1,1:1) => DeDe(lbound(DeDe,1):)
  Djl1(1:1,1:1) => DeDe(lbound(DeDe,1):)
end if

if (lbound(DeDe2,1) >= 1) then
  Dij2(1:mDij,1:mDCRij) => DeDe2(ipDDij2:ipDDij2+mDij*mDCRij-1)
  Dkl2(1:mDkl,1:mDCRkl) => DeDe2(ipDDkl2:ipDDkl2+mDkl*mDCRkl-1)
  Dik2(1:mDik,1:mDCRik) => DeDe2(ipDDik2:ipDDik2+mDik*mDCRik-1)
  Dil2(1:mDil,1:mDCRil) => DeDe2(ipDDil2:ipDDil2+mDil*mDCRil-1)
  Djk2(1:mDjk,1:mDCRjk) => DeDe2(ipDDjk2:ipDDjk2+mDjk*mDCRjk-1)
  Djl2(1:mDjl,1:mDCRjl) => DeDe2(ipDDjl2:ipDDjl2+mDjl*mDCRjl-1)
else
  ! dummy associations
  Dij2(1:1,1:1) => DeDe2(lbound(DeDe2,1):)
  Dkl2(1:1,1:1) => DeDe2(lbound(DeDe2,1):)
  Dik2(1:1,1:1) => DeDe2(lbound(DeDe2,1):)
  Dil2(1:1,1:1) => DeDe2(lbound(DeDe2,1):)
  Djk2(1:1,1:1) => DeDe2(lbound(DeDe2,1):)
  Djl2(1:1,1:1) => DeDe2(lbound(DeDe2,1):)
end if

!                                                                      *
!***********************************************************************
!                                                                      *
!           - - - - - - END PROLOGUE - - - - - -
!                                                                      *
!***********************************************************************
!                                                                      *
call mk_DCRs_and_Stabilizers(Fact,iuvwx,nDCRR,nDCRS,nDCRT,iDCRR,iDCRS,iDCRT,nSD,iSD4)
!                                                                      *
!***********************************************************************
!                                                                      *
! - - - - Loop over first set
!                                                                      *
!***********************************************************************
!                                                                      *
nOp(1) = NrOpr(0)
CoorM(:,1) = Coor(:,1)
do lDCRR=0,nDCRR-1
  nOp(2) = NrOpr(iDCRR(lDCRR))
  call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
  AeqB = EQ(CoorM(1,1),CoorM(1,2))
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! - - - - Loop over second set
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do lDCRS=0,nDCRS-1
    Coorm(:,3) = Coor(:,3)
    call OA(iDCRS(lDCRS),Coor(1:3,4),CoorM(1:3,4))
    CeqD = EQ(Coor(1,3),CoorM(1,4))
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! - - - - Loop over third set
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    do lDCRT=nDCRT-1,0,-1

      nOp(3) = NrOpr(iDCRT(lDCRT))
      nOp(4) = NrOpr(ieor(iDCRT(lDCRT),iDCRS(lDCRS)))

      iDCRTS = ieor(iDCRT(lDCRT),iDCRS(lDCRS))
      call OA(iDCRTS,Coor(1:3,4),CoorM(1:3,4))
      call OA(iDCRT(lDCRT),Coor(1:3,3),CoorM(1:3,3))

      AeqC = EQ(CoorM(1,1),CoorM(1,3))
      ABeqCD = AeqB .and. CeqD .and. AeqC
      ! No contribution to geometric derivatives from one-center integrals
      if (ABeqCD) cycle

      ! Find the proper centers to start of with the angular
      ! momentum on. If la == lb there will exist an
      ! ambiguity to which center that angular momentum should
      ! be accumulated on. In that case we will use A and C of
      ! the order as defined by the basis functions types.

      if (iAngV(1) >= iAngV(2)) then
        CoorAC(:,1) = CoorM(:,1)
      else
        CoorAC(:,1) = CoorM(:,2)
      end if
      if (iAngV(3) >= iAngV(4)) then
        CoorAC(:,2) = CoorM(:,3)
      else
        CoorAC(:,2) = CoorM(:,4)
      end if

      ! Calculate the desymmetrized two-electron density matrix in cartesian AO base.

      call Timing(dum1,Time,dum2,dum3)
      if (ldot2) call TwoDns(iAngV,iCmp,shijij,ishll,ishell,iAO,nOp,iBasi,jBasj,kBask,lBasl,Aux,nAux,Work2,nWork2,Work3,nWork3, &
                             Work4,nWork4,PSO,nPSO,Fact)

      call Timing(dum1,Time,dum2,dum3)
      CpuStat(nTwoDens) = CpuStat(nTwoDens)+Time
      !----------------------------------------------------------------*
      !
      ! Fix the control matrices for derivatives and try to use
      ! translation invariance as efficiently as possible.
      !
      !----------------------------------------------------------------*
      JfHss(:,:,:,:) = IfHss
      JfGrd(:,:) = IfGrd
      ifg(:) = .true.
      Tr(:) = .false.
      JndHss(:,:,:,:,0:nIrrep-1) = IndHss(:,:,:,:,0:nIrrep-1)
      JndGrd(:,:,0:nIrrep-1) = IndGrd(:,:,0:nIrrep-1)

      ! Delete one center that should be calculated with translation invariance

      call Translation(ifg,jfgrd,jfhss,tr,jndgrd,jndhss,coorm,nirrep,indgrd,indhss)

      if (.not. ldot) then
        JfHss(:,:,:,:) = .false.
        JndHss(:,:,:,:,:) = 0
      end if

      !----------------------------------------------------------------*
      !
      ! Loops to partition the primitives
      !
      !----------------------------------------------------------------*
      lDCR1 = NrOpr(iDCRR(lDCRR))+1
      lDCR2 = NrOpr(iDCRS(lDCRS))+1
      ix2 = iPhase(1,iDCRT(lDCRT))
      iy2 = iPhase(2,iDCRT(lDCRT))
      iz2 = iPhase(3,iDCRT(lDCRT))

      nZeta_Tot = k2Data1(lDCR1)%IndZ(nZeta+1)
      nEta_Tot = k2Data2(lDCR2)%IndZ(nEta+1)

      first = .true.
      nGr = 0
      do iZeta=1,nZeta_Tot,IncZet
        mZeta = min(IncZet,nZeta_Tot-iZeta+1)
        ! Check that subblock of contraction matrix has non-zero elements.
        if (all(Coeff2(:,:) == Zero)) cycle
        do iEta=1,nEta_Tot,IncEta
          mEta = min(IncEta,nEta_Tot-iEta+1)
          ! Check that subblock of contraction matrix has non-zero elements.
          if (all(Coeff4(:,:) == Zero)) cycle
          !------------------------------------------------------------*
          !     PRE PRESCREENING                                       *
          !------------------------------------------------------------*

          lZeta = mZeta
          lEta = mEta

          ! Decontract the 2nd order density matrix

          ! Work4->Work2  Work3:scratch
          call Timing(dum1,Time,dum2,dum3)
          if (ldot2) call Tcrtnc_h(Coeff1,nAlpha,iBasi,Coeff2,nBeta,jBasj,Coeff3,nGamma,kBask,Coeff4,nDelta,lBasl, &
                                   Work4,mab*mcd,Work3,nWork3/2,Work2, &
                                   k2Data1(lDCR1)%IndZ(iZeta:iZeta+mZeta-1),mZeta, &
                                   k2Data2(lDCR2)%IndZ(iEta:iEta+mEta-1),mEta)
          call Timing(dum1,Time,dum2,dum3)
          CPUStat(nTwoDens) = CPUStat(nTwoDens)+Time

          ! Transfer k2 data and prescreen

          ! Work2:PAO-> Work2
          ! Work3 Scratch
          call Timing(dum1,Time,dum2,dum3)
          call Screen_mck(iZeta-1,iEta-1,Work2,Work3,mab*mcd,nZeta,nEta,mZeta,mEta,lZeta,lEta, &
                          k2Data1(lDCR1),k2Data2(lDCR2), &
                          BraKet%Zeta,BraKet%ZInv,BraKet%P,BraKet%xA,BraKet%xB,BraKet%KappaAB, &
                          BraKet%Eta,BraKet%EInv,BraKet%Q,BraKet%xG,BraKet%xD,BraKet%KappaCD, &
                          BraKet%xpre,1,1,1,ix2,iy2,iz2,CutInt,PreScr,Braket%IndZet,Braket%IndEta,ldot2)
          call Timing(dum1,Time,dum2,dum3)
          CPUStat(nScreen) = CPUStat(nScreen)+Time

          if (lZeta*lEta == 0) cycle

          ! Compute integral derivative and accumulate
          ! contribution to the molecular gradient.

          ! Work2:PAO
          ! Work3:Work area  The PO integrals are stored in the begining of Work3

          call Timing(dum1,Time,dum2,dum3)

          call Rysg2(iAngV,nRys,lZeta*lEta,BraKet%xA,BraKet%xB,BraKet%xG,BraKet%xD,BraKet%Zeta,BraKet%ZInv,lZeta,BraKet%Eta, &
                     BraKet%EInv,lEta,BraKet%P,nZeta,BraKet%Q,nEta,CoorM,CoorM,CoorAC,Work3,nWork3,TERI1,ModU2,Cff2D,Work2, &
                     mab*mcd,Hess,nHess,JfGrd,JndGrd,JfHss,JndHss,nOp,iuvwx,IfG,nGr,Indx,lgrad,ldot,Tr)
          call Timing(dum1,Time,dum2,dum3)
          CPUStat(nIntegrals) = CPUStat(nIntegrals)+Time

          ! Work3 AO
          ! Work3_3  Scratch
          ! ->    Work3_2
          !------------------------------------------------------------*
          !
          ! Transform integrals to AO base
          !
          !------------------------------------------------------------*
          ip2 = nGr*mab*mcd*lZeta*lEta+1
          call Timing(dum1,Time,dum2,dum3)
          call Cntrct_mck(First,Coeff1,nAlpha,iBasi,Coeff2,nBeta,jBasj,Coeff3,nGamma,kBask,Coeff4,nDelta,lBasl,Work3,nGr*mab*mcd, &
                          Work3(ip2),nwork3-ip2,BraKet%xpre,WorkX,nWorkX,lZeta*lEta,BraKet%IndZet,nZeta,lZeta,BraKet%IndEta,nEta, &
                          lEta)
        end do
      end do

      ! Mark which derivatives should be calculated with translation invariance.

      if (nGr == 0) cycle
      do iCNT=1,4
        if (Tr(iCnt)) then
          do iCar=1,3
            l_og = .false.
            do iIrr=0,nIrrep-1
              l_og = l_og .or. (indgrd(iCar,iCnt,iIrr) /= 0)
            end do
            if (l_og) Indx(iCar,iCnt) = -1
          end do
        end if
      end do

      if (Fact /= One) then
        n = nGr*mab*mcd*nijkl
        WorkX(1:n) = Fact*WorkX(1:n)
      end if

      !----------------------------------------------------------------*
      !
      !     Transpose abcd,g,IJKL -> bcd,g,IJKL,A Work3 -> Work3_2
      !
      !----------------------------------------------------------------*

      niag = nijkl*nTri_Elem1(lb)*mcd*nGr
      call CrSph_mck(WorkX,niag,nTri_Elem1(la),RSph(ipSph(la)),la,Shells(iShlla)%Transf,Shells(iShlla)%Prjct,Work3,iCmpa)
      nw3 = niag*iCmpa
      ip2 = 1+nw3

      !----------------------------------------------------------------*
      !
      !     Transpose   bcd,g,IJKL,A -> cd,g,IJKL,AB Work3_2->Work3
      !
      !----------------------------------------------------------------*

      niag = nijkl*mcd*nGr*iCmpa
      nw3_2 = niag*jCmpb
      if (nw3+nw3_2 > nWork3) then
        write(u6,*) '1: nw3+nw3_2 > nWork3'
        call Abend()
      end if
      call CrSph_mck(Work3,niag,nTri_Elem1(lb),RSph(ipSph(lb)),lb,Shells(jShllb)%Transf,Shells(jShllb)%Prjct,Work3(ip2),jCmpb)

      !----------------------------------------------------------------*
      !
      !     Transpose  cd,g,IJKL,AB -> d,g,IJKL,ABC  Work3->Work3_2
      !
      !----------------------------------------------------------------*

      niag = nijkl*nGr*nTri_Elem1(ld)*iCmpa*jCmpb
      call CrSph_mck(Work3(ip2),niag,nTri_Elem1(lc),RSph(ipSph(lc)),lc,Shells(kShllc)%Transf,Shells(kShllc)%Prjct,Work3,kCmpc)
      if (niag*kCmpc > nw3) then
        write(u6,*) 'niag*kCmpc > nw3'
        call Abend()
      end if
      nw3 = niag*kCmpc
      ip2 = nw3+1

      !----------------------------------------------------------------*
      !
      !     Transpose   d,g,IJKL,ABC -> g,IJKL,ABCD Work3_2->Work3
      !
      !----------------------------------------------------------------*

      niag = nijkl*nGr*iCmpa*jCmpb*kCmpc
      nw3_2 = niag*lCmpd
      if (nw3+nw3_2 > nWork3) then
        write(u6,*) '2: nw3+nw3_2 > nWork3'
        call Abend()
      end if
      call CrSph_mck(Work3,niag,nTri_Elem1(ld),RSph(ipSph(ld)),ld,Shells(lShlld)%Transf,Shells(lShlld)%Prjct,Work3(ip2),lCmpd)

      !----------------------------------------------------------------*
      !
      !     Transpose g,IJKL,ABCD -> IJKL,ABCD,g Work3->Buffer
      !
      !----------------------------------------------------------------*

      niag = nijkl*iCmpa*jCmpb*kCmpc*lCmpd
      call DGetMO(Work3(ip2),nGr,nGr,niag,Fin,niag)

      ! DEBUG  (calculates gradient from transformed integrals)

      call Timing(dum1,Time,dum2,dum3)
      CPUStat(nTrans) = CPUStat(nTrans)+Time

      !----------------------------------------------------------------*
      !
      !     Send the integrals to clrbuffer for construction of
      !
      !----------------------------------------------------------------*

      call ClrBuf(idcrr(ldcrr),idcrs(ldcrs),idcrt(ldcrt),nGr,Shijij,iAngV,iCmp,iShll,iShell,iShell,iBasi,jBasj,kBask,lBasl,Dij1, &
                  Dij2,mDij,mDCRij,Dkl1,Dkl2,mDkl,mDCRkl,Dik1,Dik2,mDik,mDCRik,Dil1,Dil2,mDil,mDCRil,Djk1,Djk2,mDjk,mDCRjk,Djl1, &
                  Djl2,mDjl,mDCRjl,fin,nfin,Temp(ipFT),nFT,Temp(ipS1),nS1,Temp(ipS2),nS2,Temp(ipTemp),nTe,TwoHam,nTwo2,JndGrd, &
                  Indx,iao,iaost,iuvwx,n8,ltri,moip,nAcO,rMoin,nmoin,ntemp,Buffer,nOp,Din,Dan)

    end do

  end do

end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Twoel_Mck
