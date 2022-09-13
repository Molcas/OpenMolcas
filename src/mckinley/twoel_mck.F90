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

subroutine TwoEl_mck(Coor,iAngV,iCmp,iShell,iShll,iAO,iAOst,iStb,jStb,kStb,lStb,nRys,Data1,nData1,Data2,nData2,Pren,Prem,nAlpha, &
                     nBeta,jPrInc,nGamma,nDelta,lPrInc,Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl,Zeta,ZInv,P,rKab,nZeta, &
                     Eta,EInv,Q,rKcd,nEta,xA,xB,xG,xD,xPre,Hess,nHess,IfGrd,IndGrd,IfHss,IndHss,IfG,PSO,nPSO,Work2,nWork2,Work3, &
                     nWork3,Work4,nWork4,Aux,nAux,WorkX,nWorkX,Shijij,Dij1,Dij2,mDij,nDij,Dkl1,Dkl2,mDkl,nDkl,Dik1,Dik2,mDik,nDik, &
                     Dil1,Dil2,mDil,nDil,Djk1,Djk2,mDjk,nDjk,Djl1,Djl2,mDjl,nDjl,icmpi,Fin,nfin,Temp,nTemp,nTwo2,nFt,IndZet, &
                     IndEta,TwoHam,Buffer,nBuffer,lgrad,ldot,n8,ltri,Dan,Din,moip,naco,rMOIN,nMOIN,new_fock)
!***********************************************************************
!                                                                      *
!     Input:                                                           *
!     Data1                                                            *
!     Data2                                                            *
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

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use McKinley_global, only: CPUStat, nIntegrals, nScreen, nTrans, nTwoDens, PreScr
use Index_Functions, only: nTri_Elem1
use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: MolWgh, Shells
use Center_Info, only: dc
use Phase_Info, only: iPhase
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
#include "ndarray.fh"
integer(kind=iwp), intent(in) :: iAngV(4), iCmp(4), iShell(4), iShll(4), iAO(4), iAOst(4), iStb, jStb, kStb, lStb, nRys, nData1, &
                                 nData2, nAlpha, nBeta, jPrInc, nGamma, nDelta, lPrInc, iBasi, jBasj, kBask, lBasl, nZeta, nEta, &
                                 nHess, IndGrd(3,4,0:7), IndHss(4,3,4,3,0:7), nPSO, nWork2, nWork3, nWork4, nAux, nWorkX, mDij, &
                                 nDij, mDkl, nDkl, mDik, nDik, mDil, nDil, mDjk, nDjk, mDjl, nDjl, icmpi(4), nfin, nTemp, nTwo2, &
                                 nFt, nBuffer, moip(0:7), naco, nMOIN
real(kind=wp), intent(in) :: Coor(3,4), Data1(nZeta*nDArray+nDScalar,nData1), Data2(nEta*nDArray+nDScalar,nData2), &
                             Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj), Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl), &
                             PSO(iBasi*jBasj*kBask*lBasl,nPSO), Dij1(mDij,nDij), Dij2(mDij,nDij), Dkl1(mDkl,nDkl), &
                             Dkl2(mDkl,nDkl), Dik1(mDik,nDik), Dik2(mDik,nDik), Dil1(mDil,nDil), Dil2(mDil,nDil), Djk1(mDjk,nDjk), &
                             Djk2(mDjk,nDjk), Djl1(mDjl,nDjl), Djl2(mDjl,nDjl), Dan(*), Din(*)
real(kind=wp), intent(inout) :: Pren, Prem, Hess(nHess), WorkX(nWorkX), TwoHam(nTwo2), Buffer(nBuffer), rMOIN(nMOIN)
real(kind=wp), intent(out) :: Zeta(nZeta), ZInv(nZeta), P(nZeta,3), rKab(nZeta), Eta(nEta), EInv(nEta), Q(nEta,3), rKcd(nEta), &
                              xA(nZeta), xB(nZeta), xG(nEta), xD(nEta), xpre(nGamma*nDelta*nAlpha*nBeta), Work2(nWork2), &
                              Work3(nWork3), Work4(nWork4), Aux(nAux), Fin(nfin), Temp(nTemp)
logical(kind=iwp), intent(in) :: IfGrd(3,4), IfHss(4,3,4,3), Shijij, lgrad, ldot, n8, ltri, new_fock
logical(kind=iwp), intent(out) :: IfG(4)
integer(kind=iwp), intent(out) :: IndZet(nAlpha*nBeta), IndEta(nGamma*nDelta)
integer(kind=iwp) :: iCmpa, iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), iDCRTS, IncEta, IncZet, Indx(3,4), ip, ip2, ipFT, ipS1, ipS2, &
                     ipTemp, iShlla, iStabM(0:7), iStabN(0:7), iuvwx(4), ix2, iy2, iz2, jCmpb, JndGrd(3,4,0:7), &
                     JndHss(4,3,4,3,0:7), jShllb, kCmpc, kShllc, la, lb, lc, lCmpd, ld, lDCR1, lDCR2, lEta, LmbdR, LmbdS, LmbdT, &
                     lShlld, lStabM, lStabN, lZeta, mab, mcd, mEta, mZeta, n, nabcd, nDCRR, nDCRS, nDCRT, nEta_Tot, nGr, niag, &
                     nijkl, nOp(4), nS1, nS2, nTe, nw3, nw3_2, nZeta_Tot
real(kind=wp) :: CoorAC(3,2), CoorM(3,4), dum1, dum2, dum3, Fact, FactNd, Time, u, v, w, x
logical(kind=iwp) :: ABeq, ABeqCD, AeqB, AeqC, CDeq, CeqD, first, JfGrd(3,4), JfHss(4,3,4,3), l_og, ldot2, no_integrals, Tr(4)
integer(kind=iwp), external :: ip_abMax, ip_IndZ, ip_Z, NrOpr
logical(kind=iwp), external :: EQ, lEmpty
external :: TERI1, ModU2, Cff2D

!                                                                      *
!***********************************************************************
!                                                                      *
call TwoEl_mck_Internal(Data1,Data2)

! This is to allow type punning without an explicit interface
contains

subroutine TwoEl_mck_Internal(Data1,Data2)

  real(kind=wp), target :: Data1(nZeta*nDArray+nDScalar,nData1), Data2(nEta*nDArray+nDScalar,nData2)
  integer(kind=iwp), pointer :: iData1(:), iData2(:)
  integer(kind=iwp) :: iCar, iCNT, iEta, iIrr, iZeta, lDCRR, lDCRS, lDCRT
  !Bug in gcc 7: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94270
# ifdef _WARNING_WORKAROUND_
  interface
    subroutine Rysg2(iAnga,nRys,nT,Alpha,Beta,Gmma,Delta,Zeta,ZInv,nZeta,Eta,EInv,nEta,P,lP,Q,lQ,Coori,Coora,CoorAC,Array,nArray, &
                     Tvalue,ModU2,Cff2D,PAO,nPAO,Hess,nHess,IfGrd,IndGrd,IfHss,IndHss,nOp,iuvwx,IfG,mVec,Index_Out,lGrad,lHess,Tr)
      use Definitions, only: wp, iwp
      integer(kind=iwp), intent(in) :: iAnga(4), nRys, nT, nZeta, nEta, lP, lQ, nArray, nPAO, nHess, IndGrd(3,4,0:7), nOp(4), &
                                       iuvwx(4)
      real(kind=wp), intent(in) :: Alpha(nZeta), Beta(nZeta), Gmma(nEta), Delta(nEta), Zeta(nZeta), ZInv(nZeta), Eta(nEta), &
                                   EInv(nEta), P(lP,3), Q(lQ,3), Coori(3,4), Coora(3,4), CoorAC(3,2), PAO(nT,nPAO)
      real(kind=wp), intent(inout) :: Array(nArray), Hess(nHess)
      external :: Tvalue, ModU2, Cff2D
      logical(kind=iwp), intent(inout) :: IfGrd(3,4), IfHss(4,3,4,3), IfG(4), Tr(4)
      integer(kind=iwp), intent(inout) :: IndHss(4,3,4,3,0:7)
      integer(kind=iwp), intent(out) :: mVec, Index_Out(3,4)
      logical(kind=iwp), intent(in) :: lGrad, lHess
    end subroutine Rysg2
  end interface
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! PROLOGUE
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  nGr = 0
  ABeq = EQ(Coor(1,1),Coor(1,2))
  CDeq = EQ(Coor(1,3),Coor(1,4))
  la = iAngV(1)
  lb = iAngV(2)
  lc = iAngV(3)
  ld = iAngV(4)
  ldot2 = ldot
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
  LmbdT = 0
  nijkl = iBasi*jBasj*kBask*lBasl
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
# ifdef _WARNING_WORKAROUND_
  ! Avoid some warnings about unset output arguments
  Temp(nTemp) = One
# endif

  iuvwx(1) = dc(iStb)%nStab
  iuvwx(2) = dc(jStb)%nStab
  iuvwx(3) = dc(kStb)%nStab
  iuvwx(4) = dc(lStb)%nStab
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !           - - - - - - END PROLOGUE - - - - - -
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Find the Double Coset Representatives for center A and B

  if (nIrrep == 1) then
    nDCRR = 1
    iDCRR(0) = 0
    LmbdR = 1
  else
    call DCR(LmbdR,dc(iStb)%iStab,dc(iStb)%nStab,dc(jStb)%iStab,dc(jStb)%nStab,iDCRR,nDCRR)
  end if
  u = real(dc(iStb)%nStab,kind=wp)
  v = real(dc(jStb)%nStab,kind=wp)

  ! Find stabilizer for center A and B

  if (nIrrep == 1) then
    lStabM = 1
    iStabM(0) = 0
  else
    call Inter(dc(iStb)%iStab,dc(iStb)%nStab,dc(jStb)%iStab,dc(jStb)%nStab,iStabM,lStabM)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Find the Double Coset Representatives for center C and D.

  if (nIrrep == 1) then
    nDCRS = 1
    iDCRS(0) = 0
    LmbdS = 1
  else
    call DCR(LmbdS,dc(kStb)%iStab,dc(kStb)%nStab,dc(lStb)%iStab,dc(lStb)%nStab,iDCRS,nDCRS)
  end if
  w = real(dc(kStb)%nStab,kind=wp)
  x = real(dc(lStb)%nStab,kind=wp)

  ! Find stabilizer for center C and D

  if (nIrrep == 1) then
    lStabN = 1
    iStabN(0) = 0
  else
    call Inter(dc(kStb)%iStab,dc(kStb)%nStab,dc(lStb)%iStab,dc(lStb)%nStab,iStabN,lStabN)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Find the Double Coset Representatives for the two charge distributions.

  if (nIrrep == 1) then
    nDCRT = 1
    iDCRT(0) = 0
    LmbdT = 1
  else
    call DCR(LmbdT,iStabM,lStabM,iStabN,lStabN,iDCRT,nDCRT)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Factor due to summation over DCR

  if (MolWgh == 1) then
    Fact = real(nIrrep,kind=wp)/real(LmbdT,kind=wp)
  else if (MolWgh == 0) then
    Fact = u*v*w*x/real(nIrrep**3*LmbdT,kind=wp)
  else
    Fact = sqrt(u*v*w*x)/real(nIrrep*LmbdT,kind=wp)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  nOp(1) = NrOpr(0)
  CoorM(:,1) = Coor(:,1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! - - - - Loop over first set
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do lDCRR=0,nDCRR-1
    nOp(2) = NrOpr(iDCRR(lDCRR))
    call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
    AeqB = EQ(CoorM(1,1),CoorM(1,2))
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! - - - - Loop over second set
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    do lDCRS=0,nDCRS-1
      Coorm(:,3) = Coor(:,3)
      call OA(iDCRS(lDCRS),Coor(1:3,4),CoorM(1:3,4))
      CeqD = EQ(Coor(1,3),CoorM(1,4))
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! - - - - Loop over third set
      !                                                                *
      !*****************************************************************
      !                                                                *
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

        !--------------------------------------------------------------*
        !
        ! Loops to partition the primitives
        !
        !--------------------------------------------------------------*
        lDCR1 = NrOpr(iDCRR(lDCRR))+1
        lDCR2 = NrOpr(iDCRS(lDCRS))+1
        ix2 = iPhase(1,iDCRT(lDCRT))
        iy2 = iPhase(2,iDCRT(lDCRT))
        iz2 = iPhase(3,iDCRT(lDCRT))

        call c_f_pointer(c_loc(Data1(ip_IndZ(1,nZeta),lDCR1)),iData1,[nZeta+1])
        call c_f_pointer(c_loc(Data2(ip_IndZ(1,nEta),lDCR2)),iData2,[nEta+1])
        nZeta_Tot = iData1(nZeta+1)
        nEta_Tot = iData2(nEta+1)

        no_integrals = .true.
        first = .true.
        nGr = 0
        do iZeta=1,nZeta_Tot,IncZet
          mZeta = min(IncZet,nZeta_Tot-iZeta+1)
          ! Check that subblock of contraction matrix has non-zero elements.
          if (lEmpty(Coeff2,nBeta,nBeta,jBasj)) cycle
          do iEta=1,nEta_Tot,IncEta
            mEta = min(IncEta,nEta_Tot-iEta+1)
            ! Check that subblock of contraction matrix has non-zero elements.
            if (lEmpty(Coeff4,nDelta,nDelta,lBasl)) cycle
            Pren = Pren+real(mab*mcd*mZeta*mEta,kind=wp)
            !----------------------------------------------------------*
            !
            ! Fix the control matrixes for derivatives
            ! and try to use translation invariance as
            ! efficient as possible.
            !
            ! OBS DETTA SKALL FLYTTAS UT UR INRE LOOPEN
            !
            !----------------------------------------------------------*
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
            !----------------------------------------------------------*
            !     PRE PRESCREENING                                     *
            !----------------------------------------------------------*

            lZeta = mZeta
            lEta = mEta

            ! Decontract the 2nd order density matrix

            ! Work4->Work2  Work3:scratch
            call Timing(dum1,Time,dum2,dum3)
            if (ldot2) call Tcrtnc_h(Coeff1,nAlpha,iBasi,Coeff2,nBeta,jBasj,Coeff3,nGamma,kBask,Coeff4,nDelta,lBasl,Work4,mab*mcd, &
                                     Work3,nWork3/2,Work2,iData1(iZeta:iZeta+mZeta-1),mZeta,iData2(iEta:iEta+mEta-1),mEta)
            call Timing(dum1,Time,dum2,dum3)
            CPUStat(nTwoDens) = CPUStat(nTwoDens)+Time

            ! Transfer k2 data and prescreen

            ! Work2:PAO-> Work2
            ! Work3 Scratch
            call Timing(dum1,Time,dum2,dum3)
            call Screen_mck(Work2,Work3,mab*mcd,nZeta,nEta,mZeta,mEta,lZeta,lEta,Zeta,ZInv,P,xA,xB,rKab, &
                            Data1(ip_Z(iZeta,nZeta),lDCR1),iData1(iZeta:iZeta+mZeta-1),Data1(ip_abMax(nZeta),ldcr1),Eta,EInv,Q,xG, &
                            xD,rKcd,Data2(ip_Z(iEta,nEta),lDCR2),iData2(iEta:iEta+mEta-1),Data2(ip_abMax(nEta),ldcr2),xpre,1,1,1, &
                            ix2,iy2,iz2,CutInt,PreScr,IndZet,IndEta,ldot2)
            call Timing(dum1,Time,dum2,dum3)
            CPUStat(nScreen) = CPUStat(nScreen)+Time

            Prem = Prem+real(mab*mcd*lZeta*lEta,kind=wp)
            if (lZeta*lEta /= 0) no_integrals = .false.
            if (lZeta*lEta == 0) cycle

            ! Compute integral derivative and accumulate
            ! contribution to the molecular gradient.

            ! Work2:PAO
            ! Work3:Work area  The PO integrals are stored in the begining of Work3

            call Timing(dum1,Time,dum2,dum3)

            call Rysg2(iAngV,nRys,lZeta*lEta,xA,xB,xG,xD,Zeta,ZInv,lZeta,Eta,EInv,lEta,P,nZeta,Q,nEta,CoorM,CoorM,CoorAC,Work3, &
                       nWork3,TERI1,ModU2,Cff2D,Work2,mab*mcd,Hess,nHess,JfGrd,JndGrd,JfHss,JndHss,nOp,iuvwx,IfG,nGr,Indx,lgrad, &
                       ldot,Tr)
            call Timing(dum1,Time,dum2,dum3)
            CPUStat(nIntegrals) = CPUStat(nIntegrals)+Time

            ! Work3 AO
            ! Work3_3  Scratch
            ! ->    Work3_2
            !----------------------------------------------------------*
            !
            ! Transform integrals to AO base
            !
            !----------------------------------------------------------*
            ip2 = nGr*mab*mcd*lZeta*lEta+1
            call Timing(dum1,Time,dum2,dum3)
            call Cntrct_mck(First,Coeff1,nAlpha,iBasi,Coeff2,nBeta,jBasj,Coeff3,nGamma,kBask,Coeff4,nDelta,lBasl,Work3, &
                            nGr*mab*mcd,Work3(ip2),nwork3-ip2,xpre,WorkX,nWorkX,lZeta*lEta,IndZet,nZeta,lZeta,IndEta,nEta,lEta)
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

        if (MolWgh == 1) then
          FactNd = real(nIrrep,kind=wp)/real(LmbdT,kind=wp)
        else if (MolWgh == 0) then
          FactNd = u*v*w*x/real(nIrrep**3*LmbdT,kind=wp)
        else
          factNd = sqrt(u*v*w*x)/real(nirrep*lmbdt,kind=wp)
        end if

        if (FactNd /= One) then
          n = nGr*mab*mcd*nijkl
          WorkX(1:n) = FactNd*WorkX(1:n)
        end if

        !--------------------------------------------------------------*
        !
        !     Transpose abcd,g,IJKL -> bcd,g,IJKL,A Work3 -> Work3_2
        !
        !--------------------------------------------------------------*

        niag = nijkl*nTri_Elem1(lb)*mcd*nGr
        call CrSph_mck(WorkX,niag,nTri_Elem1(la),RSph(ipSph(la)),la,Shells(iShlla)%Transf,Shells(iShlla)%Prjct,Work3,iCmpa)
        nw3 = niag*iCmpa
        ip2 = 1+nw3

        !--------------------------------------------------------------*
        !
        !     Transpose   bcd,g,IJKL,A -> cd,g,IJKL,AB Work3_2->Work3
        !
        !--------------------------------------------------------------*

        niag = nijkl*mcd*nGr*iCmpa
        nw3_2 = niag*jCmpb
        if (nw3+nw3_2 > nWork3) then
          write(u6,*) '1: nw3+nw3_2 > nWork3'
          call Abend()
        end if
        call CrSph_mck(Work3,niag,nTri_Elem1(lb),RSph(ipSph(lb)),lb,Shells(jShllb)%Transf,Shells(jShllb)%Prjct,Work3(ip2),jCmpb)

        !--------------------------------------------------------------*
        !
        !     Transpose  cd,g,IJKL,AB -> d,g,IJKL,ABC  Work3->Work3_2
        !
        !--------------------------------------------------------------*

        niag = nijkl*nGr*nTri_Elem1(ld)*iCmpa*jCmpb
        call CrSph_mck(Work3(ip2),niag,nTri_Elem1(lc),RSph(ipSph(lc)),lc,Shells(kShllc)%Transf,Shells(kShllc)%Prjct,Work3,kCmpc)
        if (niag*kCmpc > nw3) then
          write(u6,*) 'niag*kCmpc > nw3'
          call Abend()
        end if
        nw3 = niag*kCmpc
        ip2 = nw3+1

        !--------------------------------------------------------------*
        !
        !     Transpose   d,g,IJKL,ABC -> g,IJKL,ABCD Work3_2->Work3
        !
        !--------------------------------------------------------------*

        niag = nijkl*nGr*iCmpa*jCmpb*kCmpc
        nw3_2 = niag*lCmpd
        if (nw3+nw3_2 > nWork3) then
          write(u6,*) '2: nw3+nw3_2 > nWork3'
          call Abend()
        end if
        call CrSph_mck(Work3,niag,nTri_Elem1(ld),RSph(ipSph(ld)),ld,Shells(lShlld)%Transf,Shells(lShlld)%Prjct,Work3(ip2),lCmpd)

        !--------------------------------------------------------------*
        !
        !     Transpose g,IJKL,ABCD -> IJKL,ABCD,g Work3->Buffer
        !
        !--------------------------------------------------------------*

        niag = nijkl*iCmpa*jCmpb*kCmpc*lCmpd
        call DGetMO(Work3(ip2),nGr,nGr,niag,Fin,niag)

        ! DEBUG  (calculates gradient from transformed integrals)

        call Timing(dum1,Time,dum2,dum3)
        CPUStat(nTrans) = CPUStat(nTrans)+Time

        !--------------------------------------------------------------*
        !
        !     Send the integrals to clrbuffer for construction of
        !
        !--------------------------------------------------------------*

        call ClrBuf(idcrr(ldcrr),idcrs(ldcrs),idcrt(ldcrt),nGr,Shijij,iAngV,iCmpi,iCmp,iShll,iShell,iShell,iBasi,jBasj,kBask, &
                    lBasl,Dij1,Dij2,mDij,nDij,Dkl1,Dkl2,mDkl,nDkl,Dik1,Dik2,mDik,nDik,Dil1,Dil2,mDil,nDil,Djk1,Djk2,mDjk,nDjk, &
                    Djl1,Djl2,mDjl,nDjl,fin,nfin,Temp(ipFT),nFT,Temp(ipS1),nS1,Temp(ipS2),nS2,Temp(ipTemp),nTe,TwoHam,nTwo2, &
                    JndGrd,Indx,iao,iaost,iuvwx,n8,ltri,moip,nAcO,rMoin,nmoin,ntemp,Buffer,nOp,Din,Dan,new_fock)

      end do

    end do

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  return

end subroutine Twoel_Mck_Internal

end subroutine Twoel_Mck
