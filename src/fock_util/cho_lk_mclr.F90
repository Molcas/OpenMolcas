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
! Copyright (C) Mickael G. Delcey                                      *
!***********************************************************************

subroutine CHO_LK_MCLR(DLT,DI,DA,G2,Kappa,JI,KI,JA,KA,FkI,FkA,MO_Int,QVec,Ash,CMO,CMO_Inv,nOrb,nAsh,DoAct,Fake_CMO2,LuAChoVec, &
                       LuIChoVec,iAChoVec)
!*********************************************************************
!  Author : M. G. Delcey based on cho_LK_rassi_x
!
!  Note:  this routine differs from CHO_LK_RASSI_X because it can
!         handle inactive or active matrix and 1-index transformed
!         densities
!
! *************** INACTIVE AO-BASIS FOCK MATRIX **********************
!
!   FI(ab) = 2 * sum_J  Lab,J * U(J)  -  sum_Jk  Yka,J * Xkb,J
!
!      U(J) = sum_gd  Lgd,J * DI(gd)
!
!      a,b,g,d:  AO-index
!      k:        MO-index   belonging to (Inactive)
!      v,w,x,y:  MO-indices belonging to (Active)
!
!*********************************************************************

use ChoArr, only: nBasSh, nDimRS
use ChoSwp, only: IndRed, InfVec, nnBstRSh
use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Fock_util_global, only: Deco, dmpk, Estimate, Nscreen, Update
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, G2_Type, L_Full_Type, Lab_Type, NDSBA_Type, SBA_Type
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
type(DSBA_Type), intent(in) :: DLT, DA, Kappa, Ash(2), CMO, CMO_Inv
type(DSBA_Type), intent(inout) :: DI, JI(1), KI, JA, KA
real(kind=wp), intent(in) :: G2(*)
type(DSBA_Type), intent(_OUT_) :: FkI, FkA, QVec
real(kind=wp), intent(inout) :: MO_Int(*)
integer(kind=iwp), intent(in) :: nOrb(8), nAsh(8), LuAChoVec(8), LuIChoVec(8), iAChoVec
logical(kind=iwp), intent(in) :: DoAct, Fake_CMO2
#include "warnings.h"
#include "cholesky.fh"
#include "choorb.fh"
integer(kind=iwp) :: i, ia, iab, iabg, iAdr, iAdr2, iag, iaSh, iaSkip, iASQ(8,8,8), ib, iBatch, ibcount, ibg, ibs, ibSh, ibSkip, &
                     iCase, iE, iij, ijS, ijsym, ik, ikl, iLoc, iml, Inc, ioff, ioffa, iOffAB, ioffb, iOffShb, ipG, irc, ired1, &
                     IREDC, iS, ish, iShp, iSwap, ISYM, iSyma, iSymb, iSymv, isymx, iSymy, iTmp, IVEC2, iVrs, jab, jAsh, jaSkip, &
                     jDen, jK, jK_a, jml, jmlmax, JNUM, jOffAB, JRED, JRED1, JRED2, jrs, jS, jsym, jvc, JVEC, k, kaOff(8), kAsh, &
                     kDen, kMOs, kOff(8), krs, kS, kscreen, kSym, l, l1, lAsh, LFMAX, LFULL(2), LKsh, LKshp, LREAD, ls, lSh, lSym, &
                     lvec, LWORK, MaxAct, MaxB, MaxRedT, MaxVecPerBatch, Mmax, mrs, mSh, mTvec, mTvec1, MUSED, MxB, MxBasSh, &
                     n1, n2, nA2, NAv, NAw, Nax, Nay, nBatch, nBsa, nChMo(8), nDen, nMat, nMOs, nnA, nnO, nnShl_2, nRS, NumCV, &
                     numSh1, numSh2, NUMV, NumVT, nVec, nVrs
real(kind=wp) :: Fac1, Fac2, Fact, LKThr, mydmpk, norm, SKsh, tact(2), tau, TCC1, TCC2, TCINT1, TCINT2, TCINT3, TCINT4, tcoul(2), &
                 TCR1, TCR2, TCS1, TCS2, TCT1, TCT2, TCX1, TCX2, Temp, texch(2), thrv(2), tint1(2), tint2(2), tint3(2), tintg(2), &
                 tmotr(2), TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tQmat(2), tread(2), tscrn(2), TWC1, TWC2, &
                 TWINT1, TWINT2, TWINT3, TWINT4, TWR1, TWR2, TWS1, TWS2, TWT1, TWT2, TWX1, TWX2, xFab, xtau(2), xTmp, YMax, YshMax
logical(kind=iwp) :: add, DoScreen, ReadInter, timings
integer(kind=iwp), save :: nVec_
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: Debug
#endif
character(len=50) CFmt
character :: mode, mode2
type(DSBA_Type) :: CM(2), JALT(1), QTmp(2), Tmp(2)
type(SBA_Type) :: Lpq(3)
type(NDSBA_Type) :: DiaH
type(G2_Type) :: MOScr
type(L_Full_Type) :: L_Full
type(Lab_Type) :: Lab
integer(kind=iwp), allocatable :: Indx(:,:,:), iShp_rs(:), kOffSh(:,:), nnBfShp(:,:)
real(kind=wp), allocatable :: AbsC(:), Diag(:), Drs(:), Faa(:), Fia(:), Frs(:,:), Lrs(:,:), MLk(:,:,:), SumAClk(:,:,:), &
                              SvShp(:,:), VJ(:), Ylk(:,:,:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: myJRED1, NNBSTMX, ntv0
real(kind=wp), allocatable :: DiagJ(:)
#endif
real(kind=wp), parameter :: FactCI = -Two, FactXI = Half
logical(kind=iwp), parameter :: DoRead = .false.
character(len=*), parameter :: SECNAM = 'CHO_LK_MCLR'
integer(kind=iwp), external :: Cho_LK_MaxVecPerBatch
real(kind=wp), external :: Cho_LK_ScreeningThreshold, ddot_

!***********************************************************************
#ifdef _DEBUGPRINT_
Debug = .false. ! to avoid double printing in CASSCF-debug
#endif

! Allow LT-format access to JA although it is in SQ-format
call Allocate_DT(JALT(1),nBas,nBas,nSym,aCase='TRI',Ref=JA%A0)
timings = .false.

IREDC = -1 ! unknown reduced set in core

nDen = 2 ! the two bi-orthonormal sets of orbitals
if (Fake_CMO2) nDen = 1 ! MO1 = MO2
kDen = nDen

call CWTIME(TOTCPU1,TOTWALL1) ! start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = zero ! time read/transform vectors
tcoul(:) = zero ! time for computing Coulomb
texch(:) = zero ! time for computing Exchange
tintg(:) = zero ! time for computing (tw|xy) integrals
tmotr(:) = zero ! time for the half-transf of vectors
tscrn(:) = zero ! time for screening overhead
tint1(:) = zero
tint2(:) = zero
tint3(:) = zero
tQmat(:) = zero
tact(:) = zero

! ==================================================================

! Various offsets
! ----------------
nnO = 0
nnA = 0
nA2 = 0
kOff(1) = 0
kAOff(1) = 0
MaxB = nBas(1)
MaxAct = nAsh(1)
do ISYM=2,NSYM
  MaxB = max(MaxB,nBas(iSym))
  MaxAct = max(MaxAct,nAsh(iSym))
  nnO = nnO+nOrb(iSym-1)
  nnA = nnA+nAsh(iSym-1)
  nA2 = nA2+nAsh(iSym-1)**2
  kOff(iSym) = nnO
  kAOff(iSym) = nnA
end do
nnO = nnO+nOrb(nSym)
nnA = nnA+nAsh(nSym)
nA2 = nA2+nAsh(nSym)**2

nChMO(:) = nOrb(:)

if (DoAct) then
  ioff = 0
  do ijsym=1,nsym
    do isym=1,nsym
      jsym = ieor(isym-1,ijsym-1)+1
      do kSym=1,nSym
        lSym = ieor(kSym-1,ijSym-1)+1
        iASQ(isym,jsym,kSym) = ioff
        ioff = ioff+nASh(iSym)*nAsh(jSym)*nASh(kSym)*nASh(lSym)
      end do
    end do
  end do

  ! *** memory for the Q matrices --- temporary array
  call Allocate_DT(QTmp(1),nBas,nAsh,nSym)
  call Allocate_DT(QTmp(2),nBas,nAsh,nSym)
  QTmp(1)%A0(:) = Zero
  QTmp(2)%A0(:) = Zero

  iCase = 1
  call Allocate_DT(MOScr,nAsh,nSym,iCase)
  MOScr%A0(:) = Zero
end if
!*************************************************
if (Deco) then
  call Allocate_DT(CM(1),nBas,nBas,nSym)
  call Allocate_DT(CM(2),nBas,nBas,nSym)

  call Allocate_DT(Tmp(1),nBas,nBas,nSym)
  call Allocate_DT(Tmp(2),nBas,nBas,nSym)

  do iS=1,nSym

    !* Create Cholesky orbitals from DI

    call CD_InCore(DI%SB(iS)%A2,nBas(iS),CM(1)%SB(iS)%A2,nBas(iS),nChMO(iS),1.0e-12_wp,irc)
    if (.not. Fake_CMO2) then

      !* MO transform

      call DGEMM_('T','T',nChMO(iS),nBas(iS),nBas(iS),One,CM(1)%SB(iS)%A2,nBas(iS),CMO_inv%SB(iS)%A2,nBas(iS),Zero, &
                  Tmp(2)%SB(iS)%A2,nChMO(iS))

      !* Create one-index transformed Cholesky orbitals

      call DGEMM_('N','N',nChMO(iS),nBas(iS),nBas(iS),One,Tmp(2)%SB(iS)%A2,nChMO(iS),Kappa%SB(is)%A2,nBas(iS),Zero, &
                  Tmp(1)%SB(iS)%A2,nChMO(iS))

      !* AO transform

      call DGEMM_('N','T',nBas(iS),nChMO(iS),nBas(iS),One,CMO%SB(iS)%A2,nBas(iS),Tmp(1)%SB(iS)%A2,nChMO(iS),Zero,CM(2)%SB(iS)%A2, &
                  nBas(iS))
    end if
  end do
  call Deallocate_DT(Tmp(2))
  call Deallocate_DT(Tmp(1))
else

  write(u6,*) 'Cho_LK_MCLR: this will not work'
  call Abend()

end if

!*************************************************

! Define the max number of vectors to be treated in core at once

MaxVecPerBatch = Cho_LK_MaxVecPerBatch()

! Define the screening threshold

LKThr = Cho_LK_ScreeningThreshold(-One)
mydmpk = dmpk
!mydmpk = Zero

! Vector MO transformation screening thresholds
NumVT = NumChT
call GAIGOP_SCAL(NumVT,'+')
thrv(1) = (LKThr/(max(1,nnO)*NumVT))*mydmpk**2
xtau(1) = sqrt((LKThr/max(1,nnO))*mydmpk)
xtau(2) = xtau(1) ! dummy init

if (.not. Fake_CMO2) then
  norm = sqrt(ddot_(size(Kappa%A0),Kappa%A0,1,Kappa%A0,1))
  xtau(2) = sqrt((LKThr/max(1,nnO))*mydmpk)*norm
  mydmpk = min(norm,mydmpk)
  thrv(2) = (LKThr/(max(1,nnO)*NumVT))*mydmpk**2
end if
tau = (LKThr/max(1,nnO))*mydmpk

MaxRedT = MaxRed
call GAIGOP_SCAL(MaxRedT,'+')

if (Estimate) then
  xtau(1) = xtau(1)/sqrt(real(MaxRedT,kind=wp))
  if (.not. Fake_CMO2) xtau(2) = xtau(2)/sqrt(real(MaxRedT,kind=wp))
  tau = tau/MaxRedT
end if

call mma_allocate(DIAG,NNBSTRT(1),Label='DIAG')

#ifdef _MOLCAS_MPP_
if ((nProcs > 1) .and. Update .and. Is_Real_Par()) then
  NNBSTMX = 0
  do i=1,nSym
    NNBSTMX = max(NNBSTMX,NNBSTR(i,1))
  end do
  call mma_allocate(diagJ,NNBSTMX,Label='diagJ')
  diagJ(:) = Zero
end if
#endif
! *************** Read the diagonal integrals (stored as 1st red set)
if (Update) call CHO_IODIAG(DIAG,2) ! 2 means "read"

! allocate memory for sqrt(D(a,b)) stored in full (squared) dim
call Allocate_DT(DiaH,nBas,nBas,nSym)
DiaH%A0(:) = Zero

! allocate memory for the abs(C(l)[k])
call mma_allocate(AbsC,MaxB,Label='AbsC')

! allocate memory for the Y(l)[k] vectors
call mma_allocate(Ylk,MaxB,nnO,nDen,Label='Ylk')

! allocate memory for the ML[k] lists of largest elements in significant shells
call mma_allocate(MLk,nShell,nnO,nDen,Label='MLk')

! allocate memory for the lists of  S:= sum_l abs(C(l)[k]) for each shell
call mma_allocate(SumAClk,nShell,nnO,nDen)

! allocate memory for the Index arrays
call mma_allocate(Indx,[0,nShell],[1,nnO],[1,nDen],Label='Indx')

! allocate memory for kOffSh
call mma_allocate(kOffSh,nShell,nSym,Label='kOffh')

! allocate memory for nnBfShp
nnShl_2 = nShell**2
call mma_allocate(nnBfShp,nnShl_2,nSym,Label='nnBfShp')

! allocate memory for iShp_rs
call mma_allocate(iShp_rs,nnShl_tot,Label='iShp_rs')

! allocate memory for the shell-pair Frobenius norm of the vectors
call mma_allocate(SvShp,nnShl,2,Label='SvShp')

! *** Compute Shell Offsets ( MOs and transformed vectors)

MxBasSh = 0

do iSyma=1,nSym

  LKsh = 0

  do iaSh=1,nShell ! kOffSh(iSh,iSym)

    kOffSh(iaSh,iSyma) = LKsh

    LKsh = LKsh+nBasSh(iSyma,iaSh)

    MxBasSh = max(MxBasSh,nBasSh(iSyma,iaSh))

  end do

end do

! allocate memory for the diagonal elements of the Fock matrix
call mma_allocate(Fia,MxBasSh,Label='Fia')
call mma_allocate(Faa,nShell,Label='Faa')
Fia(:) = Zero
Faa(:) = Zero

! *** Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)
do jDen=1,nDen
  do kSym=1,nSym

    do jK=1,nOrb(kSym)
      jK_a = jK+kOff(kSym)

      do iaSh=1,nShell

        SKsh = zero
        iS = kOffSh(iaSh,kSym)+1
        iE = kOffSh(iaSh,kSym)+nBasSh(kSym,iaSh)
        do ik=iS,iE
          SKsh = SKsh+CM(jDen)%SB(kSym)%A2(ik,jK)**2
        end do

        SumAClk(iaSh,jk_a,jDen) = Sksh

      end do

    end do
  end do
end do

! *** Compute Shell-pair Offsets in the K-matrix

do iSyma=1,nSym

  LKshp = 0

  do iaSh=1,nShell

    do ibSh=1,nShell

      iShp = nShell*(iaSh-1)+ibSh

      nnBfShp(iShp,iSyma) = LKShp

      LKShp = LKShp+nBasSh(iSyma,iaSh)*nBasSh(iSyma,ibSh)

    end do

  end do

end do

! *** Mapping shell pairs from the full to the reduced set

call Mk_iShp_rs(iShp_rs,nShell)

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! *************** BIG LOOP OVER VECTORS SYMMETRY *******************
do jSym=1,nSym

  iAdr = 0

  NumCV = NumCho(jSym)
  call GAIGOP_SCAL(NumCV,'max')
  if (NumCV < 1) cycle

  JNUM = 1
  call Allocate_DT(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym,Memory=LFULL)

  iLoc = 3 ! use scratch location in reduced index arrays

! ****************     MEMORY MANAGEMENT SECTION    *****************
! --------------------------------------------------------------
! compute memory needed to store at least 1 vector of JSYM
! and do all the subsequent calculations
! --------------------------------------------------------------
  mTvec1 = 0
  MxB = 0
  do l=1,nSym
    k = Mul(l,JSYM)
    Mmax = max(0,nOrb(k))
    if (Mmax > 0) MxB = max(MxB,nBas(l))
    if (DoAct) then
      mTvec1 = mTvec1+nAsh(k)*nBas(l)
    end if
  end do

  LFMAX = max(3*mTvec1,LFULL(1)) ! re-use memory for the active vec
  mTvec = nDen*max(MxB,1) ! mem for storing half-transformed vec

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

  JRED1 = InfVec(1,2,jSym)            ! red set of the 1st vec
  JRED2 = InfVec(NumCho(jSym),2,jSym) ! red set of the last vec
# ifdef _MOLCAS_MPP_
  myJRED1 = JRED1 ! first red set present on this node
  ntv0 = 0
# endif

  ! entire red sets range for parallel run
  call GAIGOP_SCAL(JRED1,'min')
  call GAIGOP_SCAL(JRED2,'max')

  kscreen = 1
  DoScreen = .true.

  do JRED=JRED1,JRED2

    call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

    if (nVrs /= 0) then ! otherwise no vectors in that (jred,jsym)

      if (nVrs < 0) then
        write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!',nVrs
        call Abend()
      end if

      call Cho_X_SetRed(irc,iLoc,JRED)
      ! set index arrays at iLoc
      if (irc /= 0) then
        write(u6,*) SECNAM//'cho_X_setred non-zero return code. rc= ',irc
        call Abend()
      end if

      IREDC = JRED

      nRS = nDimRS(JSYM,JRED)

      if (JSYM == 1) then

        call mma_allocate(Drs,nRS,Label='Drs')
        if (DoAct) then
          call mma_allocate(Frs,nRS,2,Label='Frs')
        else
          call mma_allocate(Frs,nRS,1,Label='Frs')
        end if
        Drs(:) = Zero
        Frs(:,:) = Zero

      end if

      call mma_maxDBLE(LWORK)

      nVec = min((LWORK-LFULL(2))/(nRS+mTvec+LFMAX),min(nVrs,MaxVecPerBatch))

      ! Store nVec to make sure the routine always uses the same
      if (iAChoVec == 1) nVec_ = nVec
      ReadInter = (iAChoVec == 2) .and. (nVec == nVec_)
      ! nVec /= nVec_ should happen only if lack of memory
      !ReadInter = .false.

      if (nVec < 1) then
        write(u6,*) SECNAM//': Insufficient memory for batch'
        write(u6,*) 'LWORK= ',LWORK
        write(u6,*) 'min. mem. need= ',nRS+mTvec+LFMAX
        write(u6,*) 'nRS= ',nRS
        write(u6,*) 'mTvec= ',mTvec
        write(u6,*) 'LFMAX= ',LFMAX
        write(u6,*) 'jsym= ',jsym
        call Quit(_RC_MEMORY_ERROR_)
        nBatch = -9999 ! dummy assignment
      end if

      LREAD = nRS*nVec

      call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

      if (JSYM == 1) then
        ! Transform the density to reduced storage
        add = .false.
        nMat = 1
        call swap_full2rs(irc,iLoc,nRS,nMat,JSYM,[DLT],Drs,add)
      end if

      ! BATCH over the vectors ----------------------------

      nBatch = (nVrs-1)/nVec+1

      do iBatch=1,nBatch

        if (iBatch == nBatch) then
          JNUM = nVrs-nVec*(nBatch-1)
        else
          JNUM = nVec
        end if

        JVEC = nVec*(iBatch-1)+iVrs
        IVEC2 = JVEC-1+JNUM

        call CWTIME(TCR1,TWR1)

        call CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,NUMV,IREDC,MUSED)

        if ((NUMV <= 0) .or. (NUMV /= JNUM)) return

        call CWTIME(TCR2,TWR2)
        tread(1) = tread(1)+(TCR2-TCR1)
        tread(2) = tread(2)+(TWR2-TWR1)

        if (JSYM == 1) then
          ! ************ (alpha+beta) COULOMB CONTRIBUTION  ****************
          !
          ! Contraction with the density matrix
          ! -----------------------------------
          ! V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * DI(rs)
          !==========================================================

          call CWTIME(TCC1,TWC1)

          call mma_allocate(VJ,JNUM,Label='VJ')

          call DGEMV_('T',nRS,JNUM,ONE,Lrs,nRS,Drs,1,ZERO,VJ,1)

          ! FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
          !===========================================================

          Fact = real(min(jVec-iVrs,1),kind=wp)

          call DGEMV_('N',nRS,JNUM,FactCI,Lrs,nRS,VJ,1,Fact,Frs(:,1),1)

          call CWTIME(TCC2,TWC2)
          tcoul(1) = tcoul(1)+(TCC2-TCC1)
          tcoul(2) = tcoul(2)+(TWC2-TWC1)

          call mma_deallocate(VJ)

        end if ! Coulomb contribution

        ! *************** EXCHANGE CONTRIBUTIONS  **********************
        !                                                              *
        !***************************************************************
        !***************************************************************
        !***************************************************************
        !                                                              *
        call Allocate_DT(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym)
        call Allocate_DT(Lab,JNUM,nBasSh,nBas,nShell,nSym,nDen)

        call CWTIME(TCS1,TWS1)
        ! -----------------------------------------------------------------
        ! Estimate the diagonals :   D(a,b) = sum_J (Lab,J)^2
        !
        if (Estimate) then

          call Fzero(DIAG(1+iiBstR(jSym,1)),NNBSTR(jSym,1))

          do krs=1,nRS

            mrs = iiBstR(JSYM,iLoc)+krs
            jrs = IndRed(mrs,iLoc) ! address in 1st red set

            do jvc=1,JNUM

              Diag(jrs) = Diag(jrs)+Lrs(krs,jvc)**2

            end do

          end do

        end if

        call CWTIME(TCS2,TWS2)
        tscrn(1) = tscrn(1)+(TCS2-TCS1)
        tscrn(2) = tscrn(2)+(TWS2-TWS1)

        call CWTIME(TCX1,TWX1)

        ! Reorder vectors to Full-dimensions
        !
        ! Vectors are returned in the storage LaJ,b with the restriction:
        !
        !    Sym(a) >= Sym(b)
        !
        ! and blocked in shell pairs

        call CHO_getShFull(Lrs,lread,JNUM,JSYM,IREDC,L_Full,SvShp,nnShl,iShp_rs,nnShl_tot)

        call CWTIME(TCX2,TWX2)
        texch(1) = texch(1)+(TCX2-TCX1)
        texch(2) = texch(2)+(TWX2-TWX1)

        if (DoScreen) then

          call CWTIME(TCS1,TWS1)

          ! Compute DH(a,b)=sqrt(D(a,b)) from the updated diagonals.
          !                              Only the symmetry blocks with
          !                              compound symmetry JSYM are computed
          ! ------------------------------------------------------------
          ired1 = 1 ! location of the 1st red set
          call swap_tosqrt(irc,ired1,NNBSTRT(1),JSYM,DIAH,DIAG)

          call CWTIME(TCS2,TWS2)
          tscrn(1) = tscrn(1)+(TCS2-TCS1)
          tscrn(2) = tscrn(2)+(TWS2-TWS1)

        end if

        do kSym=1,nSym

          lSym = Mul(JSYM,kSym)

          do jK=1,nOrb(kSym)

            jk_a = kOff(kSym)+jK

            Lab%A0(1:nDen*nBas(lSym)*JNUM) = Zero

            if (DoScreen) then

              call CWTIME(TCS1,TWS1)
              !--------------------------------------------------------------
              ! Setup the screening
              !--------------------------------------------------------------

              do jDen=1,nDen

                do ik=1,nBas(kSym)
                  AbsC(ik) = abs(CM(jDen)%SB(kSym)%A2(ik,jK))
                end do

                if (lSym >= kSym) then
                  ! ----------------------------------------------------------
                  ! Y(l)[k] = sum_n  DH(l,n) * |C(n)[k]|
                  !===========================================================
                  Mode = 'N'
                  n1 = nBas(lSym)
                  n2 = nBas(kSym)

                else
                  ! ----------------------------------------------------------
                  ! Y(l)[k] = sum_n  DH(n,l) * |C(n)[k]|
                  !===========================================================
                  Mode = 'T'
                  n1 = nBas(kSym)
                  n2 = nBas(lSym)

                end if

                if (n1 > 0) call DGEMV_(Mode,n1,n2,ONE,DiaH%SB(lSym,kSym)%A2,n1,AbsC,1,ZERO,Ylk(1,jK_a,jDen),1)

              end do

              ! List the shells present in Y(l)[k] by the largest element
              do jDen=1,nDen
                do ish=1,nShell
                  YshMax = zero
                  do ibs=1,nBasSh(lSym,ish)
                    YshMax = max(YshMax,Ylk(koffSh(ish,lSym)+ibs,jK_a,jDen))
                  end do
                  MLk(ish,jK_a,jDen) = YshMax
                end do
              end do

              ! Sort the lists ML[k]
              do jDen=1,nDen
                do ish=1,nShell
                  Indx(iSh,jk_a,jDen) = ish
                end do
              end do

              ! ****  The Max in the MO set 1 is used as reference
              numSh1 = 0 ! # of significant shells in MO set 1
              YMax = MLk(1,jK_a,1)
              jmlmax = 1
              do iml=2,nShell ! get the max in the MO set 1
                if (MLk(iml,jK_a,1) > YMax) then
                  YMax = MLk(iml,jK_a,1)
                  jmlmax = iml
                end if
              end do
              if (jmlmax /= 1) then ! swap positions
                xTmp = MLk(1,jK_a,1)
                iTmp = Indx(1,jk_a,1)
                MLk(1,jK_a,1) = YMax
                Indx(1,jk_a,1) = Indx(jmlmax,jk_a,1)
                MLk(jmlmax,jK_a,1) = xTmp
                Indx(jmlmax,jk_a,1) = iTmp
              end if

              ! **** Sort the list for the MO set 2   iff  MOs1 /= MOs2
              if (.not. Fake_CMO2) then
                numSh2 = 0 ! # of significant shells in MO set 2
                jml = 1
                do while (jml <= nShell)

                  YMax = MLk(jml,jK_a,2)
                  jmlmax = jml

                  do iml=jml+1,nShell ! get the max
                    if (MLk(iml,jK_a,2) > YMax) then
                      YMax = MLk(iml,jK_a,2)
                      jmlmax = iml
                    end if
                  end do

                  if (jmlmax /= jml) then ! swap positions
                    xTmp = MLk(jml,jK_a,2)
                    iTmp = Indx(jml,jk_a,2)
                    MLk(jml,jK_a,2) = YMax
                    Indx(jml,jk_a,2) = Indx(jmlmax,jk_a,2)
                    MLk(jmlmax,jK_a,2) = xTmp
                    Indx(jmlmax,jk_a,2) = iTmp
                  end if

                  ! Exact bounds (quadratic scaling of the MO transformation)
                  ! Note that in true RASSI the exchange matrix is not
                  ! positive definite.

                  if (MLk(jml,jK_a,2) >= xtau(2)) then
                    numSh2 = numSh2+1
                  else
                    jml = nShell ! exit the loop
                  end if

                  jml = jml+1

                end do

                Indx(0,jk_a,2) = numSh2
                numSh1 = 1

              else ! fake biorthonormal basis

                numSh2 = 6669666 ! dummy assignement

                if (MLk(1,jK_a,1) >= xtau(1)) numSh1 = 1

              end if

              ! **** Sort the list for the MO set 1 only if needed
              if (numSh2 > 0) then
                jml = 2 ! the 1st element has already been treated
                do while (jml <= nShell)

                  YMax = MLk(jml,jK_a,1)
                  jmlmax = jml
                  do iml=jml+1,nShell ! get the max
                    if (MLk(iml,jK_a,1) > YMax) then
                      YMax = MLk(iml,jK_a,1)
                      jmlmax = iml
                    end if
                  end do

                  if (jmlmax /= jml) then ! swap positions
                    xTmp = MLk(jml,jK_a,1)
                    iTmp = Indx(jml,jk_a,1)
                    MLk(jml,jK_a,1) = YMax
                    Indx(jml,jk_a,1) = Indx(jmlmax,jk_a,1)
                    MLk(jmlmax,jK_a,1) = xTmp
                    Indx(jmlmax,jk_a,1) = iTmp
                  end if

                  if ((.not. Fake_CMO2) .and. (MLk(jml,jK_a,1) >= xtau(1))) then
                    numSh1 = numSh1+1

                  ! Here we use a non-exact bound for the exchange matrix because a
                  ! fake rassi (MOs1=MOs2) has a positive definite exchange
                  else if (Fake_CMO2 .and. (MLk(jml,jK_a,1) >= xtau(1))) then
                    numSh1 = numSh1+1
                  else
                    jml = nShell ! exit the loop
                  end if

                  jml = jml+1

                end do
              else
                numSh1 = 0
              end if

              Indx(0,jk_a,1) = numSh1

              call CWTIME(TCS2,TWS2)
              tscrn(1) = tscrn(1)+(TCS2-TCS1)
              tscrn(2) = tscrn(2)+(TWS2-TWS1)
              !------------------------------------------------------------------
            end if ! Screening setup

            ! Transform vectors for shells in the lists ML[k]
            !
            ! Screening based on the Frobenius norm: sqrt(sum_ij A(i,j)^2)
            !
            !  || La,J[k] ||  <=  || Lab,J || * || Cb[k] ||

            call CWTIME(TCT1,TWT1)

            do jDen=1,nDen

              !MGD maybe check that NumSh1 matches the one on the record?

              do iSh=1,Indx(0,jk_a,jDen)

                iaSh = Indx(iSh,jk_a,jDen)

                Lab%Keep(iaSh,jDen) = .true.

                ibcount = 0

                do ibSh=1,nShell

                  iOffShb = kOffSh(ibSh,kSym)

                  iShp = iTri(iaSh,ibSh)

                  if (iShp_rs(iShp) <= 0) cycle

                  if (lSym >= kSym) then

                    if ((nnBstRSh(JSym,iShp_rs(iShp),iLoc)*nBasSh(lSym,iaSh)*nBasSh(kSym,ibSh) > 0) .and. &
                        (abs(SumAClk(ibSh,jK_a,jDen)*SvShp(iShp_rs(iShp),1)) >= thrv(jDen))) then

                      if (ReadInter .and. (jDen == 1)) then

                        !* Read vectors

                        ibcount = 1

                        exit
                      else

                        !* Or compute them

                        ibcount = ibcount+1

                        l1 = 1
                        if (iaSh < ibSh) l1 = 2

                        ! LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
                        ! ----------------------------------

                        call DGEMV_('N',nBasSh(lSym,iaSh)*JNUM,nBasSh(kSym,ibSh),ONE,L_Full%SPB(lSym,iShp_rs(iShp),l1)%A21, &
                                    nBasSh(lSym,iaSh)*JNUM,CM(jDen)%SB(kSym)%A2(iOffShb+1:,jK),1,ONE,Lab%SB(iaSh,lSym,jDen)%A,1)

                      end if

                    end if

                  else ! lSym < kSym
                    if ((nnBstRSh(JSym,iShp_rs(iShp),iLoc)*nBasSh(lSym,iaSh)*nBasSh(kSym,ibSh) > 0) .and. &
                        (abs(SumAClk(ibSh,jK_a,jDen)*SvShp(iShp_rs(iShp),1)) >= thrv(jDen))) then
                      ibcount = ibcount+1

                      l1 = 1
                      if (ibSh < iaSh) l1 = 2

                      ! LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
                      ! ----------------------------------

                      call DGEMV_('T',nBasSh(kSym,ibSh),JNUM*nBasSh(lSym,iaSh),One,L_Full%SPB(kSym,iShp_rs(iShp),l1)%A12, &
                                  nBasSh(kSym,ibSh),CM(jDen)%SB(kSym)%A2(iOffShb+1:,jK),1,ONE,Lab%SB(iaSh,lSym,jDen)%A,1)

                    end if
                  end if

                end do ! ibSh

                if (lSym >= kSym) then

                  if (ReadInter .and. (jDen == 1)) then

                    if (ibcount > 0) then
                      lvec = nBasSh(lSym,iaSh)*JNUM
                      call DDAFILE(LuIChoVec(Jsym),2,Lab%SB(iaSh,lSym,jDen)%A,lvec,iAdr)
                    end if

                  else

                    !* Store vectors

                    if ((iAChoVec == 1) .and. (jDen == 1) .and. (ibcount > 0)) then
                      lvec = nBasSh(lSym,iaSh)*JNUM
                      call DDAFILE(LuIChoVec(Jsym),1,Lab%SB(iaSh,lSym,jDen)%A,lvec,iAdr)

                    end if
                  end if
                end if

                ! The following re-assignement is used later on to check if the
                ! iaSh vector LaJ[k] can be neglected because identically zero

                if (ibcount == 0) Lab%Keep(iaSh,jDen) = .false.

              end do ! iSh

              do iSh=Indx(0,jk_a,jDen)+1,nshell
                iaSh = Indx(iSh,jk_a,jDen)
                Lab%Keep(iaSh,jDen) = .false.
              end do

            end do ! jDen

            call CWTIME(TCT2,TWT2)
            tmotr(1) = tmotr(1)+(TCT2-TCT1)
            tmotr(2) = tmotr(2)+(TWT2-TWT1)

            ! Prepare the J-screening

            call CWTIME(TCS1,TWS1)

            do iSh=1,Indx(0,jk_a,1)

              iaSh = Indx(iSh,jk_a,1)

              iaSkip = merge(1,0,Lab%Keep(iaSh,1))

              jaSkip = merge(1,0,Lab%Keep(iaSh,kDen))

              if (iaSkip*jaSkip == 0) cycle

              if (lSym >= kSym) then

                ! Faa,[k] = sum_J  LaJ[k2]*LaJ[k1]
                ! --------------------------------
                Inc = nBasSh(lSym,iaSh)
                n1 = 1

              else ! lSym < kSym

                ! Faa,[k] = sum_J  LJa[k2]*LJa[k1]
                ! --------------------------------
                Inc = 1
                n1 = JNUM

              end if

              Temp = Zero
              do ia=1,nBasSh(lSym,iaSh)
                Fia(ia) = DDot_(JNUM,Lab%SB(iaSh,lSym,kDen)%A(1+n1*(ia-1):),Inc,Lab%SB(iaSh,lSym,1)%A(1+n1*(ia-1):),Inc)
                Temp = max(abs(Fia(ia)),Temp)
              end do

              Faa(iaSh) = Temp

            end do

            call CWTIME(TCS2,TWS2)
            tscrn(1) = tscrn(1)+(TCS2-TCS1)
            tscrn(2) = tscrn(2)+(TWS2-TWS1)

            !--------------------------------------------------------
            ! Compute exchange matrix for the interacting shell pairs
            !--------------------------------------------------------

            call CWTIME(TCX1,TWX1)

            do lSh=1,Indx(0,jk_a,1)

              iaSh = Indx(lSh,jk_a,1)

              iaSkip = merge(1,0,Lab%Keep(iaSh,kDen))

              mSh = 1

              do while (mSh <= Indx(0,jk_a,kDen))

                ibSh = Indx(mSh,jk_a,kDen)

                ibSkip = merge(1,0,Lab%Keep(ibSh,1))

                iShp = nShell*(iaSh-1)+ibSh

                iOffShb = kOffSh(ibSh,lSym)

                iOffAB = nnBfShp(iShp,lSym)

                xFab = sqrt(abs(Faa(iaSh)*Faa(ibSh)))

                if (MLk(lSh,jK_a,1)*MLk(mSh,jK_a,kDen) < tau) then

                  mSh = Indx(0,jk_a,kDen) ! skip rest

                else if ((xFab >= tau/MaxRedT) .and. (iaSkip*ibSkip == 1)) then

                  nBsa = nBasSh(lSym,iaSh)
                  if (lSym >= kSym) then

                    !  F(a,b)[k] = F(a,b)[k] + FactXI * sum_J  X2(a,J)[k] * X1(b,J)[k]
                    ! ----------------------------------------------------------------
                    n1 = nBasSh(lSym,iaSh)
                    n2 = nBasSh(lSym,ibSh)
                    Mode = 'N'
                    Mode2 = 'T'

                  else ! lSym < kSym

                    !  F(a,b)[k] = F(a,b)[k] + FactXI * sum_J  X2(J,a)[k] * X1(J,b)[k]
                    ! ----------------------------------------------------------------
                    n1 = JNUM
                    n2 = JNUM
                    Mode = 'T'
                    Mode2 = 'N'

                  end if

                  call DGEMM_(Mode,Mode2,nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,FactXI,Lab%SB(iaSh,lSym,kDen)%A,n1, &
                              Lab%SB(ibSh,lSym,1)%A,n2,ONE,KI%SB(lSym)%A1(1+iOffAB:),nBsa)

                end if

                mSh = mSh+1 ! update shell counter

              end do

            end do

            call CWTIME(TCX2,TWX2)
            texch(1) = texch(1)+(TCX2-TCX1)
            texch(2) = texch(2)+(TWX2-TWX1)

          end do ! loop over k MOs

        end do ! loop over MOs symmetry

        call Deallocate_DT(Lab)
        call Deallocate_DT(L_Full)
        !                                                              *
        !***************************************************************
        !***************************************************************
        !***************************************************************
        !                                                              *
        DoScreen = .false. ! avoid redo screening inside batch loop

        ! Diagonals updating. It only makes sense if Nscreen > 0

        if (Update .and. (Nscreen > 0)) then

          call CWTIME(TCS1,TWS1)
          ! -----------------------------------------------------------------
          ! update the diagonals :   D(a,b) = D(a,b) - sum_J (Lab,J)^2
          !
          ! subtraction is done in the 1st reduced set
#         ifdef _MOLCAS_MPP_
          if ((nProcs > 1) .and. Is_Real_Par()) then

            do krs=1,nRS

              mrs = iiBstR(JSYM,iLoc)+krs
              jrs = IndRed(mrs,iLoc)-iiBstR(JSYM,1)

              do jvc=1,JNUM

                DiagJ(jrs) = DiagJ(jrs)+Lrs(krs,jvc)**2
              end do

            end do

          else
#         endif

            do krs=1,nRS

              mrs = iiBstR(JSYM,iLoc)+krs
              jrs = IndRed(mrs,iLoc) ! address in 1st red set

              do jvc=1,JNUM

                Diag(jrs) = Diag(jrs)-Lrs(krs,jvc)**2
              end do

            end do

#         ifdef _MOLCAS_MPP_
          end if
#         endif

          call CWTIME(TCS2,TWS2)
          tscrn(1) = tscrn(1)+(TCS2-TCS1)
          tscrn(2) = tscrn(2)+(TWS2-TWS1)

        end if

        ! ************  END EXCHANGE CONTRIBUTION  ****************

        ! ----------------------------------------------------------------
        ! First half Active transformation  Lvb,J = sum_a  C1(v,a) * Lab,J
        ! ----------------------------------------------------------------

        call CWTIME(TCINT1,TWINT1)

        ! Set up the skipping flags
        ! The memory used before for the full-dimension AO-vectors
        !     is now re-used to store half and full transformed
        !     vectors in the active space
        ! ---------------------------------------------------------
        if (DoAct) then

          iSwap = 0  ! Lvb,J are returned
          ! Lvb,J
          ! Lvi,J i general MO index
          ! L~vi,J ~ transformed index
          call Allocate_DT(Lpq(1),nAsh,nBas,nVec,JSYM,nSym,iSwap)
          call Allocate_DT(Lpq(2),nAsh,nAsh,nVec,JSYM,nSym,iSwap)
          call Allocate_DT(Lpq(3),nAsh,nAsh,nVec,JSYM,nSym,iSwap)

          !MGD should we compute only if there are active orbitals in this sym?

          !* Read vectors

          kMOs = 1
          nMOs = 1 ! Active MOs (1st set)
          if (iAChoVec == 2) then
            ioff = 0
            do i=1,nSym
              k = Mul(i,JSYM)
              lvec = nAsh(k)*nBas(i)*JNUM
              iAdr2 = (JVEC-1)*nAsh(k)*nBas(i)+ioff
              call DDAFILE(LuAChoVec(Jsym),2,Lpq(1)%SB(k)%A3,lvec,iAdr2)
              ioff = ioff+nAsh(k)*nBas(i)*NumCho(jSym)
            end do
          else
            ! Lrs * MO
            call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,ASh(1),Lpq,DoRead)
          end if

          if (irc /= 0) return

          !* Store vectors

          ioff = 0
          if (iAChoVec == 1) then
            do i=1,nSym
              k = Mul(i,JSYM)
              lvec = nAsh(k)*nBas(i)*JNUM
              iAdr2 = (JVEC-1)*nAsh(k)*nBas(i)+ioff
              call DDAFILE(LuAChoVec(Jsym),1,Lpq(1)%SB(k)%A3,lvec,iAdr2)
              ioff = ioff+nAsh(k)*nBas(i)*NumCho(jSym)
            end do
          end if
          call CWTIME(TCINT2,TWINT2)
          tint1(1) = tint1(1)+(TCINT2-TCINT1)
          tint1(2) = tint1(2)+(TWINT2-TWINT1)

          ! ----------------------------------------------------------------
          ! Active-Active transformation  Lvw,J = sum_b  Lvb,J * C2(w,b)
          ! ----------------------------------------------------------------
          do iSymb=1,nSym

            iSymv = Mul(JSYM,iSymb)
            NAv = nAsh(iSymv)
            NAw = nAsh(iSymb)

            if (NAv*Naw /= 0) then

              do JVC=1,JNUM

                call CWTIME(TCINT2,TWINT2)
                ! Lv~w
                if (.not. Fake_CMO2) then
                  call DGEMM_('N','T',NAv,Naw,NBAS(iSymb),One,Lpq(1)%SB(iSymv)%A3(:,:,JVC),NAv,Ash(2)%SB(iSymb)%A2,NAw,Zero, &
                              Lpq(3)%SB(iSymv)%A3(:,:,JVC),NAv)
                  call CWTIME(TCINT4,TWINT4)
                  tint1(1) = tint1(1)+(TCINT4-TCINT2)
                  tint1(2) = tint1(2)+(TWINT4-TWINT2)
                  ! ----------------------------------------------------------------
                  ! Formation of the Q matrix Qpx = Lpy Lv~w Gxyvw
                  ! ----------------------------------------------------------------
                  ! ~Lxy=Lv~w Gxyvw
                  !MGD probably additional nSym loop
                  call DGEMV_('N',NAv*Naw,NAv*Naw,ONE,G2,NAv*Naw,Lpq(3)%SB(iSymv)%A3(:,:,JVC),1,ZERO,Lpq(2)%SB(iSymv)%A3(:,:,JVC),1)
                  ! Qpx=Lpy ~Lxy
                  call DGEMM_('T','N',NBAS(iSymb),NAw,Nav,Two,Lpq(1)%SB(iSymv)%A3(:,:,JVC),NAv,Lpq(2)%SB(iSymv)%A3(:,:,JVC),Naw, &
                              ONE,QTmp(2)%SB(iSymb)%A2,NBAS(iSymb))
                  call CWTIME(TCINT3,TWINT3)
                  tQmat(1) = tQmat(1)+(TCINT3-TCINT4)
                  tQmat(2) = tQmat(2)+(TWINT3-TWINT4)
                end if
                call CWTIME(TCINT3,TWINT3)
                ! ----------------------------------------------------------------
                ! Active-Active transformation  Lvw,J = sum_b  Lvb,J * C2(w,b)
                ! ----------------------------------------------------------------
                ! Lvw
                call DGEMM_('N','T',NAv,Naw,NBAS(iSymb),One,Lpq(1)%SB(iSymv)%A3(:,:,JVC),NAv,Ash(1)%SB(iSymb)%A2,NAw,Zero, &
                            Lpq(2)%SB(iSymv)%A3(:,:,JVC),NAv)

                call CWTIME(TCINT2,TWINT2)
                tint1(1) = tint1(1)+(TCINT2-TCINT3)
                tint1(2) = tint1(2)+(TWINT2-TWINT3)
              end do

              ! *************** EVALUATION OF THE (TW|XY) INTEGRALS ***********
              if (iSymv > iSymb) cycle
              do isymx=1,iSymb
                iSymy = Mul(JSYM,iSymx)
                if ((iSymy > iSymx) .or. ((iSymb == isymx) .and. (iSymy > iSymv))) cycle
                Nax = nAsh(iSymx)
                Nay = nAsh(iSymy)
                if (NAx*Nay /= 0) then

                  i = 2
                  if (.not. Fake_CMO2) i = 3

                  ! (tu|v~w) = Ltu*Lv~w
                  call DGEMM_('N','T',Nav*Naw,NAx*Nay,JNUM,One,Lpq(2)%SB(iSymv)%A3,NAv*Naw,Lpq(i)%SB(iSymx)%A3,NAx*Nay,One, &
                              MOScr%SB(iSymb,iSymv,iSymx)%A2,NAv*Naw)
                end if
              end do
              call CWTIME(TCINT3,TWINT3)
              tint2(1) = tint2(1)+(TCINT3-TCINT2)
              tint2(2) = tint2(2)+(TWINT3-TWINT2)
            end if
          end do

          ! ************ EVALUATION OF THE ACTIVE FOCK MATRIX *************
          ! Coulomb term
          if (JSYM == 1) then

            call mma_allocate(VJ,JNUM,Label='VJ')
            VJ(:) = Zero

            do iSymb=1,nSym

              iSymv = Mul(JSYM,iSymb)
              NAv = nAsh(iSymv)
              NAw = nAsh(iSymb)

              if (NAv*Naw /= 0) then
                i = 3
                if (Fake_CMO2) i = 2

                call DGEMV_('T',Nav*Naw,JNUM,ONE,Lpq(i)%SB(iSymv)%A3,Nav*Naw,DA%SB(iSymB)%A2,1,ONE,VJ,1)
              end if
            end do

            call DGEMV_('N',nRS,JNUM,-FactCI,Lrs,nRS,VJ,1,One,Frs(:,2),1)

            call mma_deallocate(VJ)

          end if
          call CWTIME(TCINT2,TWINT2)
          tact(1) = tact(1)+(TCINT2-TCINT3)
          tact(2) = tact(2)+(TWINT2-TWINT3)
          ! ----------------------------------------------------------------
          ! Formation of the Q matrix Qpx = L~py Lvw Gxyvw
          ! ----------------------------------------------------------------
          do iSymb=1,nSym

            iSymv = Mul(JSYM,iSymb)
            NAv = nAsh(iSymv)
            NAw = nAsh(iSymb)

            if (NAv*Naw /= 0) then

              Lpq(3)%SB(iSymv)%A3(:,:,:) = Zero

              do iSymx=1,nSym
                iSymy = Mul(JSYM,iSymx)
                Nax = nAsh(iSymx)
                Nay = nAsh(iSymy)

                ipG = 1+iASQ(isymb,iSymv,iSymx)

                if (NAx*Nay /= 0) then
                  call DGEMM_('N','N',Nav*Naw,JNUM,NAx*Nay,One,G2(ipG),NAv*Naw,Lpq(2)%SB(iSymy)%A3,NAx*Nay,ONE, &
                              Lpq(3)%SB(iSymv)%A3,Nav*Naw)
                end if

              end do

              do JVC=1,JNUM
                call DGEMM_('T','N',NBAS(iSymb),NAw,Nav,One,Lpq(1)%SB(iSymv)%A3(:,:,JVC),NAv,Lpq(3)%SB(iSymv)%A3(:,:,JVC),Nav,ONE, &
                            QTmp(1)%SB(iSymb)%A2,NBAS(iSymb))
              end do

              !MGD check timing loops are correct (not always timed from previous reference)
            end if
          end do
          call CWTIME(TCINT3,TWINT3)
          tQmat(1) = tQmat(1)+(TCINT3-TCINT2)
          tQmat(2) = tQmat(2)+(TWINT3-TWINT2)

          call Deallocate_DT(Lpq(2))

          call Allocate_DT(Lpq(2),nAsh,nBas,nVec,JSYM,nSym,iSwap)

          ! ************ EVALUATION OF THE ACTIVE FOCK MATRIX *************
          ! Exchange term
          do iSymb=1,nSym

            iSymv = Mul(JSYM,iSymb)
            NAv = nAsh(iSymv)
            NAw = nAsh(iSymb)

            if (NAv*nBas(iSymb) /= 0) then

              do JVC=1,JNUM
                call DGEMM_('T','N',NBAS(iSymb),Nav,Nav,ONE,Lpq(1)%SB(iSymv)%A3(:,:,JVC),Nav,DA%SB(iSymb)%A2,Nav,ZERO, &
                            Lpq(2)%SB(iSymv)%A3(:,:,JVC),NBAS(iSymb))
              end do
              call CWTIME(TCINT2,TWINT2)
              tact(1) = tact(1)+(TCINT2-TCINT3)
              tact(2) = tact(2)+(TWINT2-TWINT3)

              ! ----------------------------------------------------------------
              ! First half Active transformation  L~vb,J = sum_a  ~C(v,a) * Lab,J
              ! ----------------------------------------------------------------
              if (.not. Fake_CMO2) then
                call CWTIME(TCINT2,TWINT2)

                call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,Ash(2),Lpq,DoRead)
                call CWTIME(TCINT3,TWINT3)
                tint1(1) = tint1(1)+(TCINT3-TCINT2)
                tint1(2) = tint1(2)+(TWINT3-TWINT2)
                ! ----------------------------------------------------------------
                ! Formation of the Q matrix Qpx = Lp~y Lvw Gxyvw
                ! ----------------------------------------------------------------
                do JVC=1,JNUM
                  call DGEMM_('T','N',NBAS(iSymb),NAw,Nav,One,Lpq(1)%SB(iSymb)%A3(:,:,JVC),NAv,Lpq(3)%SB(iSymv)%A3(:,:,JVC),Naw, &
                              ONE,QTmp(2)%SB(iSymb)%A2,NBAS(iSymb))
                end do
                call CWTIME(TCINT2,TWINT2)
                tQmat(1) = tQmat(1)+(TCINT2-TCINT3)
                tQmat(2) = tQmat(2)+(TWINT2-TWINT3)
              end if
            end if
          end do ! iSymb

          ! ************ EVALUATION OF THE ACTIVE FOCK MATRIX *************
          ! Exchange term
          !MGD what if Naw=0 but not NBas(b)?
          do iSymb=1,nSym

            iSymv = Mul(JSYM,iSymb)
            NAv = nAsh(iSymv)
            NAw = nAsh(iSymb)

            if (NAv*Naw /= 0) then

              do JVC=1,JNUM
                do is=1,NBAS(iSymb)
                  call DGEMV_('N',NBAS(iSymb),Nav,-FactXI,Lpq(2)%SB(iSymb)%A3(:,:,JVC),NBAS(iSymb),Lpq(1)%SB(iSymv)%A3(:,is,JVC), &
                              1,ONE,KA%SB(iSymB)%A2(:,is),1)

                end do
              end do
              call CWTIME(TCINT3,TWINT3)
              tact(1) = tact(1)+(TCINT3-TCINT2)
              tact(2) = tact(2)+(TWINT3-TWINT2)

            end if

          end do ! iSymb

          call Deallocate_DT(Lpq(3))
          call Deallocate_DT(Lpq(2))
          call Deallocate_DT(Lpq(1))

        end if ! If (DoAct)

        call CWTIME(TCINT2,TWINT2)
        tintg(1) = tintg(1)+(TCINT2-TCINT1)
        tintg(2) = tintg(2)+(TWINT2-TWINT1)

        if (irc /= 0) return

        ! ---------------- END (TW|XY) EVALUATION -----------------------

      end do ! end batch loop

      if (JSYM == 1) then
        ! backtransform fock matrix to full storage
        add = .true.
        nMat = 1
        call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,JI,Frs(:,1),add)
        if (DoAct) then
          call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,JALT,Frs(:,2),add)
        end if
      end if

      ! free memory
      call mma_deallocate(Lrs)

      if (JSYM == 1) then
        call mma_deallocate(Frs)
        call mma_deallocate(Drs)
      end if

    end if

    ! Screening control section
    DoScreen = kscreen == Nscreen

    if (.not. DoScreen) then
      kscreen = kscreen+1
    else
      kscreen = 1
    end if

#   ifdef _MOLCAS_MPP_
    if ((nProcs > 1) .and. Update .and. DoScreen .and. Is_Real_Par()) then
      call GaDsum(DiagJ,nnBSTR(JSYM,1))
      call Daxpy_(nnBSTR(JSYM,1),-One,DiagJ,1,Diag(1+iiBstR(JSYM,1)),1)
      call Fzero(DiagJ,nnBSTR(JSYM,1))
    end if
    ! Need to activate the screening to setup the contributing shell
    ! indices the first time the loop is entered .OR. whenever other nodes
    ! have performed screening in the meanwhile
    if ((nProcs > 1) .and. (.not. DoScreen) .and. (nVrs == 0) .and. Is_Real_Par()) then
      ntv0 = ntv0+1
      DoScreen = ((JRED < myJRED1) .or. (ntv0 >= Nscreen))
      if (DoScreen) ntv0 = 0
    end if
#   endif

  end do ! loop over red sets

end do ! loop over JSYM

! Accumulate Coulomb and Exchange contributions
do iSym=1,nSym

  do iaSh=1,nShell

    ioffa = kOffSh(iaSh,iSym)

    do ibSh=1,nShell

      iShp = nShell*(iaSh-1)+ibSh
      iOffAB = nnBfShp(iShp,iSym)

      iShp = nShell*(ibSh-1)+iaSh
      jOffAB = nnBfShp(iShp,iSym)

      ioffb = kOffSh(ibSh,iSym)

      do ib=1,nBasSh(iSym,ibSh)

        do ia=1,nBasSh(iSym,iaSh)

          iab = nBasSh(iSym,iaSh)*(ib-1)+ia
          !MGD warning with sym
          jab = nBasSh(iSym,ibSh)*(ia-1)+ib

          iag = ioffa+ia
          ibg = ioffb+ib

          iabg = iTri(iag,ibg)

          FkI%SB(iSym)%A2(iag,ibg) = JI(1)%SB(iSym)%A1(iabg)+KI%SB(iSym)%A1(iOffAB+iab)+KI%SB(iSym)%A1(jOffAB+jab)

          FkA%SB(iSym)%A2(iag,ibg) = JALT(1)%SB(iSym)%A1(iabg)+KA%SB(iSym)%A2(iag,ibg)+KA%SB(iSym)%A2(ibg,iag)

        end do

      end do

    end do

  end do

end do
call CWTIME(TCINT1,TWINT1)

call Deallocate_DT(JALT(1))

if (DoAct) then

  !* Compute MO integrals

  do iS=1,nSym
    do jS=1,nsym
      ijS = ieor(is-1,js-1)+1
      do kS=1,nSym
        ls = ieor(ijs-1,ks-1)+1

        do iAsh=1,nAsh(is)
          do jAsh=1,nAsh(js)
            iij = itri(iAsh+kAOff(is),jAsh+kAOff(jS))
            Fac1 = One
            if (iAsh+kAOff(is) == jAsh+kAoff(jS)) Fac1 = Two

            do kAsh=1,nAsh(ks)
              do lAsh=1,nAsh(ls)
                ikl = itri(lAsh+kAOff(lS),kAsh+kAOff(kS))
                Fac2 = One
                if (lAsh+kAOff(lS) == kAsh+kAOff(kS)) Fac2 = Two
                if (iij /= ikl) Fac2 = Fac2*Half

                ipG = itri(iij,ikl)

                MO_Int(ipG) = MO_Int(ipG)+Fac1*Fac2*MOScr%SB(iS,jS,kS)%A4(iAsh,jAsh,kAsh,lAsh)
              end do
            end do
          end do
        end do

      end do
    end do
  end do
end if

!* Transform Fock and Q matrix to MO basis

do iS=1,nSym
  jS = iS
  if (nBas(iS) /= 0) then
    if (DoAct) then
      call DGEMM_('T','N',nBas(jS),nBas(iS),nBas(iS),One,FkA%SB(iS)%A2,nBas(iS),CMO%SB(iS)%A2,nBas(iS),Zero,JA%SB(iS)%A2,nBas(jS))
      call DGEMM_('T','N',nBas(jS),nBas(jS),nBas(iS),One,JA%SB(iS)%A2,nBas(iS),CMO%SB(jS)%A2,nBas(jS),Zero,FkA%SB(iS)%A2,nBas(jS))
    end if
    call DGEMM_('T','N',nBas(jS),nBas(iS),nBas(iS),One,FkI%SB(iS)%A2,nBas(iS),CMO%SB(iS)%A2,nBas(iS),Zero,JA%SB(iS)%A2,nBas(jS))
    call DGEMM_('T','N',nBas(jS),nBas(jS),nBas(iS),One,JA%SB(iS)%A2,nBas(iS),CMO%SB(jS)%A2,nBas(jS),Zero,FkI%SB(iS)%A2,nBas(jS))
    if (DoAct) then
      if (Fake_CMO2) then
        call DGEMM_('T','N',nBas(jS),nAsh(iS),nBas(jS),One,CMO%SB(iS)%A2,nBas(jS),QTmp(1)%SB(js)%A2,nBas(jS),Zero,QVec%SB(js)%A2, &
                    nBas(jS))
      else
        call DGEMM_('T','N',nBas(jS),nAsh(iS),nBas(jS),One,CMO%SB(iS)%A2,nBas(jS),QTmp(2)%SB(jS)%A2,nBas(jS),Zero,QVec%SB(jS)%A2, &
                    nBas(jS))
        call DGEMM_('T','N',nBas(jS),nAsh(iS),nBas(jS),One,CMO%SB(iS)%A2,nBas(jS),QTmp(1)%SB(jS)%A2,nBas(jS),Zero, &
                    QTmp(2)%SB(jS)%A2,nBas(jS))
        call DGEMM_('N','N',nBas(jS),nAsh(iS),nBas(jS),-One,Kappa%SB(iS)%A2,nBas(jS),QTmp(2)%SB(jS)%A2,nBas(jS),One, &
                    QVec%SB(jS)%A2,nBas(jS))
      end if
    end if
  end if
end do

call CWTIME(TCINT2,TWINT2)
tint3(1) = tint3(1)+(TCINT2-TCINT1)
tint3(2) = tint3(2)+(TWINT2-TWINT1)

if (DoAct) then
  call Deallocate_DT(QTmp(2))
  call Deallocate_DT(QTmp(1))
  call Deallocate_DT(MOScr)
end if

call mma_deallocate(Faa)
call mma_deallocate(Fia)
call mma_deallocate(SvShp)
call mma_deallocate(iShp_rs)
call mma_deallocate(nnBfShp)
call mma_deallocate(kOffSh)
call mma_deallocate(Indx)
call mma_deallocate(SumAClk)
call mma_deallocate(MLk)
call mma_deallocate(Ylk)
call mma_deallocate(AbsC)
call Deallocate_DT(DiaH)
#ifdef _MOLCAS_MPP_
if ((nProcs > 1) .and. Update .and. Is_Real_Par()) call mma_deallocate(DiagJ)
#endif
call mma_deallocate(Diag)

if (Deco) then
  call Deallocate_DT(CM(2))
  call Deallocate_DT(CM(1))
end if

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky MCLR timing from '//SECNAM
  write(u6,CFmt) '----------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'Fock matrix construction        CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ VECTORS                              ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'COULOMB                                   ',tcoul(1),tcoul(2)
  write(u6,'(2x,A26,2f10.2)') 'EXCHANGE                                  ',tscrn(1)+tmotr(1)+texch(1),tscrn(2)+tmotr(2)+texch(2)
  write(u6,'(2x,A26,2f10.2)') '  SCREENING                               ',tscrn(1),tscrn(2)
  write(u6,'(2x,A26,2f10.2)') '  MO TRANSFORM                            ',tmotr(1),tmotr(2)
  write(u6,'(2x,A26,2f10.2)') '  FORMATION                               ',texch(1),texch(2)
  write(u6,'(2x,A26,2f10.2)') 'ACTIVE INT, Q AND FOCK MATRIX             ',tint1(1)+tint2(1)+tint3(1)+tQmat(1)+tact(1), &
                              tint1(2)+tint2(2)+tint3(2)+tQmat(2)+tact(2)
  write(u6,'(2x,A26,2f10.2)') '  MO TRANSFORM                            ',tint1(1),tint1(2)
  write(u6,'(2x,A26,2f10.2)') '  INTEGRAL                                ',tint2(1),tint2(2)
  write(u6,'(2x,A26,2f10.2)') '  Q MATRIX                                ',tQmat(1),tQmat(2)
  write(u6,'(2x,A26,2f10.2)') '  ACTIVE FOCK MATRIX                      ',tact(1),tact(2)
  write(u6,'(2x,A26,2f10.2)') '  MO BACK TRANSFORM                       ',tint3(1),tint3(2)
  write(u6,*)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
if (Debug) then ! to avoid double printing in RASSI-debug

  write(u6,'(6X,A)') 'TEST PRINT FROM '//SECNAM
  write(u6,'(6X,A)')
  write(u6,'(6X,A)') '***** INACTIVE FOCK MATRIX ***** '
  do ISYM=1,NSYM
    if (NBAS(ISYM) > 0) then
      write(u6,'(6X,A)')
      write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
      call CHO_OUTPUT(FkI%SB(ISYM)%A2,1,NBAS(ISYM),1,NBAS(ISYM),NBAS(ISYM),NBAS(ISYM),1,u6)
    end if
  end do

end if
#endif

return

end subroutine CHO_LK_MCLR
