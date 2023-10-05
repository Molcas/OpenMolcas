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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine CHO_LK_RASSI_X(DLT,MSQ,FLT,KSQ,FSQ,TUVX,Ash,nScreen,dmpk)
!*********************************************************************
!  Author : F. Aquilante
!
!  Note:  this routine differs from CHO_LK_RASSI because it can
!         handle ALSO the case where the 2 sets of MOs are different!
!         The Exchange contribution is non-symmetric and so is FI
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

use Cholesky, only: iiBstR, IndRed, InfVec, MaxRed, nBas, nBasSh, nDimRS, nnBstR, nnBstRsh, nnBstRT, nnShl, nnShl_tot, nShell, &
                    nSym, NumCho, NumChT, timings
use Fock_util_interface, only: cho_lr_MOs
use Fock_util_global, only: Deco, Estimate, Fake_CMO2, PseudoChoMOs, Update
use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Data_Structures, only: DSBA_Type, NDSBA_Type, SBA_Type, twxy_Type
use Cholesky_Structures, only: Allocate_DT, Deallocate_DT, L_Full_Type, Lab_Type
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
type(DSBA_Type), intent(in) :: DLT, Ash(2)
type(DSBA_Type), intent(inout) :: MSQ(2), FLT(1), KSQ
type(DSBA_Type), intent(_OUT_) :: FSQ
real(kind=wp), intent(_OUT_) :: TUVX(*)
integer(kind=iwp), intent(in) :: nScreen
real(kind=wp), intent(in) :: dmpk
#include "warnings.h"
#include "rassi.fh"
integer(kind=iwp) :: ia, iab, iabg, iag, iaSh, iaSkip, ib, iBatch, ibcount, ibg, ibs, ibSh, ibSkip, iCase, iE, ik, iLoc, iml, Inc, &
                     ioffa, iOffAB, ioffb, iOffShb, iOK, irc, ired1, IREDC, iS, ish, iShp, iSwap, ISYM, iSyma, iSymb, iSymv, iTmp, &
                     IVEC2, iVrs, jaSkip, jden, jK, jK_a, jml, jmlmax, JNUM, JRED, JRED1, JRED2, jrs, jSym, jvc, JVEC, k, kDen, &
                     kMOs, kOff(8), krs, kscreen, kSym, l, l1, LFMAX, LFULL(2), LKsh, LKshp, LREAD, lSh, lSym, LWORK, MaxB, &
                     MaxRedT, MaxVecPerBatch, Mmax, mrs, mSh, mTVec, mTvec1, mTvec2, MUSED, MxB, MxBasSh, myJRED2, n1, n2, &
                     nAux(8), NAv, NAw, nBatch, nBsa, nDen, nMat, nMOs, nnO, nnShl_2, nRS, NumCV, numSh1, numSh2, NUMV, NumVT, &
                     nVec, nVrs
real(kind=wp) :: Fact, LKThr, SKsh, tau, TCC1, TCC2, TCINT1, TCINT2, tcoul(2), TCR1, TCR2, TCS1, TCS2, TCT1, TCT2, TCX1, TCX2, &
                 texch(2), THRSX, thrv, tintg(2), tmotr(2), Tmp, TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), &
                 tscrn(2), TWC1, TWC2, TWINT1, TWINT2, TWR1, TWR2, TWS1, TWS2, TWT1, TWT2, TWX1, TWX2, xFab, xtau, xTmp, YMax, &
                 YshMax
logical(kind=iwp) :: add, DoScreen, DoReord
character(len=50) :: CFmt
character :: mode, mode2
type(NDSBA_Type) :: DiaH
type(DSBA_Type) :: CM(2), MO(2)
type(SBA_Type) :: Laq(2)
type(twxy_Type) :: Scr
type(L_Full_Type) :: L_Full
type(Lab_Type) :: Lab
integer(kind=iwp), allocatable :: Indx(:,:,:), iShp_rs(:), kOffSh(:,:), nnBfShp(:,:)
real(kind=wp), allocatable :: AbsC(:), Diag(:), Drs(:), Faa(:), Fia(:), Frs(:), Lrs(:,:), MLk(:,:,:), SumAClk(:,:,:), SvShp(:,:), &
                              VJ(:), Ylk(:,:,:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: i, myJRED1, NNBSTMX, ntv0
real(kind=wp), allocatable :: DiagJ(:)
#endif
character(LEN=14), parameter :: SECNAM = 'CHO_LK_RASSI_X'
logical(kind=iwp), parameter :: DoRead = .false.
real(kind=wp), parameter :: FactCI = One, FactXI = -One
integer(kind=iwp), external :: Cho_LK_MaxVecPerBatch
real(kind=wp), external :: Cho_LK_ScreeningThreshold, ddot_

!                                                                      *
!***********************************************************************
!                                                                      *
DoReord = .false.
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

! ==================================================================

! Various offsets
! ----------------
nnO = 0
kOff(1) = 0
MaxB = nBas(1)
do ISYM=2,NSYM
  MaxB = max(MaxB,nBas(iSym))
  nnO = nnO+nIsh(iSym-1)
  kOff(iSym) = nnO
end do
nnO = nnO+nIsh(nSym)
nAux(:) = nIsh(:)+nAsh(:)

!*************************************************
if (Deco) then

  call Allocate_DT(CM(1),nBas,nAux,nSym)
  call Allocate_DT(CM(2),nBas,nAux,nSym)

  if (PseudoChoMOs) then
    call cho_get_MO(iOK,nDen,nSym,nBas,nIsh,MSQ,CM)
  else
    call cho_lr_MOs(iOK,nDen,nSym,nBas,nIsh,MSQ,CM)
  end if

  if (iOK == 0) then ! point to the "generalized" Cholesky MOs
    do jden=1,nDen
      call Allocate_DT(MO(jDen),nBas,nAux,nSym,Ref=CM(jDen)%A0)
    end do
    !write(u6,*) 'Cholesky MOs used for state A'
    !if (nDen == 2) write(u6,*) 'Pseudo Cholesky MOs used for state B'
  else
    write(u6,*) '*******************************'
    write(u6,*) '*** Resort to Canonical MOs ***'
    write(u6,*) '*******************************'

    do jden=1,nDen
      call Allocate_DT(MO(jDen),nBas,nAux,nSym,Ref=MSQ(jDen)%A0)
    end do

  end if

else

  do jden=1,nDen
    call Allocate_DT(MO(jDen),nBas,nAux,nSym,Ref=MSQ(jDen)%A0)
  end do

end if
!*************************************************

! Define the max number of vectors to be treated in core at once

MaxVecPerBatch = Cho_LK_MaxVecPerBatch()

! Define the screening threshold

! threshold for max BLB matrix element
! Note: must be consistent with threshold in subroutine rasscf/rasscf_init.f
THRSX = 1.0e-4_wp
LKThr = Cho_LK_ScreeningThreshold(THRSX)
tau = (LKThr/max(1,nnO))*dmpk

MaxRedT = MaxRed
call GAIGOP_SCAL(MaxRedT,'+')

if (Estimate) tau = tau/MaxRedT

xtau = sqrt(tau)

! Vector MO transformation screening thresholds
NumVT = NumChT
call GAIGOP_SCAL(NumVT,'+')
thrv = (sqrt(LKThr/(max(1,nnO)*NumVT)))*dmpk

call mma_allocate(DIAG,NNBSTRT(1),Label='DIAG')

#ifdef _MOLCAS_MPP_
if ((nProcs > 1) .and. Update .and. Is_Real_Par()) then
  NNBSTMX = 0
  do i=1,nSym
    NNBSTMX = max(NNBSTMX,NNBSTR(i,1))
  end do
  call mma_allocate(diagJ,NNBSTMX,Label='DiagJ')
  DiagJ(:) = Zero
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
call mma_allocate(Ylk,MaxB,nno,nDen,Label='Ylk')

! allocate memory for the ML[k] lists of largest elements in significant shells
call mma_allocate(MLk,nShell,nnO,nDen,Label='MLk')

! allocate memory for the lists of  S:= sum_l abs(C(l)[k]) for each shell
call mma_allocate(SumAClk,nShell,nnO,nDen,Label='SumAClk')

! allocate memory for the Index arrays
call mma_allocate(Indx,[0,nShell],[1,nnO],[1,nDen],Label='Indx')

! allocate memory for kOffSh
call mma_allocate(kOffSh,nShell,nSym,Label='kOffSh')

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

    do jK=1,nIsh(kSym)
      jK_a = jK+kOff(kSym)

      do iaSh=1,nShell

        SKsh = zero
        iS = kOffSh(iaSh,kSym)+1
        iE = kOffSh(iaSh,kSym)+nBasSh(kSym,iaSh)
        do ik=iS,iE
          SKsh = SKsh+MO(jDen)%SB(kSym)%A2(ik,jK)**2
        end do

        SumAClk(iaSh,jK_a,jDen) = SKsh

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

! *************** BIG LOOP OVER VECTORS SYMMETRY *******************
do jSym=1,nSym

  NumCV = NumCho(jSym)
  call GAIGOP_SCAL(NumCV,'max')
  if (NumCV < 1) cycle

  JNUM = 1
  call Allocate_DT(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym,Memory=LFULL)

  iCase = 0
  call Allocate_DT(Scr,nAsh,nAsh,JSYM,nSym,iCase)

  iLoc = 3 ! use scratch location in reduced index arrays

  ! ****************     MEMORY MANAGEMENT SECTION    *****************
  ! --------------------------------------------------------------
  ! compute memory needed to store at least 1 vector of JSYM
  ! and do all the subsequent calculations
  ! --------------------------------------------------------------
  mTvec1 = 0
  mTvec2 = 0
  MxB = 0
  do l=1,nSym
    k = Mul(l,JSYM)
    Mmax = max(0,nIsh(k))
    if (Mmax > 0) MxB = max(MxB,nBas(l))
    mTvec1 = mTvec1+nAsh(k)*nBas(l)
    mTvec2 = mTvec2+nAsh(k)*nAsh(l)
  end do
  mTVec = mTVec1+mTVec2

  LFMAX = max(mTvec,LFULL(1)) ! re-use memory for the active vec
  mTvec = nDen*max(MxB,1) ! mem for storing half-transformed vec

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

  JRED1 = InfVec(1,2,jSym)            ! red set of the 1st vec
  JRED2 = InfVec(NumCho(jSym),2,jSym) ! red set of the last vec
# ifdef _MOLCAS_MPP_
  myJRED1 = JRED1 ! first red set present on this node
  ntv0 = 0
# endif
  myJRED2 = JRED2 ! last  red set present on this node

  ! entire red sets range for parallel run
  call GAIGOP_SCAL(JRED1,'min')
  call GAIGOP_SCAL(JRED2,'max')

  kscreen = 1
  DoScreen = .true.

  do JRED=JRED1,JRED2

    call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

    if (nVrs /= 0) then ! otherwise no vectors in that (jred,jsym)

      if (nVrs < 0) then
        write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
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
        call mma_allocate(Frs,nRS,Label='Frs')
        Drs(:) = Zero
        Frs(:) = Zero
      end if

      call mma_maxDBLE(LWORK)

      nVec = min((LWORK-LFULL(2))/(nRS+mTvec+LFMAX),min(nVrs,MaxVecPerBatch))

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

        call mma_allocate(Lrs,nRS,JNUM,Label='Lrs')
        Lrs(:,:) = Zero

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

          call DGEMV_('N',nRS,JNUM,FactCI,Lrs,nRS,VJ,1,Fact,Frs,1)

          call mma_deallocate(VJ)

          call CWTIME(TCC2,TWC2)
          tcoul(1) = tcoul(1)+(TCC2-TCC1)
          tcoul(2) = tcoul(2)+(TWC2-TWC1)

        end if ! Coulomb contribution

        ! *************** EXCHANGE CONTRIBUTIONS  ***********************

        call CWTIME(TCS1,TWS1)
        ! -----------------------------------------------------------------
        ! Estimate the diagonals :   D(a,b) = sum_J (Lab,J)^2

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
        !                                                              *
        !***************************************************************
        !***************************************************************
        !***************************************************************
        !                                                              *
        call Allocate_DT(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym)
        call Allocate_DT(Lab,JNUM,nBasSh,nBas,nShell,nSym,nDen)

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
          ! ----------------------------------------------------------------
          ired1 = 1 ! location of the 1st red set
          call swap_tosqrt(irc,ired1,NNBSTRT(1),JSYM,DIAH,DIAG)

          call CWTIME(TCS2,TWS2)
          tscrn(1) = tscrn(1)+(TCS2-TCS1)
          tscrn(2) = tscrn(2)+(TWS2-TWS1)

        end if

        do kSym=1,nSym

          lSym = Mul(JSYM,kSym)

          do jK=1,nIsh(kSym)

            jK_a = jK+kOff(kSym)

            Lab%A0(1:nDen*nBas(lSym)*JNUM) = Zero

            if (DoScreen) then

              call CWTIME(TCS1,TWS1)
              !--------------------------------------------------------------
              ! Setup the screening
              !--------------------------------------------------------------

              do jDen=1,nDen

                do ik=1,nBas(kSym)
                  Absc(ik) = abs(MO(jDen)%SB(kSym)%A2(ik,jK))
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
                  Indx(iSh,jK_a,jDen) = ish
                end do
              end do

              ! ****  The Max in the MO set 1 is used as reference
              numSh1 = 0  ! # of significant shells in MO set 1
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
                iTmp = Indx(1,jK_a,1)
                MLk(1,jK_a,1) = YMax
                Indx(1,jK_a,1) = Indx(jmlmax,jK_a,1)
                MLk(jmlmax,jK_a,1) = xTmp
                Indx(jmlmax,jK_a,1) = iTmp
              end if

              ! **** Sort the list for the MO set 2   iff  MOs1 /= MOs2
              if (.not. Fake_CMO2) then
                numSh2 = 0  ! # of significant shells in MO set 2
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
                    iTmp = Indx(jml,jK_a,2)
                    MLk(jml,jK_a,2) = YMax
                    Indx(jml,jK_a,2) = Indx(jmlmax,jK_a,2)
                    MLk(jmlmax,jK_a,2) = xTmp
                    Indx(jmlmax,jK_a,2) = iTmp
                  end if

                  ! Exact bounds (quadratic scaling of the MO transformation)
                  ! Note that in true RASSI the exchange matrix is not
                  ! positive definite.

                  if (MLk(jml,jK_a,2)*MLk(1,jK_a,1) >= tau) then
                    numSh2 = numSh2+1
                  else
                    jml = nShell ! exit the loop
                  end if

                  jml = jml+1

                end do

                Indx(0,jK_a,2) = numSh2
                numSh1 = 1

              else ! fake biorthonormal basis

                numSh2 = 6669666 ! dummy assignement

                if (MLk(1,jK_a,1) >= xtau) numSh1 = 1

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
                    iTmp = Indx(jml,jK_a,1)
                    MLk(jml,jK_a,1) = YMax
                    Indx(jml,jK_a,1) = Indx(jmlmax,jK_a,1)
                    MLk(jmlmax,jK_a,1) = xTmp
                    Indx(jmlmax,jK_a,1) = iTmp
                  end if

                  if ((.not. Fake_CMO2) .and. (MLk(jml,jK_a,1)*MLk(1,jK_a,kDen) >= tau)) then
                    numSh1 = numSh1+1

                  ! Here we use a non-exact bound for the exchange matrix because a
                  ! fake rassi (MOs1=MOs2) has a positive definite exchange
                  else if (Fake_CMO2 .and. (MLk(jml,jK_a,1) >= xtau)) then
                    numSh1 = numSh1+1
                  else
                    jml = nShell ! exit the loop
                  end if

                  jml = jml+1

                end do
              else
                numSh1 = 0
              end if

              Indx(0,jK_a,1) = numSh1

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

              do iSh=1,Indx(0,jK_a,jDen)

                iaSh = Indx(ish,jK_a,jDen)

                Lab%Keep(iaSh,jDen) = .true.

                ibcount = 0

                do ibSh=1,nShell

                  iOffShb = kOffSh(ibSh,kSym)

                  iShp = iTri(iaSh,ibSh)

                  if (iShp_rs(iShp) <= 0) cycle

                  if ((nnBstRSh(JSym,iShp_rs(iShp),iLoc)*nBasSh(lSym,iaSh)*nBasSh(kSym,ibSh) > 0) .and. &
                      (sqrt(abs(SumAClk(ibSh,jk_a,jDen)*SvShp(iShp_rs(iShp),1))) >= thrv)) then

                    ibcount = ibcount+1

                    if (lSym >= kSym) then

                      l1 = 1
                      if (iaSh < ibSh) l1 = 2

                      ! LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
                      ! ----------------------------------
                      Mode = 'N'
                      n1 = nBasSh(lSym,iaSh)*JNUM
                      n2 = nBasSh(kSym,ibSh)

                      call DGEMV_(Mode,n1,n2,One,L_Full%SPB(lSym,iShp_rs(iShp),l1)%A21,n1,MO(jDen)%SB(kSym)%A2(iOffShb+1:,jK),1, &
                                  ONE,Lab%SB(iaSh,lSym,jDen)%A,1)

                    else ! lSym < kSym

                      l1 = 1
                      if (ibSh < iaSh) l1 = 2

                      ! LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
                      ! ----------------------------------
                      Mode = 'T'
                      n1 = nBasSh(kSym,ibSh)
                      n2 = JNUM*nBasSh(lSym,iaSh)

                      call DGEMV_(Mode,n1,n2,One,L_Full%SPB(kSym,iShp_rs(iShp),l1)%A12,n1,MO(jDen)%SB(kSym)%A2(iOffShb+1:,jK),1, &
                                  ONE,Lab%SB(iaSh,lSym,jDen)%A,1)

                    end if

                  end if

                end do ! ibsh

                ! The following re-assignement is used later on to check if the
                ! iaSh vector LaJ[k] can be neglected because identically zero

                if (ibcount == 0) Lab%Keep(iash,jDen) = .false.

              end do ! iSh

            end do ! jDen

            call CWTIME(TCT2,TWT2)
            tmotr(1) = tmotr(1)+(TCT2-TCT1)
            tmotr(2) = tmotr(2)+(TWT2-TWT1)

            ! Prepare the J-screening

            call CWTIME(TCS1,TWS1)

            do iSh=1,Indx(0,jK_a,1)

              iaSh = Indx(iSh,jK_a,1)

              iaSkip = merge(1,0,Lab%Keep(iash,1))

              jaSkip = merge(1,0,Lab%Keep(iash,kDen))

              if (iaSkip*jaSkip == 0) cycle

              if (lSym >= kSym) then

                ! Faa,[k] = sum_J  LaJ[k1]*LaJ[k2]
                ! --------------------------------
                Inc = nBasSh(lSym,iaSh)
                n1 = 1

              else ! lSym < kSym

                ! Faa,[k] = sum_J  LJa[k1]*LJa[k2]
                ! --------------------------------
                Inc = 1
                n1 = JNUM

              end if

              Tmp = Zero
              do ia=1,nBasSh(lSym,iaSh)
                Fia(ia) = DDot_(JNUM,Lab%SB(iaSh,lSym,1)%A(1+n1*(ia-1):),Inc,Lab%SB(iaSh,lSym,kDen)%A(1+n1*(ia-1):),Inc)
                Tmp = max(abs(Fia(ia)),Tmp)
              end do

              Faa(iaSh) = Tmp

            end do

            call CWTIME(TCS2,TWS2)
            tscrn(1) = tscrn(1)+(TCS2-TCS1)
            tscrn(2) = tscrn(2)+(TWS2-TWS1)

            !--------------------------------------------------------
            ! Compute exchange matrix for the interacting shell pairs
            !--------------------------------------------------------

            call CWTIME(TCX1,TWX1)

            do lSh=1,Indx(0,jK_a,1)

              iaSh = Indx(lSh,jK_a,1)

              iaSkip = merge(1,0,Lab%Keep(iash,kDen))

              mSh = 1

              do while (mSh <= Indx(0,jK_a,kDen))

                ibSh = Indx(mSh,jK_a,kDen)

                ibSkip = merge(1,0,Lab%Keep(ibsh,1))

                iShp = nShell*(iaSh-1)+ibSh

                iOffShb = kOffSh(ibSh,lSym)

                iOffAB = nnBfShp(iShp,lSym)

                xFab = sqrt(abs(Faa(iaSh)*Faa(ibSh)))

                if (MLk(lSh,jK_a,1)*MLk(mSh,jK_a,kDen) < tau) then

                  mSh = Indx(0,jK_a,kDen) ! skip rest

                else if ((xFab >= tau/MaxRedT) .and. (iaSkip*ibSkip == 1)) then

                  nBsa = nBasSh(lSym,iaSh)
                  if (lSym >= kSym) then

                    ! F(a,b)[k] = F(a,b)[k] + FactXI * sum_J  X2(a,J)[k] * X1(b,J)[k]
                    ! ---------------------------------------------------------------

                    n1 = nBasSh(lSym,iaSh)
                    n2 = nBasSh(lSym,ibSh)
                    Mode = 'N'
                    Mode2 = 'T'

                  else ! lSym < kSym

                    ! F(a,b)[k] = F(a,b)[k] + FactXI * sum_J  X2(J,a)[k] * X1(J,b)[k]
                    ! ---------------------------------------------------------------

                    n1 = JNUM
                    n2 = JNUM
                    Mode = 'T'
                    Mode2 = 'N'

                  end if

                  call DGEMM_(Mode,Mode2,nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,FactXI,Lab%SB(iaSh,lSym,kDen)%A,n1, &
                              Lab%SB(ibSh,lSym,1)%A,n2,ONE,KSQ%SB(lSym)%A1(iOffAB+1:),nBsa)

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

        ! The memory used before for the full-dimension AO-vectors
        !     is now re-used to store half and full transformed
        !     vectors in the active space
        ! ---------------------------------------------------------
        iSwap = 0 ! Lvb,J are returned
        call Allocate_DT(Laq(1),nAsh,nBas,JNUM,JSYM,nSym,iSwap)
        call Allocate_DT(Laq(2),nAsh,nAsh,JNUM,JSYM,nSym,iSwap)

        kMOs = 1
        nMOs = 1 ! Active MOs (1st set)

        call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,Ash,Laq,DoRead)

        if (irc /= 0) return

        ! ----------------------------------------------------------------
        ! Active-Active transformation  Lvw,J = sum_b  Lvb,J * C2(w,b)
        ! ----------------------------------------------------------------
        do iSymb=1,nSym

          iSymv = Mul(JSYM,iSymb)
          NAv = nAsh(iSymv)
          NAw = nAsh(iSymb) ! iSymb=iSymw

          if (NAv*NAw /= 0) then

            do JVC=1,JNUM

              call DGEMM_('N','T',NAv,NAw,NBAS(iSymb),One,Laq(1)%SB(iSymv)%A3(:,:,JVC),NAv,Ash(kDen)%SB(iSymb)%A2,NAw,Zero, &
                          Laq(2)%SB(iSymv)%A3(:,:,JVC),NAv)

            end do

          end if

        end do

        ! *************** EVALUATION OF THE (TW|XY) INTEGRALS ***********

        DoReord = (JRED == myJRED2) .and. (iBatch == nBatch)

        call CHO_rassi_twxy(irc,Scr,Laq(2),TUVX,nAsh,JSYM,JNUM,DoReord)

        call CWTIME(TCINT2,TWINT2)
        tintg(1) = tintg(1)+(TCINT2-TCINT1)
        tintg(2) = tintg(2)+(TWINT2-TWINT1)

        if (irc /= 0) return

        ! ---------------- END (TW|XY) EVALUATION -----------------------

        call Deallocate_DT(Laq(2))
        call Deallocate_DT(Laq(1))
        call mma_deallocate(Lrs)

      end do ! end batch loop

      if (JSYM == 1) then
        ! backtransform fock matrix to full storage
        add = .true.
        nMat = 1
        call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,FLT,Frs,add)
      end if

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
      DoScreen = (JRED < myJRED1) .or. (ntv0 >= Nscreen)
      if (DoScreen) ntv0 = 0
    end if
#   endif

  end do ! loop over red sets

  call Deallocate_DT(Scr)

end do ! loop over JSYM

! Accumulate Coulomb and Exchange contributions
do iSym=1,nSym

  do iaSh=1,nShell

    ioffa = kOffSh(iaSh,iSym)

    do ibSh=1,nShell

      iShp = nShell*(iaSh-1)+ibSh

      iOffAB = nnBfShp(iShp,iSym)

      ioffb = kOffSh(ibSh,iSym)

      do ib=1,nBasSh(iSym,ibSh)

        do ia=1,nBasSh(iSym,iaSh)

          iab = nBasSh(iSym,iaSh)*(ib-1)+ia

          iag = ioffa+ia
          ibg = ioffb+ib

          iabg = iTri(iag,ibg)

          FSQ%SB(iSym)%A2(iag,ibg) = FLT(1)%SB(iSym)%A1(iabg)+KSQ%SB(iSym)%A1(iOffAB+iab)

        end do

      end do

    end do

  end do

end do

call mma_deallocate(Fia)
call mma_deallocate(Faa)
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
do jden=1,nDen
  call Deallocate_DT(MO(jDen))
end do

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky RASSI timing from '//SECNAM
  write(u6,CFmt) '----------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'Fock matrix construction        CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ VECTORS                              ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'COULOMB                                   ',tcoul(1),tcoul(2)
  write(u6,'(2x,A26,2f10.2)') 'EXCHANGE                                  ',texch(1),texch(2)
  write(u6,'(2x,A26,2f10.2)') '(TW|XY) INTEGRALS                         ',tintg(1),tintg(2)
  write(u6,*)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
!if (Debug) then ! to avoid double printing in RASSI-debug
write(u6,'(6X,A)') 'TEST PRINT FROM '//SECNAM
write(u6,'(6X,A)')
write(u6,'(6X,A)') '***** INACTIVE FOCK MATRIX ***** '
do ISYM=1,NSYM
  if (NBAS(ISYM) > 0) then
    write(u6,'(6X,A)')
    write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
    !call CHO_OUTPUT(FSQ%SB(ISYM)%A2,1,NBAS(ISYM),1,NBAS(ISYM),NBAS(ISYM),NBAS(ISYM),1,u6)
    call Chk4Nan(nBas(iSym)**2,FSQ%SB(ISYM)%A2,iErr)
    if (iErr /= 0) then
      write(u6,*) 'CHO_LK_RASSI_X FSQ corrupted!'
      call Abend()
    end if
  end if
end do
! endif
#endif

return

end subroutine CHO_LK_RASSI_X
