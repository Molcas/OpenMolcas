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

subroutine CHO_LK_SCF(rc,nDen,FLT,KLT,nForb,nIorb,Porb,PLT,FactXI,nScreen,dmpk,dFmat)
!*********************************************************************
!  Author : F. Aquilante
!
! *************** INACTIVE AO-BASIS FOCK MATRIX **********************
!
!   F(ab) = sum_J  Lab,J * V(J)  -  sum_Jk  Lka,J * Lkb,J
!
!   Exchange term computed using the Local-exchange (LK) algorithm
!
!*********************************************************************
!
!      V(J) = sum_gd  Lgd,J * Ptot(gd)
!
!      a,b,g,d:  AO-index
!      k:        MO-index   belonging to (Frozen+Inactive)
!
!*********************************************************************

use ChoArr, only: nBasSh, nDimRS
use ChoSwp, only: IndRed, InfVec, nnBstRSh
use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Fock_util_global, only: Estimate, Update
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, L_Full_Type, Lab_Type, NDSBA_Type
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: nDen, nForb(8,nDen), nIorb(8,nDen), nScreen
type(DSBA_Type), intent(inout) :: FLT(nDen), KLT(nDen)
type(DSBA_Type), intent(in) :: Porb(nDen), PLT(nDen)
real(kind=wp), intent(in) :: FactXI, dmpk, dFmat
#include "chotime.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "warnings.h"
integer(kind=iwp) :: i, i1, ia, iab, iabg, iag, iaSh, iaSkip, ib, iBatch, ibcount, ibg, ibs, ibSh, ibSkip, iE, ik, iLoc, iml, Inc, &
                     ioffa, iOffAB, ioffb, iOffShb, irc, ired1, IREDC, iS, ish, iShp, ISYM, iSyma, iTmp, IVEC2, iVrs, jDen, jK, &
                     jK_a, jml, jmlmax, JNUM, JRED, JRED1, JRED2, jrs, jSym, jvc, JVEC, k, kOff(8,2), krs, kscreen, kSym, l, &
                     LFULL(2), LKsh, LKshp, LREAD, lSh, lSym, LWORK, MaxB, MaxRedT, MaxVecPerBatch, mDen, Mmax, mrs, mSh, mTvec, &
                     MUSED, MxBasSh, n1, n2, nBatch, nBs, nMat, nnO, nRS, nT1, nT2, NumCV, numSh, NUMV, NumVT, nVec, nVrs, nOrb(8,2)
real(kind=wp) :: Fact, fcorr, LKThr, SKsh, tau(2), TCC1, TCC2, tcoul(2), TCR1, TCR2, TCS1, TCS2, TCT1, TCT2, TCX1, TCX2, texch(2), &
                 thrv(2), tmotr(2), Tmp, TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), tscrn(2), TWC1, TWC2, &
                 TWR1, TWR2, TWS1, TWS2, TWT1, TWT2, TWX1, TWX2, xFab, xTmp, YMax, YshMax
logical(kind=iwp) :: add, DoScreen
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: Debug
#endif
character(len=50) :: CFmt
character :: mode, mode2
type(NDSBA_type) :: DiaH
type(L_Full_Type) :: L_Full
type(Lab_Type) :: Lab
integer(kind=iwp), allocatable :: Indx(:,:), iShp_rs(:), kOffSh(:,:), nnBfShp(:,:)
real(kind=wp), allocatable :: AbsC(:), Diag(:), Drs(:), Faa(:), Fia(:), Frs(:), Lrs(:,:), MLk(:,:), SumAClk(:,:), SvShp(:,:), &
                              VJ(:), Ylk(:,:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: myJRED1, NNBSTMX, ntv0
real(kind=wp), allocatable :: DiagJ(:)
#endif
real(kind=wp), parameter :: FactCI = One
character(len=*), parameter :: SECNAM = 'CHO_LK_SCF'
integer(kind=iwp), external :: Cho_LK_MaxVecPerBatch
real(kind=wp), external :: Cho_LK_ScreeningThreshold, ddot_

!***********************************************************************
#ifdef _DEBUGPRINT_
Debug = .false. ! to avoid double printing in SCF-debug
#endif

IREDC = -1 ! unknown reduced set

iLoc = 3 ! use scratch location in reduced index arrays

if ((nDen /= 1) .and. (nDen /= 2)) then
  write(u6,*) SECNAM//'Invalid parameter nDen= ',nDen
  call abend()
end if

call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = zero ! time read/transform vectors
tcoul(:) = zero ! time for computing Coulomb
texch(:) = zero ! time for computing Exchange
tmotr(:) = zero ! time for the half-transf of vectors
tscrn(:) = zero ! time for screening overhead

! ==================================================================

! Various offsets
! ----------------
MaxB = nBas(1)
do ISYM=2,NSYM
  MaxB = max(MaxB,nBas(iSym))
end do

!*************************************************
nnO = 0
nT1 = 0
do jDen=1,nDen

  do ISYM=1,NSYM

    nOrb(iSym,jDen) = nForb(iSym,jDen)+nIorb(iSym,jDen)

    kOff(iSym,jDen) = nnO

    nnO = nnO+nOrb(iSym,jDen)

  end do

  nT1 = nT1+nnO*(2-jDen) ! tot # alpha orbitals

end do

nT2 = nnO-nT1

! Define max number of vectors to be treated in core at once

MaxVecPerBatch = Cho_LK_MaxVecPerBatch()

! Define the basic screening threshold

LKThr = Cho_LK_ScreeningThreshold(dFmat)

! Adjust the damping according to the Abs(Max offDiag FMOmat)
!tbp, may 2013: adjustment moved to Cho_LK_ScreeningThreshold
fcorr = dmpk
!tbp if (dFmat > zero) then
!tbp   if (dFmat < 1.0e3_wp*LKThr) then
!tbp     fcorr = dmpk*1.0e-2_wp
!tbp   end if
!tbp   if (dFmat <= LKThr) fcorr = fcorr*1.0e-2_wp
!tbp end if

tau(1) = (LKThr/max(1,nT1))*fcorr ! screening alpha Fock matrix
tau(2) = (LKThr/max(1,nT2))*fcorr ! screening beta Fock matrix

MaxRedT = MaxRed
call GAIGOP_SCAL(MaxRedT,'+')
if (Estimate) then
  do i=1,nDen
    tau(i) = tau(i)/MaxRedT
  end do
end if

!xtau(1) = sqrt(tau(1))
!xtau(2) = sqrt(tau(2))

! Vector MO transformation screening thresholds
NumVT = NumChT
call GAIGOP_SCAL(NumVT,'+')
thrv(1) = (sqrt(LKThr/(max(1,nT1)*NumVT)))*fcorr
thrv(2) = (sqrt(LKThr/(max(1,nT2)*NumVT)))*fcorr

call mma_allocate(DIAG,NNBSTRT(1),Label='DIAG')
DIAG(:) = Zero

#ifdef _MOLCAS_MPP_
if ((nProcs > 1) .and. Update .and. Is_Real_Par()) then
  NNBSTMX = 0
  do i=1,nSym
    NNBSTMX = max(NNBSTMX,NNBSTR(i,1))
  end do
  call mma_allocate(DiagJ,NNBSTMX,Label='DiagJ')
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
call mma_allocate(Ylk,MaxB,nnO,Label='Ylk')

! allocate memory for the ML[k] list of largest elements in significant shells
call mma_allocate(MLk,nShell,nnO,Label='MLk')

! allocate memory for the list of  S:= sum_l abs(C(l)[k]) for each shell
call mma_allocate(SumAClk,nShell,nnO,Label='SumAClk')

! allocate memory for the Index array
call mma_allocate(Indx,[0,nShell],[1,nnO],Label='Indx')

! allocate memory for kOffSh
call mma_allocate(kOffSh,nShell,nSym,Label='kOffSh')

! allocate memory for nnBfShp
call mma_allocate(nnBfShp,nnShl_tot,nSym,Label='nnBfShp')

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

! allocate memory for the Diagonal of the Fock matrix
call mma_allocate(Fia,MxBasSh,Label='Fia')
call mma_allocate(Faa,nShell,Label='Faa')
Fia(:) = Zero
Faa(:) = Zero

! *** Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)
do jDen=1,nDen
  do kSym=1,nSym
    do jK=1,nOrb(kSym,jDen)
      jK_a = jK+kOff(kSym,jDen)

      do iaSh=1,nShell

        SKsh = zero
        iS = kOffSh(iaSh,kSym)+1
        iE = kOffSh(iaSh,kSym)+nBasSh(kSym,iaSh)
        do ik=iS,iE
          SKsh = SKsh+POrb(jDen)%SB(kSym)%A2(ik,jK)**2
        end do

        SumAClk(iaSh,jK_a) = SKsh

      end do

    end do
  end do
end do

! *** Compute Shell-pair Offsets in the K-matrix

do iSyma=1,nSym

  LKshp = 0

  do iaSh=1,nShell

    do ibSh=1,iaSh

      iShp = iaSh*(iaSh-1)/2+ibSh

      nnBfShp(iShp,iSyma) = LKShp

      LKShp = LKShp+nBasSh(iSyma,iaSh)*nBasSh(iSyma,ibSh)-(1-min((iaSh-ibSh),1))*nBasSh(iSyma,iaSh)*(nBasSh(iSyma,iaSh)-1)/2

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

  ! ****************     MEMORY MANAGEMENT SECTION    *****************
  ! --------------------------------------------------------------
  ! compute memory needed to store at least 1 vector of JSYM
  ! half-transformed for a given MO-index k  --> La,[k]
  ! --------------------------------------------------------------
  mTvec = 0 ! mem for storing the half-transformed vec

  do l=1,nSym
    k = Mul(l,JSYM)
    Mmax = 0
    do jDen=1,nDen
      Mmax = max(Mmax,nForb(k,jDen)+nIorb(k,jDen))
    end do
    if (Mmax > 0) mTvec = max(mTvec,nBas(l))
  end do

  mTvec = max(mTvec,1)

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
        write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
        write(u6,*) 'nVrs=',nVrs
        call Abend()
      end if

      call Cho_X_SetRed(irc,iLoc,JRED)
      ! set index arrays at iLoc
      if (irc /= 0) then
        write(u6,*) SECNAM//'cho_X_setred non-zero return code.   rc= ',irc
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

      nVec = min((LWORK-LFULL(2))/(nRS+mTvec+LFULL(1)),min(nVrs,MaxVecPerBatch))

      if (nVec < 1) then
        write(u6,*) SECNAM//': Insufficient memory for batch'
        write(u6,*) ' LWORK= ',LWORK
        write(u6,*) ' jsym= ',jsym
        write(u6,*) ' min. mem. need= ',nRS+mTvec+LFULL(1)
        write(u6,*) ' nRS = ',nRS
        write(u6,*) ' mTvec = ',mTvec
        write(u6,*) ' LFULL = ',LFULL
        call Quit(_RC_MEMORY_ERROR_)
        nBatch = -9999 ! dummy assignment
      end if

      LREAD = nRS*nVec

      call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

      if (JSYM == 1) then
        ! Transform the density to reduced storage
        add = .false.
        nMat = 1
        call swap_full2rs(irc,iLoc,nRS,nMat,JSYM,PLT,Drs,add)
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

        if ((NUMV <= 0) .or. (NUMV /= JNUM)) then
          rc = 77
          return
        end if

        call CWTIME(TCR2,TWR2)
        tread(1) = tread(1)+(TCR2-TCR1)
        tread(2) = tread(2)+(TWR2-TWR1)

        if (JSYM == 1) then
          ! ************ (alpha+beta) COULOMB CONTRIBUTION  ****************
          !
          ! Contraction with the density matrix
          ! -----------------------------------
          ! V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * P(rs)
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

        end if  ! Coulomb contribution

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
        mDen = 1
        call Allocate_DT(Lab,JNUM,nBasSh,nBas,nShell,nSym,mDen)

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

        do jDen=1,nDen

          do kSym=1,nSym

            lSym = Mul(JSYM,kSym)

            do jK=1,nOrb(kSym,jDen)
              jK_a = jK+kOff(kSym,jDen)

              Lab%A0(1:nBas(lSym)*JNUM) = Zero

              if (DoScreen) then

                call CWTIME(TCS1,TWS1)
                !--------------------------------------------------------------
                ! Setup the screening
                !--------------------------------------------------------------
                do ik=1,nBas(kSym)
                  AbsC(ik) = abs(POrb(jDen)%SB(kSym)%A2(ik,jK))
                end do

                if (lSym >= kSym) then
                  !-----------------------------------------------------------
                  ! Y(l)[k] = sum_n  DH(l,n) * |C(n)[k]|
                  !===========================================================
                  Mode = 'N'
                  n1 = nBas(lSym)
                  n2 = nBas(kSym)

                else
                  !-----------------------------------------------------------
                  ! Y(l)[k] = sum_n  DH(n,l) * |C(n)[k]|
                  !===========================================================
                  Mode = 'T'
                  n1 = nBas(kSym)
                  n2 = nBas(lSym)

                end if

                if (n1 > 0) call DGEMV_(Mode,n1,n2,ONE,DiaH%SB(lSym,kSym)%A2,n1,AbsC,1,ZERO,Ylk(1,jK_a),1)

                ! List the shells present in Y(l)[k] by the largest element
                do ish=1,nShell
                  YshMax = zero
                  do ibs=1,nBasSh(lSym,ish)
                    YshMax = max(YshMax,Ylk(koffSh(ish,lSym)+ibs,jK_a))
                  end do
                  MLk(ish,jK_a) = YshMax
                end do

                ! Sort the list ML[k]
                numSh = 0 ! # of significant shells
                jml = 1

                do ish=1,nShell
                  Indx(ish,jK_a) = ish
                end do

                do while (jml <= nShell)

                  YMax = MLk(jml,jK_a)
                  jmlmax = jml

                  do iml=jml+1,nShell ! get the max
                    if (MLk(iml,jK_a) > YMax) then
                      YMax = MLk(iml,jK_a)
                      jmlmax = iml
                    end if
                  end do

                  if (jmlmax /= jml) then ! swap positions
                    xTmp = MLk(jml,jK_a)
                    iTmp = Indx(jml,jK_a)
                    MLk(jml,jK_a) = YMax
                    Indx(jml,jK_a) = Indx(jmlmax,jK_a)
                    MLk(jmlmax,jK_a) = xTmp
                    Indx(jmlmax,jK_a) = iTmp
                  end if

                  ! Exact bounds (quadratic scaling of the MO transformation)
                  !tbp, may 2013: reintroduce exact bound
                  if (MLk(jml,jK_a)*MLk(1,jK_a) >= tau(jDen)) then
                  ! Here we use a non-exact bound for the exchange matrix
                  ! The positive definiteness of the exchange matrix
                  ! combined with the structure of the density matrix makes this
                  ! bound acceptable and likely to be almost exact for what concerns
                  ! the exchange energy
                  !tbp if (MLk(jml,jK_a) >= xtau(jDen)) then
                    numSh = numSh+1
                  else
                    jml = nShell  ! exit the loop
                  end if

                  jml = jml+1

                end do

                Indx(0,jk_a) = numSh

                call CWTIME(TCS2,TWS2)
                tscrn(1) = tscrn(1)+(TCS2-TCS1)
                tscrn(2) = tscrn(2)+(TWS2-TWS1)
                !------------------------------------------------------------------
              end if ! Screening setup

              ! Transform vectors for shells in the list ML[k]
              !
              ! Screening based on the Frobenius norm: sqrt(sum_ij  A(i,j)^2)
              !
              !   || La,J[k] ||  <=  || Lab,J || * || Cb[k] ||

              call CWTIME(TCT1,TWT1)

              do iSh=1,Indx(0,jk_a)

                iaSh = Indx(iSh,jK_a)

                Lab%Keep(iaSh,1) = .true.

                ibcount = 0

                do ibSh=1,nShell

                  iOffShb = kOffSh(ibSh,kSym)

                  iShp = iTri(iaSh,ibSh)

                  if (iShp_rs(iShp) <= 0) cycle

                  if ((nnBstRSh(JSym,iShp_rs(iShp),iLoc)*nBasSh(lSym,iaSh)*nBasSh(kSym,ibSh) > 0) .and. &
                      (sqrt(abs(SumAClk(ibSh,jK_a)*SvShp(iShp_rs(iShp),1))) >= thrv(jDen))) then

                    ibcount = ibcount+1

                    if (lSym >= kSym) then

                      i1 = 1
                      if (iaSh < ibSh) i1 = 2

                      ! LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
                      ! ----------------------------------
                      Mode = 'N'
                      n1 = nBasSh(lSym,iaSh)*JNUM
                      n2 = nBasSh(kSym,ibSh)

                      call DGEMV_(Mode,n1,n2,ONE,L_Full%SPB(lSym,iShp_rs(iShp),i1)%A21,n1,POrb(jDen)%SB(kSym)%A2(iOffShb+1:,jK),1, &
                                  ONE,Lab%SB(iaSh,lSym,1)%A,1)

                    else ! lSym < kSym

                      i1 = 1
                      if (ibSh < iaSh) i1 = 2

                      ! LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
                      ! ----------------------------------
                      Mode = 'T'
                      n1 = nBasSh(kSym,ibSh)
                      n2 = JNUM*nBasSh(lSym,iaSh)

                      call DGEMV_(Mode,n1,n2,ONE,L_Full%SPB(kSym,iShp_rs(iShp),i1)%A12,n1,POrb(jDen)%SB(kSym)%A2(iOffShb+1:,jK),1, &
                                  ONE,Lab%SB(iaSh,lSym,1)%A,1)

                    end if

                  end if

                end do

                ! The following re-assignement is used later on to check if the
                ! iaSh vector LaJ[k] can be neglected because identically zero

                if (ibcount == 0) Lab%Keep(iaSh,1) = .false.

              end do

              call CWTIME(TCT2,TWT2)
              tmotr(1) = tmotr(1)+(TCT2-TCT1)
              tmotr(2) = tmotr(2)+(TWT2-TWT1)

              ! Prepare the J-screening

              call CWTIME(TCS1,TWS1)

              do iSh=1,Indx(0,jk_a)

                iaSh = Indx(iSh,jK_a)

                if (.not. Lab%Keep(iaSh,1)) cycle

                if (lSym >= kSym) then

                  ! Faa,[k] = sum_J  (LaJ[k])**2
                  ! ------------------------------
                  Inc = nBasSh(lSym,iaSh)
                  n1 = 1

                else ! lSym < kSym

                  ! Faa,[k] = sum_J  (LJa[k])**2
                  ! -----------------------------
                  Inc = 1
                  n1 = JNUM

                end if

                Tmp = Zero
                do ia=1,nBasSh(lSym,iaSh)
                  Fia(ia) = DDot_(JNUM,Lab%SB(iash,lSym,1)%A(1+n1*(ia-1):),Inc,Lab%SB(iash,lSym,1)%A(1+n1*(ia-1):),Inc)
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

              do lSh=1,Indx(0,jK_a)

                iaSh = Indx(lSh,jK_a)

                iaSkip = merge(1,0,Lab%Keep(iaSh,1))

                mSh = 1

                do while (mSh <= Indx(0,jK_a))

                  ibSh = Indx(mSh,jK_a)

                  ibSkip = merge(1,0,Lab%Keep(ibSh,1))

                  iShp = iTri(iaSh,ibSh)

                  iOffShb = kOffSh(ibSh,lSym)

                  iOffAB = nnBfShp(iShp,lSym)

                  xFab = sqrt(abs(Faa(iaSh)*Faa(ibSh)))

                  if (MLk(lSh,jK_a)*MLk(mSh,jK_a) < tau(jDen)) then

                    mSh = Indx(0,jK_a)  ! skip the rest

                  else if ((iaSh == ibSh) .and. (xFab >= tau(jDen)/MaxRedT) .and. (iaSkip == 1)) then

                    if (lSym >= kSym) then

                      ! F(a,b)[k] = F(a,b)[k] - FactXI * sum_J  L(a,J)[k] * L(b,J)[k]
                      ! --------------------------------------------------------------
                      Mode = 'N'
                      Mode2 = 'T'
                      n1 = nBasSh(lSym,iaSh)
                      nBs = nBasSh(lSym,iaSh)

                    else ! lSym < kSym

                      ! F(a,b)[k] = F(a,b)[k] - FactXI * sum_J  L(J,a)[k] * L(J,b)[k]
                      ! --------------------------------------------------------------
                      Mode = 'T'
                      Mode2 = 'N'
                      n1 = JNUM
                      nBs = nBasSh(lSym,iaSh)

                    end if

                    call DGEMM_Tri(Mode,Mode2,nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,-FActXI,Lab%SB(iaSh,lSym,1)%A,n1, &
                                   Lab%SB(ibSh,lSym,1)%A,n1,ONE,KLT(jDen)%SB(lSym)%A1(iOffAB+1:),nBs)

                  else if ((iaSh > ibSh) .and. (xFab >= tau(jDen)/MaxRedT) .and. (iaSkip*ibSkip == 1)) then

                    if (lSym >= kSym) then

                      ! F(a,b)[k] = F(a,b)[k] - FactXI * sum_J  L(a,J)[k] * L(b,J)[k]
                      ! --------------------------------------------------------------
                      Mode = 'N'
                      Mode2 = 'T'
                      nBs = nBasSh(lSym,iaSh)
                      n1 = nBasSh(lSym,iaSh)
                      n2 = nBasSh(lSym,ibSh)

                    else ! lSym < kSym

                      ! F(a,b)[k] = F(a,b)[k] - FactXI * sum_J  L(J,a)[k] * L(J,b)[k]
                      ! --------------------------------------------------------------
                      Mode = 'T'
                      Mode2 = 'N'
                      n1 = JNUM
                      n2 = JNUM
                      nBs = nBasSh(lSym,iaSh)

                    end if

                    call DGEMM_(Mode,Mode2,nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,-FactXI,Lab%SB(iaSh,lSym,1)%A,n1, &
                                Lab%SB(ibSh,lSym,1)%A,n2,ONE,KLT(jDen)%SB(lSym)%A1(iOffAB+1:),nBs)

                  end if

                  mSh = mSh+1 ! update shell counter

                end do

              end do

              call CWTIME(TCX2,TWX2)
              texch(1) = texch(1)+(TCX2-TCX1)
              texch(2) = texch(2)+(TWX2-TWX1)

            end do ! loop over k MOs

          end do ! loop over MOs symmetry

        end do ! loop over densities

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

        ! --------------------------------------------------------------------
        ! --------------------------------------------------------------------

      end do ! end batch loop

      if (JSYM == 1) then
        ! backtransform fock matrix to full storage
        add = .true.
        nMat = 1
        do jDen=1,nDen
          call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,FLT(jDen),Frs,add)
        end do
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
      DoScreen = (JRED < myJRED1) .or. (ntv0 >= Nscreen)
      if (DoScreen) ntv0 = 0
    end if
#   endif

  end do ! loop over red sets

end do ! loop over JSYM

! Accumulate Coulomb and Exchange contributions
do jDen=1,nDen

  do iSym=1,nSym

    do iaSh=1,nShell

      ioffa = kOffSh(iaSh,iSym)

      ! ibSh < iaSh
      ! -----------
      do ibSh=1,iaSh-1

        iShp = iaSh*(iaSh-1)/2+ibSh

        iOffAB = nnBfShp(iShp,iSym)

        ioffb = kOffSh(ibSh,iSym)

        do ib=1,nBasSh(iSym,ibSh)

          do ia=1,nBasSh(iSym,iaSh)

            iab = nBasSh(iSym,iaSh)*(ib-1)+ia

            iag = ioffa+ia
            ibg = ioffb+ib

            iabg = iTri(iag,ibg)

            FLT(jDen)%sb(iSym)%A1(iabg) = FLT(jDen)%sb(iSym)%A1(iabg)+KLT(jDen)%SB(iSym)%A1(iOffAB+iab)

          end do

        end do

      end do

      ! ibSh = iaSh
      ! -----------
      iShp = iaSh*(iaSh+1)/2

      iOffAB = nnBfShp(iShp,iSym)

      do ib=1,nBasSh(iSym,iaSh)

        do ia=ib,nBasSh(iSym,iaSh)

          iab = ia*(ia-1)/2+ib

          iag = ioffa+ia
          ibg = ioffa+ib

          iabg = iag*(iag-1)/2+ibg

          FLT(jDen)%sb(iSym)%A1(iabg) = FLT(jDen)%sb(iSym)%A1(iabg)+KLT(jDen)%SB(iSym)%A1(iOffAB+iab)

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

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky SCF timing from '//SECNAM
  write(u6,CFmt) '------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'Fock matrix construction        CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ VECTORS                              ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'COULOMB                                   ',tcoul(1),tcoul(2)
  write(u6,'(2x,A26,2f10.2)') 'SCREENING OVERHEAD                        ',tscrn(1),tscrn(2)
  write(u6,'(2x,A26,2f10.2)') 'MO HALF-TRANSFORM VECTORS                 ',tmotr(1),tmotr(2)
  write(u6,'(2x,A26,2f10.2)') 'EXCHANGE                                  ',texch(1),texch(2)
  write(u6,*)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
if (Debug) then ! to avoid double printing in SCF-debug

  write(u6,'(6X,A)') 'TEST PRINT FROM '//SECNAM
  write(u6,'(6X,A)')
  write(u6,'(6X,A)') '***** FOCK MATRIX AO-BASIS ***** '
  do jDen=1,nDen
    if (nDen == 2) then
      if (jden == 1) write(u6,'(6X,A)') '******** ALPHA SPIN ******** '
      if (jden == 2) write(u6,'(6X,A)') '******** BETA SPIN ********* '
    end if
    do ISYM=1,NSYM
      if (NBAS(ISYM) > 0) then
        write(u6,'(6X,A)')
        write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES: ',ISYM
        call TRIPRT('','',FLT(jDen)%SB(iSym)%A1,NBAS(ISYM))
      end if
    end do
  end do

end if
#endif

rc = 0

return

end subroutine CHO_LK_SCF
