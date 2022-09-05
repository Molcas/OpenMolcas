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

subroutine CHO_FOCK_RASSI_X(DLT,MO1,MO2,FLT,FSQ,TUVX)
!*********************************************************************
!  Author : F. Aquilante
!
!  Note:  this routine differs from CHO_FOCK_RASSI because it can
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

use ChoArr, only: nDimRS
use ChoSwp, only: InfVec
use Symmetry_Info, only: Mul
use Fock_util_global, only: Fake_CMO2
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, SBA_Type, twxy_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
type(DSBA_Type), intent(in) :: DLT, MO1(2), MO2(2)
type(DSBA_Type), intent(inout) :: FLT(1), FSQ
real(kind=wp), intent(_OUT_) :: TUVX(*)
#include "chotime.fh"
#include "rassi.fh"
#include "cholesky.fh"
#include "choorb.fh"
integer(kind=iwp) :: i, ia, iabt, ib, iBatch, iCase, iLoc, irc, iREDC, iSkip(8), iSwap, iSym, iSyma, iSymb, iSymk, iSymv, IVEC2, &
                     iVrs, jDen, JNUM, JRED, JRED1, JRED2, jSym, JVC, JVEC, k, kDen, kMOs, l, LREAD, LWORK, mDen, mTTvec, mTvec, &
                     MUSED, NAv, NAw, nBatch, nDen, NK, nMOs, nRS, NUMV, nVec, nVrs, rc
real(kind=wp) :: Fact, TCC1, TCC2, TCINT1, TCINT2, tcoul(2), TCR1, TCR2, TCR3, TCR4, TCR7, TCX1, TCX2, texch(2), tintg(2), TOTCPU, &
                 TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), TWC1, TWC2, TWINT1, TWINT2, TWR1, TWR2, TWR3, TWR4, &
                 TWR7, TWX1, TWX2
logical(kind=iwp) :: add, DoReord
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: Debug
#endif
character(len=50) :: CFmt
type(SBA_Type), target :: Laq(2)
type(twxy_Type) :: Scr
real(kind=wp), allocatable :: Drs(:), Frs(:), Lrs(:,:)
real(kind=wp), pointer :: VJ(:) => null()
real(kind=wp), parameter :: FactCI = One, FactXI = -One
character(len=*), parameter :: SECNAM = 'CHO_FOCK_RASSI_X'
logical(kind=iwp), parameter :: DoRead = .false.

!*************************************************
#ifdef _DEBUGPRINT_
Debug = .false. ! to avoid double printing in CASSCF-debug
#endif
DoReord = .false.
IREDC = -1 ! unknown reduced set in core

nDen = 2
if (Fake_CMO2) nDen = 1 ! MO1 = MO2
kDen = nDen

call CWTIME(TOTCPU1,TOTWALL1) ! start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = zero ! time read/transform vectors
tcoul(:) = zero ! time for computing Coulomb
texch(:) = zero ! time for computing Exchange
tintg(:) = zero ! time for computing (tw|xy) integrals

! *************** BIG LOOP OVER VECTORS SYMMETRY *******************
do jSym=1,nSym

  if (NumCho(jSym) < 1) cycle

  iCase = 0 ! twxy
  call Allocate_DT(Scr,nAsh,nAsh,JSYM,nSym,iCase)

  iLoc = 3 ! use scratch location in reduced index arrays

  ! ****************     MEMORY MANAGEMENT SECTION    *****************
  ! --------------------------------------------------------------
  ! compute memory needed to store at least 1 vector of JSYM
  ! and do all the subsequent calculations
  ! --------------------------------------------------------------
  mTvec = 0  ! mem for storing the half-transformed vec
  mTTvec = 0  ! mem for Lvb,J and Lvw,J
  do l=1,nSym
    k = Mul(l,JSYM)
    mTvec = mTvec+nDen*nBas(l)*nIsh(k)
    mTTvec = mTTvec+(nBas(l)+nAsh(l))*nAsh(k)
  end do

  mTvec = max(mTvec,mTTvec,1)

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

  JRED1 = InfVec(1,2,jSym)            ! red set of the 1st vec
  JRED2 = InfVec(NumCho(jSym),2,jSym) ! red set of the last vec

  do JRED=JRED1,JRED2

    call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

    if (nVrs == 0) cycle ! no vectors in that (jred,jsym)

    if (nVrs < 0) then
      write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
      call abend()
    end if

    call Cho_X_SetRed(irc,iLoc,JRED) ! set index arrays at iLoc
    if (irc /= 0) then
      write(u6,*) SECNAM//'cho_X_setred non-zero return code.   rc= ',irc
      call abend()
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

    nVec = min(LWORK/(nRS+mTvec),nVrs)

    if (nVec < 1) then
      write(u6,*) SECNAM//': Insufficient memory for batch'
      write(u6,*) 'LWORK= ',LWORK
      write(u6,*) 'min. mem. need= ',nRS+mTvec
      write(u6,*) 'jsym= ',jsym
      rc = 33
      call Abend()
      nBatch = -9999 ! dummy assignment
    end if

    LREAD = nRS*nVec

    call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

    if (JSYM == 1) then
      ! Transform the density to reduced storage
      add = .false.
      mDen = 1
      call swap_full2rs(irc,iLoc,nRS,mDen,JSYM,[DLT],Drs,add)
    end if

    ! BATCH over the vectors ----------------------------

    nBatch = (nVrs-1)/nVec+1

    do iBatch=1,nBatch

      if (iBatch == nBatch) then
        JNUM = nVrs-nVec*(nBatch-1)
      else
        JNUM = nVec
      end if

      iSwap = 2 ! LpJ,b are returned
      do jDen=1,nDen
        call Allocate_DT(Laq(jDen),nIsh,nBas,nVec,JSYM,nSym,iSwap)
      end do

      JVEC = nVec*(iBatch-1)+iVrs
      IVEC2 = JVEC-1+JNUM

      call CWTIME(TCR1,TWR1)

      call CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,NUMV,IREDC,MUSED)

      if ((NUMV <= 0) .or. (NUMV /= JNUM)) then
        rc = 77
        write(u6,*) 'return code = ',rc
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
        ! V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * DI(rs)
        !==========================================================

        call CWTIME(TCC1,TWC1)

        VJ(1:JNUM) => Laq(1)%A0(1:JNUM)

        call DGEMV_('T',nRS,JNUM,ONE,Lrs,nRS,Drs,1,ZERO,VJ,1)

        ! FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
        !===========================================================

        Fact = real(min(jVec-iVrs,1),kind=wp)

        call DGEMV_('N',nRS,JNUM,FactCI,Lrs,nRS,VJ,1,Fact,Frs,1)

        call CWTIME(TCC2,TWC2)
        tcoul(1) = tcoul(1)+(TCC2-TCC1)
        tcoul(2) = tcoul(2)+(TWC2-TWC1)

        VJ => null()

      end if ! Coulomb contribution

      ! *************** EXCHANGE CONTRIBUTIONS  ***********************

      call CWTIME(TCR3,TWR3)

      kMOs = 1
      nMOs = nDen

      ! Set up the skipping flags
      ! ---------------------------------------------------------
      do i=1,nSym

        k = Mul(i,JSYM)
        iSkip(k) = min(1,nIsh(k)*NBAS(i))

      end do
      ! -------------------------------------------------------------

      ! *********************** HALF-TRANSFORMATION  ****************

      call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,MO1,Laq,DoRead)

      call CWTIME(TCR4,TWR4)
      tread(1) = tread(1)+(TCR4-TCR3)
      tread(2) = tread(2)+(TWR4-TWR3)

      if (irc /= 0) then
        rc = irc
        write(u6,*) 'CHO_X_getVtra failed! '
        return
      end if

      call CWTIME(TCX1,TWX1)

      do iSyma=1,nSym

        iSymk = Mul(JSYM,iSyma)

        ! ---------------------------------------------------------------------
        ! *** Compute the InActive exchange matrix
        !
        !     FI(ab) = FI(ab) + FactXI * sum_Jk  X(2)kJ,a * X(1)kJ,b
        ! ---------------------------------------------------------------------
        NK = nIsh(iSymk)

        if (iSkip(iSymk) /= 0) then

          call DGEMM_('T','N',NBAS(iSyma),NBAS(iSyma),NK*JNUM,FactXI,Laq(kDen)%SB(iSymk)%A3,NK*JNUM,Laq(1)%SB(iSymk)%A3,NK*JNUM, &
                      One,FSQ%SB(iSyma)%A2,NBAS(iSyma))

        end if

        ! --------------------------------------------------------------------
      end do  !loop over MOs symmetries

      call CWTIME(TCX2,TWX2)
      texch(1) = texch(1)+(TCX2-TCX1)
      texch(2) = texch(2)+(TWX2-TWX1)

      do jDen=1,nDen
        call Deallocate_DT(Laq(jDen))
      end do

      ! ************  END EXCHANGE CONTRIBUTION  ****************

      iSwap = 0  ! Lvb,J are returned
      call Allocate_DT(Laq(1),nAsh,nBas,nVec,JSYM,nSym,iSwap)
      call Allocate_DT(Laq(2),nAsh,nAsh,nVec,JSYM,nSym,iSwap)

      ! ----------------------------------------------------------------
      ! First half Active transformation  Lvb,J = sum_a  C1(v,a) * Lab,J
      ! ----------------------------------------------------------------

      call CWTIME(TCR7,TWR7)

      ! Set up the skipping flags
      ! ---------------------------------------------------------
      do i=1,nSym

        k = Mul(i,JSYM)
        iSkip(k) = min(1,NBAS(i)*nAsh(k))

      end do

      kMOs = 1
      nMOs = 1 ! Active MOs (1st set)

      call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,MO2,Laq,DoRead)

      if (irc /= 0) then
        rc = irc
        return
      end if

      ! ----------------------------------------------------------------
      ! Active-Active transformation  Lvw,J = sum_b  Lvb,J * C2(w,b)
      ! ----------------------------------------------------------------
      do iSymb=1,nSym

        iSymv = Mul(JSYM,iSymb)
        NAv = nAsh(iSymv)
        NAw = nAsh(iSymb) ! iSymb=iSymw

        if (NAv*NAw /= 0) then

          do JVC=1,JNUM

            call DGEMM_('N','T',NAv,NAw,NBAS(iSymb),One,Laq(1)%SB(iSymv)%A3(:,:,JVC),NAv,MO2(kDen)%SB(iSymb)%A2,NAw,Zero, &
                        Laq(2)%SB(iSymv)%A3(:,:,JVC),NAv)

          end do

        end if

      end do

      ! *************** EVALUATION OF THE (TW|XY) INTEGRALS ***********

      call CWTIME(TCINT1,TWINT1)

      DoReord = (JRED == JRED2) .and. (iBatch == nBatch)

      call CHO_rassi_twxy(irc,Scr,Laq(2),TUVX,nAsh,JSYM,JNUM,DoReord)

      call CWTIME(TCINT2,TWINT2)
      tintg(1) = tintg(1)+(TCINT2-TCINT1)
      tintg(2) = tintg(2)+(TWINT2-TWINT1)

      if (irc /= 0) then
        rc = irc
        return
      end if

      ! ---------------- END (TW|XY) EVALUATION -----------------------

      call Deallocate_DT(Laq(2))
      call Deallocate_DT(Laq(1))
    end do ! end batch loop

    if (JSYM == 1) then
      ! backtransform fock matrix to full storage
      add = .true.
      mDen = 1
      call swap_rs2full(irc,iLoc,nRS,mDen,JSYM,FLT,Frs,add)
    end if

    ! free memory
    call mma_deallocate(Lrs)

    if (JSYM == 1) then
      call mma_deallocate(Frs)
      call mma_deallocate(Drs)
    end if

  end do ! loop over red sets

  call Deallocate_DT(Scr)

end do ! loop over JSYM

! Accumulate Coulomb and Exchange contributions
do iSym=1,nSym

  do ia=1,nBas(iSym)
    do ib=1,ia-1
      iabt = ia*(ia-1)/2+ib
      FSQ%SB(iSym)%A2(ib,ia) = FSQ%SB(iSym)%A2(ib,ia)+FLT(1)%SB(iSym)%A1(iabt)
      FSQ%SB(iSym)%A2(ia,ib) = FSQ%SB(iSym)%A2(ia,ib)+FLT(1)%SB(iSym)%A1(iabt)
    end do
    iabt = ia*(ia+1)/2
    FSQ%SB(iSym)%A2(ia,ia) = FSQ%SB(iSym)%A2(ia,ia)+FLT(1)%SB(iSym)%A1(iabt)
  end do

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

  write(u6,'(2x,A26,2f10.2)') 'READ/TRANSFORM VECTORS                    ',tread(1),tread(2)
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
if (Debug) then ! to avoid double printing in RASSI-debug

  write(u6,'(6X,A)') 'TEST PRINT FROM '//SECNAM
  write(u6,'(6X,A)')
  do ISYM=1,NSYM
    if (NBAS(ISYM) > 0) then
      write(u6,'(6X,A)') '***** INACTIVE FOCK MATRIX ***** '
      write(u6,'(6X,A)')
      write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
      call CHO_OUTPUT(FSQ%SB(ISYM)%A2,1,NBAS(ISYM),1,NBAS(ISYM),NBAS(ISYM),NBAS(ISYM),1,u6)
    end if
  end do

end if
#endif

rc = 0

return

end subroutine CHO_FOCK_RASSI_X
