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

subroutine CHO_FMCSCF(rc,FLT,nForb,nIorb,nAorb,FactXI,DLT,DoActive,POrb,nChM,W_PWXY,CMO,ExFac)
!***********************************************************************
!  Author : F. Aquilante
!
!**************** INACTIVE AO-BASIS FOCK MATRIX ************************
!
! FI(ab) = 2 * sum_J  Lab,J * V(J)  -  sum_Jk  Lka,J * Lkb,J
!
!****************   ACTIVE AO-BASIS FOCK MATRIX ************************
!
! FA(ab) = sum_J  Lab,J * U(J)  -  0.5 * sum_Jw  Lwa,J * Lwb,J
!
!****************   (WA|XY) integrals     ******************************
!
! (WA|XY) = sum_J  L(wa,J) * L(xy,J)
!
!***********************************************************************
!
! V(J) = sum_gd  Lgd,J * DI(gd)
! U(J) = sum_gd  Lgd,J * DA(gd)
!
! a,b,g,d:  AO-index
! k:        MO-index   belonging to (Frozen+Inactive)
! u,w,x,y:  MO-indices belonging to (Active)
!
!***********************************************************************

use ChoArr, only: nDimRS
use ChoSwp, only: InfVec
use Symmetry_Info, only: Mul
use Data_structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, SBA_Type, twxy_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: rc
type(DSBA_Type), intent(inout) :: FLT(2)
integer(kind=iwp), intent(in) :: nForb(8), nIorb(8), nAorb(8), nChM(8)
real(kind=wp), intent(in) :: FactXI, ExFac
type(DSBA_Type), intent(in) :: DLT(2), POrb(3), CMO
logical(kind=iwp), intent(in) :: DoActive
real(kind=wp), intent(_OUT_) :: W_PWXY(*)
#include "chotime.fh"
#include "cholesky.fh"
#include "choorb.fh"
integer(kind=iwp) :: i, iBatch, iCase, iLoc, irc, IREDC, iSkip(8), iSwap, iSyma, iSymb, iSymk, iSymv, iSymw, IVEC2, iVrs, JNUM, &
                     JRED, JRED1, JRED2, jSym, JVC, JVEC, k, kMOs, l, LREAD, LWORK, mTvec, mTvec1, mTvec2, mTvec3, mTvec4, MUSED, &
                     NAch, nAux(8), nAv, nAw, nBatch, nDen, NK, nMOs, nnA(8,8), nPorb(8), nRS, NUMV, nVec, nVrs
real(kind=wp) :: TCC1, TCC2, TCINT1, TCINT2, tcoul(2), TCR1, TCR2, TCR3, TCR4, TCR5, TCR6, TCR7, TCR8, TCX1, TCX2, TCX3, TCX4, &
                 texch(2), tintg(2), TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), TWC1, TWC2, TWINT1, TWINT2, &
                 TWR1, TWR2, TWR3, TWR4, TWR5, TWR6, TWR7, TWR8, TWX1, TWX2, TWX3, TWX4, xfac
logical(kind=iwp) :: add, DoRead, DoTraInt
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: Debug
#endif
character(len=50) :: CFmt
type(SBA_Type), target :: Laq(3), Lxy
type(twxy_type) :: Scr
real(kind=wp), allocatable :: Lrs(:,:), Drs(:,:), Frs(:,:)
real(kind=wp), pointer :: VJ(:) => null()
real(kind=wp), parameter :: FactCI = One, FactCA = One, FactXA = -Half
character(len=*), parameter :: SECNAM = 'CHO_FMCSCF'

!***********************************************************************
#ifdef _DEBUGPRINT_
Debug = .false. ! to avoid double printing in CASSCF-debug
#endif
if (ExFac /= One) then
  write(u6,*) 'WARNING: if you are running MCPDFT calculations'
  write(u6,*) 'and end up with this message, you are in trouble.'
  write(u6,*) 'Tweaks are needed!! Please contact'
  write(u6,*) 'Giovanni Li Manni'
  write(u6,*) 'giovannilimanni@gmail.com'
end if
DoRead = .false.
DoTraInt = .false.
IREDC = -1 ! unknown reduced set in core

nDen = 1
if (DoActive) nDen = 2

call CWTIME(TOTCPU1,TOTWALL1) ! start clock for total time

do i=1,2          ! 1 --> CPU   2 --> Wall
  tread(i) = zero ! time read/transform vectors
  tcoul(i) = zero ! time for computing Coulomb
  texch(i) = zero ! time for computing Exchange
  tintg(i) = zero ! time for computing (pu|vx) integrals
end do

! Define MOs in the Primary space
! ( Frozen + Inactive + Active )
! -------------------------------
do i=1,nSym
  nPorb(i) = nForb(i)+nIorb(i)+nAorb(i)
end do

call set_nnA(nSym,nAorb,nnA)

! ==================================================================

iLoc = 3 ! use scratch location in reduced index arrays

! *************** BIG LOOP OVER VECTORS SYMMETRY *******************
do jSym=1,nSym

  if (NumCho(jSym) < 1) cycle

  iCase = 1 ! (wa|xy)
  call Allocate_DT(Scr,nAorb,nBas,JSYM,nSym,iCase)

  ! --- Set up the skipping flags --------
  ! -------------------------------------------------------------
  do i=1,nSym
    k = Mul(i,JSYM)
    iSkip(i) = min(1,nBas(i)*nBas(k)) ! skip Lik vector
    iSkip(i) = iSkip(i)*(nPorb(i)+nChM(i))
  end do
  ! -------------------------------------------------------------

  ! ****************     MEMORY MANAGEMENT SECTION    *****************
  ! --------------------------------------------------------------
  ! compute memory needed to store at least 1 vector of JSYM
  ! and do all the subsequent calculations
  ! --------------------------------------------------------------
  mTvec1 = 0 ! mem for storing the half-transformed vec
  mTvec2 = 0 ! mem for storing the half-transformed vec
  mTvec3 = 0 ! mem for storing the half-transformed vec
  mTvec4 = 0 ! mem for storing the half-transformed vec

  nAux(:) = nForb(:)+nIorb(:)
  do l=1,nSym
    k = Mul(l,JSYM)
    mTvec1 = mTvec1+nBas(k)*nAux(l)
    mTvec2 = mTvec2+nBas(k)*nChM(l)
    mTvec3 = mTvec3+nBas(k)*nAorb(l)
    mTvec4 = mTvec4+nnA(k,l)
  end do

  mTvec = max(mTvec1,mTvec2,mTvec3+mTvec4,1)

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
      write(u6,*) SECNAM//'cho_X_setred non-zero return code.    rc= ',irc
      call abend()
    end if

    IREDC = JRED

    nRS = nDimRS(JSYM,JRED)

    if (JSYM == 1) then

      if (DoActive) then
        call mma_allocate(Drs,nRS,2,Label='Drs')
        call mma_allocate(Frs,nRS,2,Label='Frs')
      else
        call mma_allocate(Drs,nRS,1,Label='Drs')
        call mma_allocate(Frs,nRS,1,Label='Frs')
      end if

    end if

    call mma_maxDBLE(LWORK)

    nVec = min(LWORK/(nRS+mTvec),nVrs)

    if (nVec < 1) then
      write(u6,*) SECNAM//': Insufficient memory for batch'
      write(u6,*) 'LWORK= ',LWORK
      write(u6,*) 'min. mem. need= ',nRS+mTvec
      write(u6,*) 'reading ',nRS,' and transforming to ',mTvec
      write(u6,*) 'of jsym= ',jsym,' and JRED= ',JRED
      rc = 33
      call Abend()
      nBatch = -9999 ! dummy assignment
    end if

    LREAD = nRS*nVec

    call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

    if (JSYM == 1) then
      ! Transform the density to reduced storage
      add = .false.
      call swap_full2rs(irc,iLoc,nRS,nDen,JSYM,DLT,Drs,add)
    end if

    ! BATCH over the vectors ----------------------------

    nBatch = (nVrs-1)/nVec+1

    do iBatch=1,nBatch
      iSwap = 2 ! LpJ,b are returned
      call Allocate_DT(Laq(1),nAux,nBas,nVec,JSYM,nSym,iSwap)
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

        ! ************ INACTIVE COULOMB CONTRIBUTION  ****************
        !
        ! Contraction with the density matrix
        ! -----------------------------------
        ! V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * DI(rs)
        !==========================================================

        call CWTIME(TCC1,TWC1)

        VJ(1:JNUM) => Laq(1)%A0(1:JNUM)

        call DGEMV_('T',nRS,JNUM,ONE,Lrs,nRS,Drs(:,1),1,ZERO,VJ,1)

        ! FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
        !=============================================================

        xfac = real(min(jVec-iVrs,1),kind=wp)

        call DGEMV_('N',nRS,JNUM,FactCI,Lrs,nRS,VJ,1,xfac,Frs(:,1),1)

        if (DoActive) then
          ! ************ ACTIVE COULOMB CONTRIBUTION  ****************
          !
          ! Computing the intermediate vector V(J)
          ! -----------------------------------
          ! U{#J} <- U{#J}  +  sum_rs  L(rs,{#J}) * DA(rs)
          !==========================================================

          call DGEMV_('T',nRS,JNUM,ONE,Lrs,nRS,Drs(:,2),1,ZERO,VJ,1)

          ! FA(rs){#J} <- FA(rs){#J} + FactCA * sum_J L(rs,{#J})*U{#J}
          !===========================================================

          xfac = real(min(jVec-iVrs,1),kind=wp)

          call DGEMV_('N',nRS,JNUM,FactCA,Lrs,nRS,VJ,1,xfac,Frs(:,2),1)

        end if

        call CWTIME(TCC2,TWC2)
        tcoul(1) = tcoul(1)+(TCC2-TCC1)
        tcoul(2) = tcoul(2)+(TWC2-TWC1)

        VJ => null()

      end if  ! Coulomb (jsym=1)

      ! ************ BEGIN EXCHANGE CONTRIBUTIONS  ****************

      ! ***************** INACTIVE HALF-TRANSFORMATION  ****************

      kMOs = 1  ! inactive MOs
      nMOs = 1

      call CWTIME(TCR3,TWR3)

      call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,POrb,Laq,DoRead)

      call CWTIME(TCR4,TWR4)
      tread(1) = tread(1)+(TCR4-TCR3)
      tread(2) = tread(2)+(TWR4-TWR3)

      if (irc /= 0) then
        rc = irc
        return
      end if

      call CWTIME(TCX1,TWX1)

      do iSyma=1,nSym

        iSymk = Mul(JSYM,iSyma)

        ! ---------------------------------------------------------------------
        ! *** Compute only the LT part of the InActive exchange matrix ********
        !
        ! FI(ab) = FI(ab) + FactXI * sum_Jk  LkJ,a * LkJ,b
        ! ---------------------------------------------------------------------
        NK = nForb(iSymk)+nIorb(iSymk)

        if (iSkip(iSymk)*NK /= 0) then

          call DGEMM_TRI('T','N',nBas(iSyma),nBas(iSyma),NK*JNUM,FactXI,Laq(1)%SB(iSymk)%A3,NK*JNUM,Laq(1)%SB(iSymk)%A3,NK*JNUM, &
                         One,FLT(1)%SB(iSyma)%A1,nBas(iSyma))

          !write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYMA
          !call TRIPRT('FI: ',' ',FLT(1)%SB(iSyma)%A1,nBas(iSyma))

        end if

        ! --------------------------------------------------------------------
      end do  !loop over MOs symmetries

      call CWTIME(TCX2,TWX2)
      texch(1) = texch(1)+(TCX2-TCX1)
      texch(2) = texch(2)+(TWX2-TWX1)

      call Deallocate_DT(Laq(1))

      if (DoActive) then
        iSwap = 2  ! LxJ,b are returned
        call Allocate_DT(Laq(2),nChM,nBas,nVec,JSYM,nSym,iSwap)
        ! *********************** "CHOLESKY" HALF-TRANSFORMATION  ****************
        ! ----------------------------------------------------------------
        ! Using "Cholesky MOs" obtained by cholesky decomposing DA
        ! ----------------------------------------------------------------

        call CWTIME(TCR5,TWR5)

        kMOs = 2 ! Cholesky MOs
        nMOs = 2

        call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,POrb,Laq,DoRead)

        call CWTIME(TCR6,TWR6)
        tread(1) = tread(1)+(TCR6-TCR5)
        tread(2) = tread(2)+(TWR6-TWR5)

        if (irc /= 0) then
          rc = irc
          return
        end if

        call CWTIME(TCX3,TWX3)

        ! ---------------------------------------------------------------------
        ! *** Compute only the LT part of the Active exchange matrix **********
        !
        !     FA(ab) = FA(ab) + FactXA * sum_wJ  LwJ,a * LwJ,b
        ! ---------------------------------------------------------------------
        do iSyma=1,nSym

          iSymw = Mul(JSYM,iSyma)

          NAch = nChM(iSymw)

          if (iSkip(iSymw)*NAch /= 0) then

            call DGEMM_TRI('T','N',nBas(iSyma),nBas(iSyma),NAch*JNUM,FactXA,Laq(2)%SB(iSymw)%A3,NAch*JNUM,Laq(2)%SB(iSymw)%A3, &
                           NAch*JNUM,One,FLT(2)%SB(iSyma)%A1,nBas(iSyma))

            !write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYMA
            !call TRIPRT('FA: ',' ',FLT(2)%SB(iSyma)%A1,nBas(iSyma))

          end if

          ! --------------------------------------------------------------------
        end do ! loop over MOs symmetries

        call CWTIME(TCX4,TWX4)
        texch(1) = texch(1)+(TCX4-TCX3)
        texch(2) = texch(2)+(TWX4-TWX3)

        call Deallocate_DT(Laq(2))
      end if ! Do Active Exchange
      ! ************  END EXCHANGE CONTRIBUTIONS  ****************

      ! ----------------------------------------------------------------
      ! First half Active transformation  Lvb,J = sum_a  C(v,a) * Lab,J
      ! ----------------------------------------------------------------
      ! Lvw,J, LT-storage for the diagonal symmetry blocks
      iSwap = 4
      call Allocate_DT(Lxy,nAorb,nAorb,nVec,JSYM,nSym,iSwap)

      iSwap = 0 ! Lvb,J are returned
      call Allocate_DT(Laq(3),nAorb,nBas,nVec,JSYM,nSym,iSwap)

      call CWTIME(TCR7,TWR7)

      kMOs = 3 ! Active MOs
      nMOs = 3 ! Active MOs

      call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,POrb,Laq,DoRead)

      if (irc /= 0) then
        rc = irc
        return
      end if

      ! ----------------------------------------------------------------
      ! Active-Active transformation  Lvw,J = sum_b  Lvb,J * C(w,b)
      ! ----------------------------------------------------------------
      if (JSYM == 1) then ! Lvw,J in LT-storage

        do iSyma=1,nSym

          NAv = nAorb(iSyma)

          if (NAv /= 0) then

            do JVC=1,JNUM

              call DGEMM_Tri('N','T',NAv,NAv,nBas(iSyma),One,Laq(3)%SB(iSyma)%A3(:,:,JVC),NAv,POrb(3)%SB(iSyma)%A2,NAv,Zero, &
                             Lxy%SB(iSyma)%A2(:,JVC),NAv)

            end do

          end if

        end do

      else

        ! ----------------------------------------------------------------
        ! Active-Active transformation  Lvw,J = sum_b  Lvb,J * C(w,b)
        ! ----------------------------------------------------------------
        do iSymb=1,nSym

          iSymv = Mul(JSYM,iSymb)
          NAv = nAorb(iSymv)
          NAw = nAorb(iSymb) ! iSymb=iSymw

          if ((NAv*NAw /= 0) .and. (iSymv < iSymb)) then

            do JVC=1,JNUM

              call DGEMM_('N','T',NAv,NAw,nBas(iSymb),One,Laq(3)%SB(iSymv)%A3(:,:,JVC),NAv,POrb(3)%SB(iSymb)%A2,NAw,Zero, &
                          Lxy%SB(iSymv)%A2(:,JVC),NAv)

            end do

          end if

        end do

      end if ! jSym >=< 1

      call CWTIME(TCR8,TWR8)
      tread(1) = tread(1)+(TCR8-TCR7)
      tread(2) = tread(2)+(TWR8-TWR7)

      ! *************** EVALUATION OF THE (WA|XY) INTEGRALS ***********

      call CWTIME(TCINT1,TWINT1)

      DoTraInt = (JRED == JRED2) .and. (iBatch == nBatch)

      call CHO_eval_waxy(irc,Scr,Laq(3),Lxy,W_PWXY,nAorb,JSYM,JNUM,DoTraInt,CMO)

      call CWTIME(TCINT2,TWINT2)
      tintg(1) = tintg(1)+(TCINT2-TCINT1)
      tintg(2) = tintg(2)+(TWINT2-TWINT1)

      if (irc /= 0) then
        rc = irc
        return
      end if

      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------
      call Deallocate_DT(Lxy)
      call Deallocate_DT(Laq(3))

    end do  ! end batch loop

    if (JSYM == 1) then
      ! backtransform fock matrix in full storage
      add = .true.
      call swap_rs2full(irc,iLoc,nRS,nDen,JSYM,FLT,Frs,add)
    end if

    ! free memory
    call mma_deallocate(Lrs)

    if (JSYM == 1) then
      call mma_deallocate(Drs)
      call mma_deallocate(Frs)
    end if

  end do ! loop over red sets

  call Deallocate_DT(Scr)

end do ! loop over JSYM

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky CASSCF timing from '//SECNAM
  write(u6,CFmt) '----------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'Fock matrix construction        CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ/TRANSFORM VECTORS                    ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'COULOMB                                   ',tcoul(1),tcoul(2)
  write(u6,'(2x,A26,2f10.2)') 'EXCHANGE                                  ',texch(1),texch(2)
  write(u6,'(2x,A26,2f10.2)') '(WA|XY) INTEGRALS                         ',tintg(1),tintg(2)
  write(u6,*)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
if (Debug) then ! to avoid double printing in CASSCF-debug

  write(u6,'(6X,A)') 'TEST PRINT FROM '//SECNAM
  write(u6,'(6X,A)')
  write(u6,'(6X,A)') '***** INACTIVE FOCK MATRIX ***** '
  do ISYM=1,NSYM
    if (NBAS(ISYM) > 0) then
      write(u6,'(6X,A)')
      write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
      call TRIPRT('','',FLT(1)%SB(ISYM)%A1,NBAS(ISYM))
    end if
  end do
  if (DoActive) then
    write(u6,'(6X,A)')
    write(u6,'(6X,A)') '***** ACTIVE FOCK MATRIX ***** '
    do ISYM=1,NSYM
      if (NBAS(ISYM) > 0) then
        write(u6,'(6X,A)')
        write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
        call TRIPRT('','',FLT(2)%SB(ISYM)%A1,NBAS(ISYM))
      end if
    end do
  end if

end if
#endif

rc = 0

return

end subroutine CHO_FMCSCF
