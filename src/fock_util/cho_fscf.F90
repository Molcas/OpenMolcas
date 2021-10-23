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

subroutine CHO_FSCF(rc,nDen,FLT,nForb,nIorb,Porb,DLT,ExFac)
!*********************************************************************
!  Author : F. Aquilante
!
! *************** INACTIVE AO-BASIS FOCK MATRIX **********************
!
!   F(ab) = sum_J  Lab,J * V(J)  -  sum_Jk  Lka,J * Lkb,J
!
!*********************************************************************
!
!      V(J) = sum_gd  Lgd,J * Dtot(gd)
!
!      a,b,g,d:  AO-index
!      k:        MO-index   belonging to (Frozen+Inactive)
!
!*********************************************************************

use ChoArr, only: nDimRS
use ChoSwp, only: InfVec
use Data_structures, only: DSBA_Type, SBA_Type
use Data_structures, only: Allocate_SBA, Deallocate_SBA

implicit real*8(a-h,o-z)
integer rc, nDen
integer iSkip(8)
real*8 tread(2), tcoul(2), texch(2)
real*8 FactCI, FactXI, ExFac
type(DSBA_Type) Porb(nDen), DLT(nDen), FLT(nDen)
integer nForb(8,nDen), nIorb(8,nDen)
#ifdef _DEBUGPRINT_
logical Debug
#endif
logical DoRead
character*50 CFmt
character(LEN=8), parameter :: SECNAM = 'CHO_FSCF'
#include "chotime.fh"
#include "real.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"
real*8, parameter :: xone = -one
logical add
character*6 mode
real*8, allocatable :: Lrs(:,:), Drs(:), Frs(:)
real*8, allocatable :: VJ(:)
integer :: nAux(8)
type(SBA_Type) Laq(2)
!***********************************************************************
!Statement function
MulD2h(i,j) = ieor(i-1,j-1)+1
!***********************************************************************

#ifdef _DEBUGPRINT_
Debug = .false. ! to avoid double printing in CASSCF-debug
#endif

FactCI = one
FactXI = xone*ExFac

DoRead = .false.
IREDC = -1 ! unknown reduced set

if ((nDen /= 1) .and. (nDen /= 2)) then
  write(6,*) SECNAM//'Invalid parameter nDen= ',nDen
  call abend()
end if

call CWTIME(TOTCPU1,TOTWALL1) ! start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = zero ! time read/transform vectors
tcoul(:) = zero ! time for computing Coulomb
texch(:) = zero ! time for computing Exchange

! ==================================================================

iLoc = 3 ! use scratch location in reduced index arrays

! *************** BIG LOOP OVER VECTORS SYMMETRY *******************
do jSym=1,nSym

  if (NumCho(jSym) < 1) goto 1000

  ! ****************     MEMORY MANAGEMENT SECTION    *****************
  ! --------------------------------------------------------------
  ! compute memory needed to store at least 1 vector of JSYM
  ! and do all the subsequent calculations
  ! --------------------------------------------------------------
  mTvec = 0  ! mem for storing the half-transformed vec

  do l=1,nSym
    k = Muld2h(l,JSYM)
    Mmax = 0
    do jDen=1,nDen
      Mmax = max(Mmax,nForb(k,jDen)+nIorb(k,jDen))
    end do
    mTvec = mTvec+nBas(l)*Mmax
  end do

  mTvec = max(mTvec,1)

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

  JRED1 = InfVec(1,2,jSym)            ! red set of the 1st vec
  JRED2 = InfVec(NumCho(jSym),2,jSym) ! red set of the last vec

  do JRED=JRED1,JRED2

    call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

    if (nVrs == 0) goto 999 ! no vectors in that (jred,jsym)

    if (nVrs < 0) then
      write(6,*) SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
      call abend()
    end if

    call Cho_X_SetRed(irc,iLoc,JRED) ! set index arrays at iLoc
    if (irc /= 0) then
      write(6,*) SECNAM//'cho_X_setred non-zero return code.   rc= ',irc
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
      write(6,*) SECNAM//': Insufficient memory for batch'
      write(6,*) 'LWORK= ',LWORK
      write(6,*) 'min. mem. need= ',nRS+mTvec
      write(6,*) 'jsym= ',jsym
      rc = 33
      call Abend()
      nBatch = -9999 ! dummy assignment
    end if

    LREAD = nRS*nVec

    call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

    if (JSYM == 1) then
      ! Transform the density to reduced storage
      mode = 'toreds'
      add = .false.
      nMat = 1
      call Swap_rs2full(irc,iLoc,nRS,nMat,JSYM,DLT,Drs,mode,add)
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
        ! V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * DI(rs)
        !==========================================================

        call CWTIME(TCC1,TWC1)

        call mma_allocate(VJ,JNUM,Label='VJ')

        call DGEMV_('T',nRS,JNUM,ONE,Lrs,nRS,Drs,1,ZERO,VJ,1)

        ! FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
        !===============================================================

        Fact = dble(min(jVec-iVrs,1))

        call DGEMV_('N',nRS,JNUM,FactCI,Lrs,nRS,VJ,1,Fact,Frs,1)

        call CWTIME(TCC2,TWC2)
        tcoul(1) = tcoul(1)+(TCC2-TCC1)
        tcoul(2) = tcoul(2)+(TWC2-TWC1)

        call mma_deallocate(VJ)

      end if ! Coulomb contribution

      iSwap = 2 ! LpJ,b are returned
      ! *************** EXCHANGE CONTRIBUTIONS  ***********************

      do jDen=1,nDen

        nAux(:) = nForb(:,jDen)+nIorb(:,jDen)
        call Allocate_SBA(Laq(jDen),nAux,nBas,nVec,JSYM,nSym,iSwap)

        call CWTIME(TCR3,TWR3)

        kMOs = jDen ! 1--> alpha  2-->beta MOs
        nMOs = jDen

        ! Set up the skipping flags
        ! -------------------------------------------------------------
        do i=1,nSym

          k = Muld2h(i,JSYM)
          iSkip(k) = min(1,nBas(i)*(nForb(k,jDen)+nIorb(k,jDen)))

        end do
        ! -------------------------------------------------------------

        ! *********************** HALF-TRANSFORMATION  ****************

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

          iSymk = MulD2h(JSYM,iSyma)

          ! ---------------------------------------------------------------------
          ! *** Compute only the LT part of the InActive exchange matrix ********
          !
          !     FI(ab) = FI(ab) + FactXI * sum_Jk  LkJ,a * LkJ,b
          ! ---------------------------------------------------------------------
          NK = nForb(iSymk,jDen)+nIorb(iSymk,jDen)

          if (iSkip(iSymk) /= 0) then

            call DGEMM_TRI('T','N',nBas(iSyma),nBas(iSyma),NK*JNUM,FactXI,Laq(jDen)%SB(iSymk)%A3,NK*JNUM,Laq(jDen)%SB(iSymk)%A3, &
                           NK*JNUM,One,FLT(jDen)%SB(iSyma)%A1,nBas(iSyma))

          end if

          ! --------------------------------------------------------------------
        end do  !loop over MOs symmetries

        call CWTIME(TCX2,TWX2)
        texch(1) = texch(1)+(TCX2-TCX1)
        texch(2) = texch(2)+(TWX2-TWX1)

        call Deallocate_SBA(Laq(jDen))
      end do   ! loop over densities

      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------

    end do ! end batch loop

    if (JSYM == 1) then
      ! backtransform fock matrix to full storage
      mode = 'tofull'
      add = .true.
      nMat = 1
      do iDen=1,nDen
        call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,FLT(iDen),Frs,mode,add)
      end do
    end if

    ! free memory
    call mma_deallocate(Lrs)

    if (JSYM == 1) then
      call mma_deallocate(Frs)
      call mma_deallocate(Drs)
    end if

999 continue

  end do ! loop over red sets

1000 continue

end do ! loop over JSYM

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(6,*)
  write(6,CFmt) 'Cholesky SCF timing from '//SECNAM
  write(6,CFmt) '------------------------------------'
  write(6,*)
  write(6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(6,CFmt) 'Fock matrix construction        CPU       WALL   '
  write(6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(6,'(2x,A26,2f10.2)') 'READ/TRANSFORM VECTORS                    ',tread(1),tread(2)
  write(6,'(2x,A26,2f10.2)') 'COULOMB                                   ',tcoul(1),tcoul(2)
  write(6,'(2x,A26,2f10.2)') 'EXCHANGE                                  ',texch(1),texch(2)
  write(6,*)
  write(6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(6,*)

end if

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
if (Debug) then ! to avoid double printing in SCF-debug

  write(6,'(6X,A)') 'TEST PRINT FROM '//SECNAM
  write(6,'(6X,A)')
  write(6,'(6X,A)') '***** FOCK MATRIX AO-BASIS ***** '
  do jDen=1,nDen
    if (nDen == 2) then
      if (jden == 1) write(6,'(6X,A)') '******** ALPHA SPIN ******** '
      if (jden == 2) write(6,'(6X,A)') '******** BETA SPIN ********* '
    end if
    do ISYM=1,NSYM
      if (NBAS(ISYM) > 0) then
        write(6,'(6X,A)')
        write(6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
        call TRIPRT('','',FLT(jDen)%SB(ISYM)%A1,NBAS(ISYM))
      end if
    end do
  end do

end if
#endif

rc = 0

return

end subroutine CHO_FSCF
