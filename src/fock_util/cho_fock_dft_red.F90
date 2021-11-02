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

subroutine CHO_FOCK_DFT_RED(irc,DLT,FLT)
!********************************************************
!
! Author:  F. Aquilante
!
! Coulomb fock matrix  (level 2 BLAS)
!
! --- F(ab) = 2 * sum_J  Lab,J * sum_gd  D(gd) * Lgd,J
!
!********************************************************

use ChoArr, only: nDimRS
use ChoSwp, only: InfVec
use Data_Structures, only: DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: irc
type(DSBA_Type), intent(in) :: DLT
type(DSBA_Type), intent(inout) :: FLT(1)
#include "chotime.fh"
#include "cholesky.fh"
#include "choorb.fh"
integer(kind=iwp) :: i, iBatch, iLoc, IVEC2, iVrs, JNUM, JRED, JRED1, JRED2, JSYM, JVEC, LREAD, LWork, MUSED, nBatch, nDen, nRS, &
                     NUMV, nVec, nVrs
real(kind=wp) :: FactC, TCC1, TCC2, tcoul(2), TCR1, TCR2, TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), TWC1, &
                 TWC2, TWR1, TWR2, xfac
logical(kind=iwp) :: add
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: Debug
#endif
character(len=50) :: CFmt
real(kind=wp), allocatable :: Drs(:), Frs(:), Lrs(:,:), VJ(:)
character(len=*), parameter :: SECNAM = 'CHO_FOCK_DFT_RED'

#ifdef _DEBUGPRINT_
Debug = .true.
#endif

FactC = one

! For Coulomb only, the vectors symmetry is restricted to 1
JSYM = 1
if (NumCho(JSYM) < 1) return

call CWTIME(TOTCPU1,TOTWALL1) ! start clock for total time

do i=1,2          ! 1 --> CPU   2 --> Wall
  tread(i) = zero ! time for reading the vectors
  tcoul(i) = zero ! time for computing Coulomb
end do

iLoc = 3 ! use scratch location in reduced index arrays

JRED1 = InfVec(1,2,jSym)            ! red set of the 1st vec
JRED2 = InfVec(NumCho(jSym),2,jSym) ! red set of the last vec

do JRED=JRED1,JRED2

  ! Memory management section -----------------------------
  call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

  if (nVrs == 0) cycle

  if (nVrs < 0) then
    write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs < 0. STOP!!'
    call abend()
  end if

  call Cho_X_SetRed(irc,iLoc,JRED) ! set index arrays at iLoc
  if (irc /= 0) then
    write(u6,*) SECNAM//'cho_X_setred non-zero return code. rc= ',irc
    call abend()
  end if

  nRS = nDimRS(JSYM,JRED)

  call mma_allocate(Drs,nRS,Label='Drs')
  call mma_allocate(Frs,nRS,Label='Frs')
  Drs(:) = Zero
  Frs(:) = Zero

  call mma_maxDBLE(LWork)

  nVec = min(LWORK/(nRS+1),nVrs)

  if (nVec < 1) then
    write(u6,*) SECNAM//': Insufficient memory for batch'
    write(u6,*) 'LWORK= ',LWORK
    write(u6,*) 'min. mem. need= ',nRS+1
    irc = 33
    call Abend()
    nBatch = -9999 ! dummy assignment
  end if

  LREAD = nRS*nVec

  call mma_allocate(Lrs,nRS,nVec,Label='Lrs')
  call mma_allocate(VJ,nVec,Label='VJ')

  ! Transform the density to reduced storage
  add = .false.
  nDen = 1
  call swap_full2rs(irc,iLoc,nRS,nDen,JSYM,[DLT],Drs,add)

  ! BATCH over the vectors in JSYM=1 ----------------------------

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

    call CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,NUMV,JRED,MUSED)

    if ((NUMV <= 0) .or. (NUMV /= JNUM)) then
      irc = 77
      return
    end if

    call CWTIME(TCR2,TWR2)
    tread(1) = tread(1)+(TCR2-TCR1)
    tread(2) = tread(2)+(TWR2-TWR1)

    ! ************ BEGIN COULOMB CONTRIBUTION  ****************
    !
    !-Computing the intermediate vector V(J)
    !
    ! Contraction with the density matrix
    ! -----------------------------------
    ! V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * D(rs)
    !==========================================================

    call CWTIME(TCC1,TWC1)

    call DGEMV_('T',nRS,JNUM,ONE,Lrs,nRS,Drs,1,ZERO,VJ,1)

    ! Frs{#J} <- Frs{#J} + sum_J L(rs,{#J})*V{#J}
    !======================================================

    xfac = real(min(jVec-iVrs,1),kind=wp)

    call DGEMV_('N',nRS,JNUM,FactC,Lrs,nRS,VJ,1,xfac,Frs,1)

    call CWTIME(TCC2,TWC2)
    tcoul(1) = tcoul(1)+(TCC2-TCC1)
    tcoul(2) = tcoul(2)+(TWC2-TWC1)

  end do  !end batch loop

  if (nVrs > 0) then
    ! backtransform fock matrix in full storage
    add = JRED > JRED1
    call swap_rs2full(irc,iLoc,nRS,nDen,JSYM,FLT,Frs,add)
  end if

  ! free memory
  call mma_deallocate(VJ)
  call mma_deallocate(Lrs)
  call mma_deallocate(Frs)
  call mma_deallocate(Drs)

end do ! loop over red sets

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky SCF timing from '//SECNAM
  write(u6,CFmt) '-----------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'Fock matrix construction        CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ VECTORS                              ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'COULOMB                                   ',tcoul(1),tcoul(2)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
if (Debug) then ! to avoid double printing in SCF-debug

  write(u6,'(6X,A)') 'TEST PRINT FROM '//SECNAM
  write(u6,'(6X,A)')
  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NB > 0) then
      write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
      call TRIPRT('Coulomb Fmat',' ',FLT(1)%SB(ISYM)%A1,NB)
    end if
  end do

end if
#endif

irc = 0

return

end subroutine CHO_FOCK_DFT_RED
