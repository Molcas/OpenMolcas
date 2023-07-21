!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine CHO_MCA_CALCINT_2(ISHLAB)
!
! Purpose: calculate qualified integral columns from
!          shell pair distribution (**|ISHLA ISHLB).
!
! Version 2: avoid storage of full shell quadruple in interface to
!            seward; get qualified directly!

use ChoArr, only: iSP2F, MySP
use ChoSwp, only: nnBstRSh
use Constants
use stdalloc

implicit real*8(a-h,o-z)
#include "cholesky.fh"
#include "choprint.fh"
character(len=17), parameter :: SECNAM = 'CHO_MCA_CALCINT_2'
logical DOINTS
logical, parameter :: LOCDBG = .false.
integer, parameter :: INFINT = INF_INT, INFIN2 = INF_IN2
integer NAB(8)
real*8, allocatable :: IntCol(:)

#ifdef _DEBUGPRINT_
call mma_maxDBLE(MEM_START)
#endif

! Initializations.
! ----------------

call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)

NAB(1) = NQUAL(1)-IOFFQ(1)
do ISYM=2,NSYM
  NAB(ISYM) = NQUAL(ISYM)-IOFFQ(ISYM)
end do

IOFF_COL(1) = 0
LCOL = NNBSTR(1,2)*NAB(1)
do ISYM=2,NSYM
  IOFF_COL(ISYM) = LCOL
  LCOL = LCOL+NNBSTR(ISYM,2)*NAB(ISYM)
end do

XXSHL = dble(NNSHL)
XSKIP = Zero

if (IPRINT >= INFINT) write(LUPRI,*)

! Allocate memory and initialize:
! qualified columns in reduced set,
! max. shell quadruple.
! ---------------------------------

call mma_allocate(IntCol,LCOL,Label='IntCol')
IntCol(:) = Zero

! Set mapping from shell pair AB to qualified columns.
! ----------------------------------------------------

IRC = 0
ILOC = 2
call CHO_P_SETSHP2Q(IRC,ILOC,ISHLAB,NAB)
if (IRC /= 0) then
  write(LUPRI,*) SECNAM,': CHO_SETSHP2Q returned ',IRC
  call CHO_QUIT('Error termination in '//SECNAM,IRC)
end if

! Set memory used by seward.
! --------------------------

call mma_maxDBLE(LINT)
call XSETMEM_INTS(LINT)

! Loop over shell quadruples.
! ---------------------------

do ISHLCD=1,NNSHL

  ! Set left shell pair index.
  ! --------------------------

  ISCD = MYSP(ISHLCD)
  call CHO_INVPCK(ISP2F(ISCD),ISHLC,ISHLD,.true.)

  ! Find out if this shell pair (CD) contributes to
  ! current reduced set.
  ! -----------------------------------------------

  ISYM = 1
  DOINTS = (NAB(ISYM) > 0) .and. (NNBSTRSH(ISYM,ISHLCD,2) > 0)
  do while ((ISYM < NSYM) .and. (.not. DOINTS))
    ISYM = ISYM+1
    DOINTS = (NAB(ISYM) > 0) .and. (NNBSTRSH(ISYM,ISHLCD,2) > 0)
  end do

  if (DOINTS) then

    ! Print message.
    ! --------------

    if (IPRINT >= INFINT) write(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)') 'Invoking Seward for shell quadruple (',ISHLC,ISHLD,'|',ISHLA, &
                                                                   ISHLB,')'

    ! Set mapping from shell pair CD to reduced set.
    ! ----------------------------------------------

    IRC = 0
    ILOC = 2
    call CHO_SETSHP2RS(IRC,ILOC,ISHLCD,NAB)
    if (IRC /= 0) then
      write(LUPRI,*) SECNAM,': CHO_SETSHP2RS returned ',IRC
      call CHO_QUIT('Error termination in '//SECNAM,IRC)
    end if

    ! Calculate integrals.
    ! --------------------

    call CHO_TIMER(C1,W1)
    call CHO_MCA_INT_1(ISCD,ISHLAB,IntCol,LCOL,LOCDBG .or. (IPRINT >= 100))
    call CHO_TIMER(C2,W2)
    TINTEG(1,1) = TINTEG(1,1)+C2-C1
    TINTEG(2,1) = TINTEG(2,1)+W2-W1

  else

    ! Update skip counter.
    ! --------------------

    XSKIP = XSKIP+One

    ! Print message.
    ! --------------

    if (IPRINT >= INFINT) write(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)') 'NOTICE: skipping shell quadruple    (',ISHLC,ISHLD,'|',ISHLA, &
                                                                   ISHLB,')'

  end if

end do

! Write the columns to disk.
! --------------------------

call CHO_TIMER(C1,W1)
do ISYM=1,NSYM
  LTOT = NNBSTR(ISYM,2)*NAB(ISYM)
  if (LTOT > 0) then
    IOPT = 1
    KOFF = 1+IOFF_COL(ISYM)
    IADR = NNBSTR(ISYM,2)*IOFFQ(ISYM)
    call DDAFILE(LUSEL(ISYM),IOPT,IntCol(KOFF),LTOT,IADR)
  end if
end do
call CHO_TIMER(C2,W2)
TINTEG(1,2) = TINTEG(1,2)+C2-C1
TINTEG(2,2) = TINTEG(2,2)+W2-W1

! Free memory: both memory used by seward and used here.
! ------------------------------------------------------

call XRLSMEM_INTS()
call mma_deallocate(IntCol)

! Print skip statistics.
! ----------------------

if (IPRINT >= INFIN2) then
  PCT = 1.0d2*XSKIP/XXSHL
  write(LUPRI,'(A,F7.2,A)') 'Skipped',PCT,'% of rows (shell pairs) in this distribution'
end if

#ifdef _DEBUGPRINT_
call mma_maxDBLE(MEM_END)
LEAK = MEM_END-MEM_START
if (LEAK /= 0) then
  write(LUPRI,'(//,A,A,I9)') SECNAM,': Memory leak:',LEAK
  call CHO_QUIT('Memory leak detected in '//SECNAM,104)
end if
#endif

end subroutine CHO_MCA_CALCINT_2
