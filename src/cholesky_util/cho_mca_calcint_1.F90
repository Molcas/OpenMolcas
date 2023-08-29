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

subroutine CHO_MCA_CALCINT_1(ISHLAB)
!
! Purpose: calculate qualified integral columns from
!          shell pair distribution (**|ISHLA ISHLB).
!
! Version 1: store full shell quadruple.

use Index_Functions, only: nTri_Elem
use Cholesky, only: iiBstR, iiBstRSh, IndRed, INF_IN2, INF_INT, iOff_Col, iOffQ, IPRINT, iQuAB, iSP2F, nBstSh, LuPri, LuSel, &
                    nnBstR, nnBstRSh, nnShl, nQual, nSym, TINTEG
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ISHLAB
integer(kind=iwp) :: IAB, IADR, ICD, IOPT, ISHLA, ISHLB, ISHLC, ISHLCD, ISHLD, ISYM, JAB, JCD, JCD0, JCDS, KAB, KOFF, KOFF1, &
                     KOFF2, L4SH, L4SHMX, LCOL, LINT, LTOT, MAXCD, NAB(8), NUMAB, NUMCD
real(kind=wp) :: C1, C2, PCT, W1, W2, XSKIP, XXSHL
logical(kind=iwp) :: DOINTS
real(kind=wp), allocatable :: Int4Sh(:), IntCol(:)
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_MCA_CALCINT_1'

#ifdef _DEBUGPRINT_
call mma_maxDBLE(LLEAK)
MEM_START = LLEAK
#endif

! Initializations.
! ----------------

call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)

NAB(1) = 0 ! dummy initialization
NAB(1:NSYM) = NQUAL(1:NSYM)-IOFFQ(1:NSYM)

IOFF_COL(1) = 0
LCOL = NNBSTR(1,2)*NAB(1)
do ISYM=2,NSYM
  IOFF_COL(ISYM) = LCOL
  LCOL = LCOL+NNBSTR(ISYM,2)*NAB(ISYM)
end do

if (ISHLA == ISHLB) then
  NUMAB = nTri_Elem(NBSTSH(ISHLA))
else
  NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
end if
MAXCD = 0
do ISHLCD=1,NNSHL
  ISYM = 1
  DOINTS = (NAB(ISYM) > 0) .and. (NNBSTRSH(ISYM,ISHLCD,2) > 0)
  do while ((ISYM < NSYM) .and. (.not. DOINTS))
    ISYM = ISYM+1
    DOINTS = (NAB(ISYM) > 0) .and. (NNBSTRSH(ISYM,ISHLCD,2) > 0)
  end do
  if (DOINTS) then
    call CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.true.)
    if (ISHLC == ISHLD) then
      NUMCD = nTri_Elem(NBSTSH(ISHLC))
    else
      NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
    end if
    MAXCD = max(MAXCD,NUMCD)
  end if
end do
L4SHMX = MAXCD*NUMAB

XXSHL = real(NNSHL,kind=wp)
XSKIP = Zero

if (IPRINT >= INF_INT) write(LUPRI,*)

! Allocate memory and initialize:
! qualified columns in reduced set,
! max. shell quadruple.
! ---------------------------------

call mma_allocate(Int4Sh,L4SHMX,Label='Int4Sh')
call mma_allocate(IntCol,LCOL,Label='IntCol')
IntCol(:) = Zero

! Set memory used by seward.
! --------------------------

call mma_maxDBLE(LINT)
call XSETMEM_INTS(LINT)

! Loop over shell quadruples.
! ---------------------------

do ISHLCD=1,NNSHL

  ! Set left shell pair index.
  ! --------------------------

  call CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.true.)
  if (ISHLC == ISHLD) then
    NUMCD = nTri_Elem(NBSTSH(ISHLC))
  else
    NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
  end if

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

    if (IPRINT >= INF_INT) write(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)') 'Invoking Seward for shell quadruple (',ISHLC,ISHLD,'|',ISHLA, &
                                                                    ISHLB,')'

    ! Calculate integrals.
    ! --------------------

    call CWTIME(C1,W1)
    L4SH = NUMCD*NUMAB
    Int4Sh(1:L4SH) = Zero
    call CHO_MCA_INT_1(ISHLCD,ISHLAB,Int4SH,L4SH,LOCDBG .or. (IPRINT >= 100))
    call CWTIME(C2,W2)
    TINTEG(1,1) = TINTEG(1,1)+C2-C1
    TINTEG(2,1) = TINTEG(2,1)+W2-W1

    ! Extract columns in reduced set.
    ! IAB: index AB within full shell pair.
    ! JAB: index AB within current reduced set.
    ! KAB: index AB within qualifieds.
    ! -----------------------------------------

    do ISYM=1,NSYM
      do KAB=1,NAB(ISYM)

        JAB = IQUAB(IOFFQ(ISYM)+KAB,ISYM)
        IAB = INDRED(INDRED(JAB,2),1)

        do JCD0=1,NNBSTRSH(ISYM,ISHLCD,2)

          JCDS = IIBSTRSH(ISYM,ISHLCD,2)+JCD0
          JCD = IIBSTR(ISYM,2)+JCDS
          ICD = INDRED(INDRED(JCD,2),1)

          KOFF1 = IOFF_COL(ISYM)+NNBSTR(ISYM,2)*(KAB-1)+JCDS
          KOFF2 = NUMCD*(IAB-1)+ICD

          IntCol(KOFF1) = Int4Sh(KOFF2)

        end do

      end do
    end do

  else

    ! Update skip counter.
    ! --------------------

    XSKIP = XSKIP+One

    ! Print message.
    ! --------------

    if (IPRINT >= INF_INT) write(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)') 'NOTICE: skipping shell quadruple    (',ISHLC,ISHLD,'|',ISHLA, &
                                                                    ISHLB,')'

  end if

end do

! Write the columns to disk.
! --------------------------

call CWTIME(C1,W1)
do ISYM=1,NSYM
  LTOT = NNBSTR(ISYM,2)*NAB(ISYM)
  if (LTOT > 0) then
    IOPT = 1
    KOFF = 1+IOFF_COL(ISYM)
    IADR = NNBSTR(ISYM,2)*IOFFQ(ISYM)
    call DDAFILE(LUSEL(ISYM),IOPT,IntCol(KOFF),LTOT,IADR)
  end if
end do
call CWTIME(C2,W2)
TINTEG(1,2) = TINTEG(1,2)+C2-C1
TINTEG(2,2) = TINTEG(2,2)+W2-W1

! Free memory: both memory used by seward and used here.
! ------------------------------------------------------

call XRLSMEM_INTS()
call mma_deallocate(IntCol)
call mma_deallocate(Int4Sh)

! Print skip statistics.
! ----------------------

if (IPRINT >= INF_IN2) then
  PCT = 1.0e2_wp*XSKIP/XXSHL
  write(LUPRI,'(A,F7.2,A)') 'Skipped',PCT,'% of rows (shell pairs) in this distribution'
end if

#ifdef _DEBUGPRINT_
call mma_maxDBLE(LLEAK)
MEM_END = LLEAK
LEAK = MEM_END-MEM_START
if (LEAK /= 0) then
  write(LUPRI,'(//,A,A,I9)') SECNAM,': Memory leak:',LEAK
  call CHO_QUIT('Memory leak detected in '//SECNAM,104)
end if
#endif

end subroutine CHO_MCA_CALCINT_1
