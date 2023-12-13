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

subroutine CHO_MCA_CALCINT_3(XINT,LINT,ISHLAB)
!
! Purpose: calculate qualified integral columns from
!          shell pair distribution (**|ISHLA ISHLB).
!
! Version 3: avoid storage of full shell quadruple in interface to
!            seward; get qualified directly as in Version 2!
!            Changes from Version 2:
!            - addressing of qualified columns
!            - integrals returned in core (no I/O)

use Cholesky, only: INF_IN2, INF_INT, IPRINT, iSP2F, LuPri, NCOLAB, nnBstRSh, nnShl, nSym, TINTEG
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LINT, ISHLAB
real(kind=wp), intent(inout) :: XINT(LINT)
integer(kind=iwp) :: i, ILOC, IRC, ISHLA, ISHLB, ISHLC, ISHLCD, ISHLD, ISYM, NAB(8)
real(kind=wp) :: C1, C2, PCT, W1, W2, XSKIP, XXSHL
logical(kind=iwp) :: DOINTS
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_MCA_CALCINT_3'

! Initializations.
! ----------------

call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)

XXSHL = real(NNSHL,kind=wp)
XSKIP = Zero

if (IPRINT >= INF_INT) write(LUPRI,*)

! Set mapping from shell pair AB to qualified columns.
! ----------------------------------------------------

IRC = 0
ILOC = 2
call CHO_SETSHP2Q_2(IRC,ILOC,ISHLAB,NAB)
if (IRC /= 0) then
  write(LUPRI,*) SECNAM,': CHO_SETSHP2Q_2 returned ',IRC
  call CHO_QUIT('Error termination in '//SECNAM,IRC)
end if

! Print.
! ------

if (IPRINT >= INF_IN2) then
  NCOLAB = sum(NAB(1:NSYM))
  write(LUPRI,'(/,A,I5,1X,I5,A,I9,A)') 'Calculating shell pair (**|',ISHLA,ISHLB,'):',NCOLAB,' columns have been qualified'
  write(LUPRI,'(80A)') ('=',i=1,77)
end if

! Loop over shell quadruples.
! ---------------------------

do ISHLCD=1,NNSHL

  ! Set left shell pair index.
  ! --------------------------

  call CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.true.)

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

    call CWTIME(C1,W1)
    call CHO_MCA_INT_1(ISHLCD,ISHLAB,XINT,LINT,LOCDBG .or. (IPRINT >= 100))
    call CWTIME(C2,W2)
    TINTEG(1,1) = TINTEG(1,1)+C2-C1
    TINTEG(2,1) = TINTEG(2,1)+W2-W1

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

! Print skip statistics.
! ----------------------

if (IPRINT >= INF_IN2) then
  PCT = 1.0e2_wp*XSKIP/XXSHL
  write(LUPRI,'(A,F7.2,A)') 'Skipped',PCT,'% of rows (shell pairs) in this distribution'
  call XFLUSH(LUPRI)
end if

end subroutine CHO_MCA_CALCINT_3
