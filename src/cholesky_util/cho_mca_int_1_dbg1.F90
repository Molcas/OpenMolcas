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

subroutine CHO_MCA_INT_1_DBG1(DIAG,IRED)
!
! Purpose: test diagonal, reduced set IRED. Note that the
!          diagonal *must* be the original diagonal stored
!          in reduced set 1.

use Index_Functions, only: nTri_Elem
use Cholesky, only: IFCSew, iiBstR, iiBstRSh, IndRed, IndRSh, iSP2F, LuPri, Mx2Sh, nBstSh, nnBstRSh, nnBstrT, nnShl, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: DIAG(*)
integer(kind=iwp), intent(in) :: IRED
integer(kind=iwp) :: IAB, IERR, ISHLA, ISHLAB, ISHLB, ISYM, JAB, JAB1, JAB2, JSHLAB, KAB, KABAB, LINT, LINT1, LSEW, NERR, NTST, &
                     NUMAB
real(kind=wp) :: DIFF
real(kind=wp), allocatable :: xINT(:)
logical(kind=iwp), parameter :: PRTINT = .false.
character(len=*), parameter :: SECNAM = 'CHO_MCA_INT_1_DBG1'

write(LUPRI,*)
write(LUPRI,*)
write(LUPRI,*) SECNAM,': testing diagonal, reduced set ',IRED
write(LUPRI,*)

! Force computation of full shell quadruple.
! ------------------------------------------

if (IFCSEW /= 1) then
  write(LUPRI,*) SECNAM,': WARNING: resetting IFCSEW from ',IFCSEW,' to 1.'
  IFCSEW = 1
end if

LINT1 = MX2SH*MX2SH
call mma_allocate(xINT,LINT1,Label='INT')
call mma_maxDBLE(LSEW)
call XSETMEM_INTS(LSEW)

NERR = 0
NTST = 0
do ISHLAB=1,NNSHL

  ! Allocate memory for shell quadruple (AB|AB).
  ! --------------------------------------------

  call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
  if (ISHLB == ISHLA) then
    NUMAB = nTri_Elem(NBSTSH(ISHLA))
  else
    NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
  end if
  LINT = NUMAB*NUMAB

  ! Calculate integrals.
  ! --------------------

  xINT(1:LINT) = Zero
  call CHO_MCA_INT_1(ISHLAB,ISHLAB,xINT,LINT,PRTINT)

  ! Look up all diagonal elements in DIAG and compare to
  ! values just calculated.
  ! ----------------------------------------------------

  IERR = 0
  if (IRED == 1) then

    do ISYM=1,NSYM

      JAB1 = IIBSTR(ISYM,1)+IIBSTRSH(ISYM,ISHLAB,1)+1
      JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,1)-1

      do JAB=JAB1,JAB2   ! loop over elements in diagonal

        if ((JAB < 1) .or. (JAB > NNBSTRT(1))) then
          write(LUPRI,*) SECNAM,': JAB = ',JAB
          write(LUPRI,*) SECNAM,': should be between 1 and ',NNBSTRT(1)
          call CHO_QUIT(SECNAM//': index error (IRED=1)',103)
        end if

        JSHLAB = INDRSH(JAB)
        if (JSHLAB /= ISP2F(ISHLAB)) then
          write(LUPRI,*) SECNAM,': test is meaningless!'
          write(LUPRI,*) SECNAM,': JSHLAB must equal ','ISP2F(ISHLAB)'
          write(LUPRI,*) SECNAM,': JSHLAB,ISP2F(ISHLAB): ',JSHLAB,ISP2F(ISHLAB)
          call CHO_QUIT(SECNAM//': shell bug (IRED=1)',103)
        end if

        IAB = INDRED(JAB,1)
        if ((IAB < 1) .or. (IAB > NUMAB)) then
          write(LUPRI,*) SECNAM,': IAB = ',IAB
          write(LUPRI,*) SECNAM,': should be between 1 and ',NUMAB
          call CHO_QUIT(SECNAM//': index error (IRED=1)',103)
        end if

        KABAB = NUMAB*(IAB-1)+IAB
        DIFF = DIAG(JAB)-xINT(KABAB)
        if (abs(DIFF) > 1.0e-14_wp) then
          write(LUPRI,*) SECNAM,': ISHLA,ISHLB,JAB,IAB,DIFF: ',ISHLA,ISHLB,JAB,IAB,DIFF
          IERR = IERR+1
        end if

        NTST = NTST+1

      end do

    end do

    write(LUPRI,*) SECNAM,': ISHLA,ISHLB,#errors: ',ISHLA,ISHLB,IERR

    NERR = NERR+IERR

  else if ((IRED == 2) .or. (IRED == 3)) then

    do ISYM=1,NSYM

      JAB1 = IIBSTR(ISYM,IRED)+IIBSTRSH(ISYM,ISHLAB,IRED)+1
      JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,IRED)-1

      do JAB=JAB1,JAB2   ! loop over elements in diagonal

        if ((JAB < 1) .or. (JAB > NNBSTRT(IRED))) then
          write(LUPRI,*) SECNAM,': JAB = ',JAB
          write(LUPRI,*) SECNAM,': should be between 1 and ',NNBSTRT(IRED)
          call CHO_QUIT(SECNAM//': index error (IRED>1)',103)
        end if

        JSHLAB = INDRSH(INDRED(JAB,IRED))
        if (JSHLAB /= ISP2F(ISHLAB)) then
          write(LUPRI,*) SECNAM,': test is meaningless!'
          write(LUPRI,*) SECNAM,': JSHLAB must equal ISP2F(ISHLAB)'
          write(LUPRI,*) SECNAM,': JSHLAB,ISP2F(ISHLAB): ',JSHLAB,ISP2F(ISHLAB)
          call CHO_QUIT(SECNAM//': shell bug (IRED>1)',103)
        end if

        KAB = INDRED(JAB,IRED)  ! index in red. set 1
        if ((KAB < 1) .or. (KAB > NNBSTRT(1))) then
          write(LUPRI,*) SECNAM,': KAB = ',KAB
          write(LUPRI,*) SECNAM,': should be between 1 and ',NNBSTRT(1)
          call CHO_QUIT(SECNAM//': index error (IRED>1)',103)
        end if

        IAB = INDRED(KAB,1)
        if ((IAB < 1) .or. (IAB > NUMAB)) then
          write(LUPRI,*) SECNAM,': IAB = ',IAB
          write(LUPRI,*) SECNAM,': should be between 1 and ',NUMAB
          call CHO_QUIT(SECNAM//': index error (IRED>1)',103)
        end if

        KABAB = NUMAB*(IAB-1)+IAB
        DIFF = DIAG(KAB)-xINT(KABAB)
        if (abs(DIFF) > 1.0e-14_wp) then
          write(LUPRI,*) SECNAM,': ISHLA,ISHLB,JAB,IAB,DIFF: ',ISHLA,ISHLB,JAB,IAB,DIFF
          IERR = IERR+1
        end if

        NTST = NTST+1

      end do

    end do

    write(LUPRI,*) SECNAM,': ISHLA,ISHLB,#errors: ',ISHLA,ISHLB,IERR

    NERR = NERR+IERR

  else

    call CHO_QUIT(SECNAM//': IRED out of bounds!',104)

  end if

end do

call XRLSMEM_INTS()
call mma_deallocate(xINT)

write(LUPRI,*) '***END OF ',SECNAM,': #tests: ',NTST,' #errors: ',NERR

end subroutine CHO_MCA_INT_1_DBG1
