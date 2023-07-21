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

subroutine CHO_MCA_DBGINT_S(ISHLQ,NSHLQ,PRTLAB)
!
! Purpose: regenerate and check a specified set of integrals.
!
! NOTE: this is *not* meant for production calculations, only for
!       debugging, as:
!       1) full Cholesky vectors are read
!       2) calculations are performed in full (no use of red. sets
!          apart from first)

use ChoArr, only: nBstSh
use stdalloc

implicit real*8(a-h,o-z)
integer ISHLQ(4,NSHLQ)
logical PRTLAB
#include "cholesky.fh"
#include "choorb.fh"
character(len=16), parameter :: SECNAM = 'CHO_MCA_DBGINT_S'
integer, external :: CHO_F2SP
real*8 XLBAS(8)
character(len=8) LABEL
real*8, allocatable :: INT1(:), WRK(:)
! Statement functions
MULD2H(I,J) = ieor(I-1,J-1)+1
ITRI(I,J) = max(I,J)*(max(I,J)-3)/2+I+J

! Return if nothing specified.
! ----------------------------

if (NSHLQ < 1) return

! Force computation of full shell quadruple.
! ------------------------------------------

if (IFCSEW /= 1) then
  write(LUPRI,*) SECNAM,': WARNING: resetting IFCSEW from ',IFCSEW,' to 1.'
  write(LUPRI,*) SECNAM,': memory demands are significantly increased by this!'
  IFCSEW = 1
end if

! Initializations.
! ----------------

GLMAX = 0.0d0
GLMIN = 1.0d15
GLRMS = 0.0d0
XTCMP = 0.0d0
XPECT = 0.0d0

! Make first reduced set the current reduced set.
! -----------------------------------------------

call CHO_RSCOPY(1,2)

! Allocate memory for largest integral quadruple.
! -----------------------------------------------

LINTMX = MX2SH*MX2SH
call mma_allocate(INT1,LINTMX,Label='INT1')

! Allocate max. memory.
! ---------------------

call mma_maxDBLE(LWRK)
call mma_allocate(WRK,LWRK/2,Label='WRK')
call XSETMEM_INTS(LWRK/2)

! Print header.
! -------------

call CHO_HEAD('Integral Error Analysis','=',80,LUPRI)
write(LUPRI,'(/,A,/,A)') '    C     D     A     B   Abs. Min.    Abs. Max.      RMS', &
                         '--------------------------------------------------------------'

! Loop over specified shell quadruples.
! -------------------------------------

do I=1,NSHLQ

  ISHLC = ISHLQ(1,I)
  ISHLD = ISHLQ(2,I)
  ISHLA = ISHLQ(3,I)
  ISHLB = ISHLQ(4,I)

  if ((ISHLC > 0) .and. (ISHLD > 0) .and. (ISHLA > 0) .and. (ISHLB > 0)) then

    if (ISHLD == ISHLC) then
      NUMCD = NBSTSH(ISHLC)*(NBSTSH(ISHLC)+1)/2
    else
      NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
    end if
    if (ISHLB == ISHLA) then
      NUMAB = NBSTSH(ISHLA)*(NBSTSH(ISHLA)+1)/2
    else
      NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
    end if
    LINT1 = NUMCD*NUMAB ! actual space needed for (CD|AB)

    ! Compute expected number of comparisons.
    ! ---------------------------------------

    XPECTL = dble(LINT1)
    XPECT = XPECT+XPECTL

    ! Calculate shell quadruple (CD|AB).
    ! ----------------------------------

    ISHLCD = CHO_F2SP(ITRI(ISHLC,ISHLD))
    ISHLAB = CHO_F2SP(ITRI(ISHLA,ISHLB))
    if ((ISHLAB < 1) .or. (ISHLCD < 1)) call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)

    INT1(1:LINT1) = 0.0d0
    call CHO_MCA_INT_1(ISHLCD,ISHLAB,INT1,LINT1,.false.)

    ! Calculate integrals from Cholesky vectors.
    ! ------------------------------------------

    call CHO_DBGINT_CHO(INT1,NUMCD,NUMAB,WRK,LWRK/2,ERRMAX,ERRMIN,ERRRMS,NCMP,ISHLCD,ISHLAB)

    ! Write report.
    ! -------------

    if (NCMP < 1) then
      write(LUPRI,'(4(I5,1X),5X,A)') ISHLC,ISHLD,ISHLA,ISHLB,' !!! nothing compared !!! '
    else
      XCMP = dble(NCMP)
      XTCMP = XTCMP+XCMP
      RMS = sqrt(ERRRMS/XCMP)
      if (PRTLAB) then
        call CHO_INTCHK_ID_OF(LABEL,I,-1)
        write(LUPRI,'(4(I5,1X),1P,3(D12.4,1X),A,A,A)') ISHLC,ISHLD,ISHLA,ISHLB,ERRMIN,ERRMAX,RMS,'(',LABEL,')'
      else
        write(LUPRI,'(4(I5,1X),1P,3(D12.4,1X))') ISHLC,ISHLD,ISHLA,ISHLB,ERRMIN,ERRMAX,RMS
      end if
    end if

    if (abs(ERRMAX) > abs(GLMAX)) GLMAX = ERRMAX
    if (abs(ERRMIN) < abs(GLMIN)) GLMIN = ERRMIN
    GLRMS = GLRMS+ERRRMS

  end if

end do

! Print end of table.
! -------------------

write(LUPRI,'(A)') '--------------------------------------------------------------'
if (XTCMP < 1.0d0) then
  write(LUPRI,'(A,23X,A)') 'Total:',' !!! nothing compared !!! '
else
  GLRMS = sqrt(GLRMS/XTCMP)
  write(LUPRI,'(A,18X,1P,3(D12.4,1X))') 'Total:',GLMIN,GLMAX,GLRMS
end if
write(LUPRI,'(A)') '--------------------------------------------------------------'

! Release all memory allocated here (and release seward memory).
! --------------------------------------------------------------

call XRLSMEM_INTS()
call mma_deallocate(WRK)
call mma_deallocate(INT1)

! Check total number of comparisons.
! ----------------------------------

do ISYM=1,NSYM
  XLBAS(ISYM) = dble(NBAS(ISYM))
end do
XNINT = 0.0d0
do ISYM=1,NSYM
  XXLBST = 0.0d0
  do ISYMB=1,NSYM
    ISYMA = MULD2H(ISYMB,ISYM)
    if (ISYMA == ISYMB) then
      XXLBST = XXLBST+XLBAS(ISYMA)*(XLBAS(ISYMA)+1.0d0)/2.0d0
    else if (ISYMA > ISYMB) then
      XXLBST = XXLBST+XLBAS(ISYMA)*XLBAS(ISYMB)
    end if
  end do
  XNINT = XNINT+XXLBST*(XXLBST+1.0d0)/2.0d0
end do

if (abs(XTCMP-XPECT) > 1.0D-15) then
  write(LUPRI,'(/,A)') 'WARNING: not all integrals checked:'
else
  write(LUPRI,*)
end if
write(LUPRI,'(A,1P,D20.10)') 'Total number of integral comparisons    :',XTCMP
write(LUPRI,'(A,1P,D20.10)') 'Total number expected (full shell pairs):',XPECT
write(LUPRI,'(A,1P,D20.10)') 'Total number of unique integrals        :',XNINT

end subroutine CHO_MCA_DBGINT_S
