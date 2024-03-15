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

subroutine CHO_MCA_DBGINT_A()
!
! Purpose: regenerate and check all integrals (or the number
!          of columns specified in input).
!
! NOTE: this is *not* meant for production calculations, only for
!       debugging, as:
!       1) full Cholesky vectors are read
!       2) calculations are performed in full (no use of red. sets
!          apart from first)
!       3) full integral symmetry not used
!          (only partial particle permutation symmetry)

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Cholesky, only: IFCSEW, iSP2F, LuPri, MX2SH, NBAS, nBstSh, nCol_chk, nnShl, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ISAB1, ISAB2, ISCD1, ISCD2, ISHLA, ISHLAB, ISHLB, ISHLC, ISHLCD, ISHLD, ISYM, ISYMA, ISYMB, LINT1, LINTMX, &
                     LWRK, NCMP, NUMAB, NUMCD
real(kind=wp) :: ERRMAX, ERRMIN, ERRRMS, GLMAX, GLMIN, GLRMS, RMS, XLBAS(8), XNINT, XPECT, XPECTL, XTCMP, XXLBST
real(kind=wp), allocatable :: Int1(:), Wrk(:)
character(len=*), parameter :: SECNAM = 'CHO_MCA_DBGINT_A'

! Force computation of full shell quadruple.
! ------------------------------------------

if (IFCSEW /= 1) then
  write(LUPRI,*) SECNAM,': WARNING: resetting IFCSEW from ',IFCSEW,' to 1.'
  write(LUPRI,*) SECNAM,': memory demands are significantly increased by this!'
  IFCSEW = 1
end if

! Initializations.
! ----------------

GLMAX = Zero
GLMIN = 1.0e15_wp
GLRMS = Zero
XTCMP = Zero
XPECT = Zero

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

! Loop over shell quadruples.
! ---------------------------

ISAB1 = 1
if (NCOL_CHK > 0) then
  ISAB2 = min(NCOL_CHK,NNSHL)
else
  ISAB2 = NNSHL
end if

do ISHLAB=ISAB1,ISAB2

  call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
  if (ISHLB == ISHLA) then
    NUMAB = nTri_Elem(NBSTSH(ISHLA))
  else
    NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
  end if

  if (NCOL_CHK > 0) then
    ISCD1 = 1
    ISCD2 = NNSHL
  else
    ISCD1 = ISHLAB
    ISCD2 = NNSHL
  end if

  do ISHLCD=ISCD1,ISCD2

    call CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.true.)
    if (ISHLD == ISHLC) then
      NUMCD = nTri_Elem(NBSTSH(ISHLC))
    else
      NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
    end if
    LINT1 = NUMCD*NUMAB ! actual space needed for (CD|AB)

    ! Compute expected number of comparisons.
    ! ---------------------------------------

    XPECTL = real(LINT1,kind=wp)
    XPECT = XPECT+XPECTL

    ! Calculate shell quadruple (CD|AB).
    ! ----------------------------------

    INT1(1:LINT1) = Zero
    call CHO_MCA_INT_1(ISHLCD,ISHLAB,INT1,LINT1,.false.)

    ! Calculate integrals from Cholesky vectors.
    ! ------------------------------------------

    call CHO_DBGINT_CHO(INT1,NUMCD,NUMAB,WRK,LWRK/2,ERRMAX,ERRMIN,ERRRMS,NCMP,ISHLCD,ISHLAB)

    ! Write report.
    ! -------------

    if (NCMP < 1) then
      write(LUPRI,'(4(I5,1X),5X,A)') ISHLC,ISHLD,ISHLA,ISHLB,' !!! nothing compared !!! '
    else
      RMS = sqrt(ERRRMS/real(NCMP,kind=wp))
      write(LUPRI,'(4(I5,1X),3(ES12.4,1X))') ISHLC,ISHLD,ISHLA,ISHLB,ERRMIN,ERRMAX,RMS
    end if

    if (abs(ERRMAX) > abs(GLMAX)) GLMAX = ERRMAX
    if (abs(ERRMIN) < abs(GLMIN)) GLMIN = ERRMIN
    GLRMS = GLRMS+ERRRMS
    if (NCMP > 0) XTCMP = XTCMP+real(NCMP,kind=wp)

  end do

end do

! Print end of table.
! -------------------

write(LUPRI,'(A)') '--------------------------------------------------------------'
if (XTCMP < One) then
  write(LUPRI,'(A,23X,A)') 'Total:',' !!! nothing compared !!! '
else
  GLRMS = sqrt(GLRMS/XTCMP)
  write(LUPRI,'(A,18X,3(ES12.4,1X))') 'Total:',GLMIN,GLMAX,GLRMS
end if
write(LUPRI,'(A)') '--------------------------------------------------------------'

! Release all memory allocated here (and release seward memory).
! --------------------------------------------------------------

call XRLSMEM_INTS()
call mma_deallocate(WRK)
call mma_deallocate(INT1)

! Check total number of comparisons.
! ----------------------------------

XLBAS(1:NSYM) = real(NBAS(1:NSYM),kind=wp)
XNINT = Zero
do ISYM=1,NSYM
  XXLBST = Zero
  do ISYMB=1,NSYM
    ISYMA = MUL(ISYMB,ISYM)
    if (ISYMA == ISYMB) then
      XXLBST = XXLBST+XLBAS(ISYMA)*(XLBAS(ISYMA)+One)*Half
    else if (ISYMA > ISYMB) then
      XXLBST = XXLBST+XLBAS(ISYMA)*XLBAS(ISYMB)
    end if
  end do
  XNINT = XNINT+XXLBST*(XXLBST+One)*Half
end do

if (abs(XTCMP-XPECT) > 1.0e-15_wp) then
  write(LUPRI,'(/,A)') 'WARNING: not all integrals checked:'
else
  write(LUPRI,*)
end if
write(LUPRI,'(A,ES20.10)') 'Total number of integral comparisons    :',XTCMP
write(LUPRI,'(A,ES20.10)') 'Total number expected (full shell pairs):',XPECT
write(LUPRI,'(A,ES20.10)') 'Total number of unique integrals        :',XNINT

end subroutine CHO_MCA_DBGINT_A
