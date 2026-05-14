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

subroutine CHO_DIACHO(DIAG,ISYM,WRK,LWRK)
!
! Purpose: update (i.e. subtract contributions from vectors on disk)
!          of symmetry block ISYM of diagonal in red. set 1.
!          This emulates the actual procedure during decomposition.

use Cholesky, only: Cho_DecAlg, iiBstR, IndRed, InfVec, nnBstR, nSys_call, NumCho, SCDIAG
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ISYM, LWRK
real(kind=wp), intent(inout) :: Diag(*), WRK(LWRK)
integer(kind=iwp) :: IAB, IABG, ILOC, IREDC, IVEC1, JAB, JRED, JVEC, KAB, KOFFV, MUSED, NCONV, NNEG, NNEGT, NSCALL, NVRD
real(kind=wp) :: DMX, XM, XMAX, XMIN
logical(kind=iwp) :: SCDIAG_SAVE
character(len=*), parameter :: SECNAM = 'CHO_DIACHO'

! Return if nothing to do.
! ------------------------

if (NNBSTR(ISYM,1) < 1) return
if (NUMCHO(ISYM) < 1) return

! Save read-call counter.
! -----------------------

NSCALL = NSYS_CALL

! Set pointer to scratch for reduced set indices.
! -----------------------------------------------

ILOC = 3

! Set up rs1 indices at location ILOC.
! Set IREDC to identify this.
! ------------------------------------

call CHO_RSCOPY(1,ILOC)
IREDC = 1

! Start read buffer batch loop.
! -----------------------------

IVEC1 = 1
do while (IVEC1 <= NUMCHO(ISYM))

  ! Read as many vectors as possible into buffer (entire WRK).
  ! ----------------------------------------------------------

  NVRD = 0
  MUSED = 0
  call CHO_VECRD(WRK,LWRK,IVEC1,NUMCHO(ISYM),ISYM,NVRD,IREDC,MUSED)
  if (NVRD < 1) call CHO_QUIT('Insufficient scratch space for read in '//SECNAM,101)

  ! Initialize vector offset.
  ! -------------------------

  KOFFV = 0

  ! Loop over vectors in core.
  ! --------------------------

  do JVEC=1,NVRD

    ! Set index arrays for current reduced set (if not already set).
    ! --------------------------------------------------------------

    JRED = INFVEC(IVEC1+JVEC-1,2,ISYM)
    if (JRED /= IREDC) then
      if (JRED == 1) then
        call CHO_RSCOPY(1,ILOC)
      else
        call CHO_GETRED(JRED,ILOC,.false.)
        call CHO_SETREDIND(ILOC)
      end if
      IREDC = JRED
    end if

    ! Compute contributions to diagonal.
    ! Zero the diagonal element associated with this vector.
    ! ------------------------------------------------------

    do JAB=1,NNBSTR(ISYM,ILOC)
      IAB = INDRED(IIBSTR(ISYM,ILOC)+JAB,ILOC) ! address in rs1
      KAB = KOFFV+JAB ! vector address
      DIAG(IAB) = DIAG(IAB)-WRK(KAB)*WRK(KAB)
    end do
    IABG = INFVEC(IVEC1+JVEC-1,1,ISYM)
    call CHO_P_ZERODIAG_RST(DIAG,ISYM,IABG)

    ! Check diagonal.
    ! ---------------

    if (CHO_DECALG == 4) then
      SCDIAG_SAVE = SCDIAG
      SCDIAG = .false. ! do NOT screen
      DMX = One
      call CHO_CHKDIA_A4(DIAG,DMX,ISYM,NNEG,NNEGT,NCONV,XMAX,XMIN,XM)
      SCDIAG = SCDIAG_SAVE
    else
      call CHO_CHKDIA(DIAG,ISYM,XMIN,XMAX,XM,NNEGT,NNEG,NCONV)
    end if

    ! Update vector offset.
    ! ---------------------

    KOFFV = KOFFV+NNBSTR(ISYM,ILOC)

  end do

  ! Check memory.
  ! -------------

  if (KOFFV /= MUSED) call CHO_QUIT('Memory error detected in '//SECNAM,101)

  ! Update vector counter.
  ! ----------------------

  IVEC1 = IVEC1+NVRD

end do

! Restore read-call counter.
! --------------------------

NSYS_CALL = NSCALL

end subroutine CHO_DIACHO
