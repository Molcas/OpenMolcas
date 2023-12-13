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

subroutine CHO_REOVC1(IRS2F,N,LRDIM,WRK,LWRK)
!
! Purpose: reorder Cholesky vectors on disk to full storage.

use Symmetry_Info, only: Mul
use Cholesky, only: iiBstR, LuPri, NABPK, NNBST, nnBstR, nSys_call, nSym, NumCho
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, LRDIM, IRS2F(N,LRDIM), LWRK
real(kind=wp), intent(inout) :: WRK(LWRK)
integer(kind=iwp) :: I, IAB, IBATCH, ICOUNT, IOFF(8,8), IRS, ISYM, ISYMA, ISYMB, IVEC, IVEC1, KCHO1, KCHO2, KOFF, KOFF1, KREAD, &
                     LOFF, LREAD, MINMEM, NBATCH, NSCALL, NUMV, NVEC
character(len=*), parameter :: SECNAM = 'CHO_REOVC1'

if (N < 3) call CHO_QUIT('Dimension error in '//SECNAM,104)

! Save read-call counter.
! -----------------------

NSCALL = NSYS_CALL

! Make rs1 the "current" reduced set (for reading).
! -------------------------------------------------

call CHO_RSCOPY(1,2)

! Loop over Cholesky vector symmetries.
! -------------------------------------

do ISYM=1,NSYM
  if (NUMCHO(ISYM) > 0) then

    ! Open files.
    ! -----------

    call CHO_OPFVEC(ISYM,1)

    ! Set up vector batch.
    ! --------------------

    MINMEM = NNBSTR(ISYM,2)+NNBST(ISYM)
    if (MINMEM < 1) then
      write(LUPRI,*) SECNAM,': MINMEM = ',MINMEM
      call CHO_QUIT('NNBST error in '//SECNAM,104)
      NVEC = 0
    else
      NVEC = min(LWRK/MINMEM,NUMCHO(ISYM))
    end if

    if (NVEC < 1) then
      write(LUPRI,*) SECNAM,': NVEC   = ',NVEC
      write(LUPRI,*) SECNAM,': LWRK   = ',LWRK
      write(LUPRI,*) SECNAM,': MINMEM = ',MINMEM
      write(LUPRI,*) SECNAM,': NUMCHO = ',NUMCHO(ISYM)
      write(LUPRI,*) SECNAM,': ISYM   = ',ISYM
      call CHO_QUIT('Batch error in '//SECNAM,101)
      NBATCH = 0
    else
      NBATCH = (NUMCHO(ISYM)-1)/NVEC+1
    end if

    ! Start batch loop over vectors.
    ! ------------------------------

    do IBATCH=1,NBATCH

      if (IBATCH == NBATCH) then
        NUMV = NUMCHO(ISYM)-NVEC*(NBATCH-1)
      else
        NUMV = NVEC
      end if
      IVEC1 = NVEC*(IBATCH-1)+1

      ! Read batch of reduced vectors.
      ! ------------------------------

      KCHO1 = 1
      KREAD = KCHO1+NNBSTR(ISYM,2)*NUMV
      LREAD = LWRK-KREAD+1
      call CHO_GETVEC(WRK(KCHO1),NNBSTR(ISYM,2),NUMV,IVEC1,ISYM,WRK(KREAD),LREAD)

      ! Reorder.
      ! --------

      KCHO2 = KREAD
      IOFF(:,:) = 0
      ICOUNT = KCHO2-1
      do ISYMB=1,NSYM
        ISYMA = MUL(ISYMB,ISYM)
        if (ISYMA >= ISYMB) then
          IOFF(ISYMA,ISYMB) = ICOUNT
          IOFF(ISYMB,ISYMA) = ICOUNT
          ICOUNT = ICOUNT+NABPK(ISYMA,ISYMB)*NUMV
        end if
      end do

      WRK(KCHO2:KCHO2+NNBST(ISYM)*NUMV-1) = Zero
      do IVEC=1,NUMV
        KOFF1 = KCHO1+NNBSTR(ISYM,2)*(IVEC-1)-1
        do IRS=1,NNBSTR(ISYM,2)
          I = IIBSTR(ISYM,2)+IRS
          ISYMA = IRS2F(1,I)
          ISYMB = IRS2F(2,I)
          IAB = IRS2F(3,I)
          KOFF = KOFF1+IRS
          LOFF = IOFF(ISYMA,ISYMB)+NABPK(ISYMA,ISYMB)*(IVEC-1)+IAB
          WRK(LOFF) = WRK(KOFF)
        end do
      end do

      ! Write full vectors to disk.
      ! ---------------------------

      do ISYMB=1,NSYM
        ISYMA = MUL(ISYMB,ISYM)
        if (ISYMA >= ISYMB) then
          KOFF = IOFF(ISYMA,ISYMB)+1
          call CHO_WRFVEC(WRK(KOFF),ISYMA,ISYMB,IVEC1,NUMV)
        end if
      end do

    end do

    ! Close files.
    ! ------------

    call CHO_OPFVEC(ISYM,2)

  end if
end do

! Restore read-call counter.
! --------------------------

NSYS_CALL = NSCALL

end subroutine CHO_REOVC1
