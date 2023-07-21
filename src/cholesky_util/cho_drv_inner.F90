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
! Copyright (C) 2003-2010, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine CHO_DRV_Inner(IRETURN)
!
! Thomas Bondo Pedersen, 2003-2010.
!
! Purpose: driver for the Cholesky decomposition of two-electron
!          integrals. On entry, the integral program must have been
!          initialized and the relevant index info (#irreps, basis
!          functions, shells, etc.) must have been set up.
!
! NOTE: this is the "old" version prior to the parallel two-step
!       algorithm. It is still used for the "old" algorithms.
!
! Return codes, IRETURN:
!
!    0 -- successful execution
!    1 -- decomposition failed
!    2 -- memory has been out of bounds

use Para_Info, only: nProcs, Is_Real_Par
use ChoSwp, only: Diag, Diag_G, Diag_Hidden, Diag_G_Hidden
use ChoSubScr, only: Cho_SScreen
use stdalloc

implicit real*8(a-h,o-z)
#include "cholesky.fh"
#include "choprint.fh"
logical, parameter :: LOCDBG = .false.
logical, parameter :: SKIP_PRESCREEN = .false., ALLOC_BKM = .true.
logical LCONV
character(len=8), parameter :: SECNAM = 'CHO_DRV_'
real*8, parameter :: DUMTST = 0.123456789d0, DUMTOL = 1.0D-15
real*8, allocatable :: Check(:)
real*8, allocatable :: KWRK(:)
integer, allocatable :: KIRS1F(:)

#ifdef _DEBUGPRINT_
call CHO_PRTMAXMEM('CHO_DRV_ [ENTER]')
#endif

! Start overall timing.
! ---------------------

if (IPRINT >= INF_TIMING) call CHO_TIMER(TCPU0,TWALL0)

! Set return code.
! ----------------

IRETURN = 0

! Make a dummy allocation.
! ------------------------

call mma_allocate(Check,1,Label='Check')
Check(1) = DUMTST

! INITIALIZATION.
! ===============

ISEC = 1
if (IPRINT >= INF_TIMING) call CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
call CHO_INIT(SKIP_PRESCREEN,ALLOC_BKM)
call CHO_GASYNC()
if (IPRINT >= INF_TIMING) then
  call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
  call CHO_PRTTIM('Cholesky initialization',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),1)
end if
#ifdef _DEBUGPRINT_
call CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_INIT]')
IRC = 0
call CHO_DUMP(IRC,LUPRI)
if (IRC /= 0) then
  write(LUPRI,*) SECNAM,': CHO_DUMP returned ',IRC
  call CHO_QUIT('[1] Error detected in CHO_DUMP',103)
end if
#endif

! GET DIAGONAL.
! =============

ISEC = 2
if (IPRINT >= INF_TIMING) then
  call CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
  write(LUPRI,'(/,A)') '***** Starting Cholesky diagonal setup *****'
  call CHO_FLUSH(LUPRI)
end if
call CHO_GETDIAG(LCONV)
call CHO_GASYNC()
if (IPRINT >= INF_TIMING) then
  call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
  call CHO_PRTTIM('Cholesky diagonal setup',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),1)
end if
#ifdef _DEBUGPRINT_
call CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_GETDIAG]')
IRC = 0
call CHO_DUMP(IRC,LUPRI)
if (IRC /= 0) then
  write(LUPRI,*) SECNAM,': CHO_DUMP returned ',IRC
  call CHO_QUIT('[2] Error detected in CHO_DUMP',103)
end if
call CHO_PRINTLB() ! print vector dimension on each node
#endif

! DECOMPOSITION.
! ==============

ISEC = 3
if (LCONV) then
  if (RSTCHO) then
    write(LUPRI,'(//,10X,A,A,A,//)') '***** ',SECNAM,': restarted calculation converged. *****'
    call FZERO(TIMSEC(1,ISEC),4)
  else
    write(LUPRI,'(A,A)') SECNAM,': logical error: converged but not restart?!?!'
    call CHO_QUIT('Error in '//SECNAM,103)
  end if
else
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
    write(LUPRI,'(/,A)') '***** Starting Cholesky decomposition *****'
    call CHO_FLUSH(LUPRI)
  end if
  call CHO_P_SETADDR()

  if (CHO_SSCREEN) call CHO_SUBSCR_INIT()

  call CHO_DECDRV(Diag)
  call CHO_GASYNC()
  if (CHO_DECALG == 2) then
    ! generate vectors from map
    call CHO_P_OPENVR(2) ! close files
    call CHO_P_OPENVR(1) ! re-open (problem on dec-alpha)
    if (IPRINT >= INF_TIMING) then
      call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
      call CHO_PRTTIM('Cholesky map generation',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),2)
    end if
    IRC = 0
    call CHO_X_GENVEC(IRC,Diag)
    call CHO_GASYNC()
    if (IRC /= 0) then
      write(LUPRI,'(A,A)') SECNAM,': decomposition failed!'
      write(LUPRI,'(A,A,I9)') SECNAM,': CHO_X_GENVEC returned ',IRC
      IRETURN = 1
      call CHO_QUIT('Error',104)
    end if
    if (IPRINT >= INF_TIMING) then
      call CHO_TIMER(TC,TW)
      call CHO_PRTTIM('Cholesky vector generation',TC,TIMSEC(2,ISEC),TW,TIMSEC(4,ISEC),2)
    end if
  end if
  if (CHO_SSCREEN) call CHO_SUBSCR_FINAL()
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
    call CHO_PRTTIM('Cholesky decomposition',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),1)
  end if
end if
#ifdef _DEBUGPRINT_
call CHO_PRTMAXMEM('CHO_DRV_ [AFTER DECOMPOSITION]')
#endif

! CHECK DIAGONAL.
! ===============

ISEC = 4
if (LCONV) then
  call FZERO(TIMSEC(1,ISEC),4)
else
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
    write(LUPRI,'(/,A)') '***** Starting Cholesky diagonal check *****'
    call CHO_FLUSH(LUPRI)
  end if
  call mma_maxDBLE(LWRK)
  call mma_allocate(KWRK,LWRK,Label='KWRK')
  call CHO_RESTART(Diag,KWRK,LWRK,.true.,LCONV)
  call CHO_GASYNC()
  call mma_deallocate(KWRK)
  if (.not. LCONV) then
    write(LUPRI,'(A,A)') SECNAM,': Decomposition failed!'
    IRETURN = 1
    call CHO_QUIT('Decomposition failed!',104)
  end if
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
    call CHO_PRTTIM('Cholesky diagonal check',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),1)
  end if
# ifdef _DEBUGPRINT_
  call CHO_PRTMAXMEM('CHO_DRV_ [AFTER DIAGONAL CHECK]')
# endif
end if

! PARALLEL RUNS: WRITE GLOBAL DIAGONAL TO DISK.
! =============================================

call CHO_P_WRDIAG()

! CHECK INTEGRALS.
! ================

ISEC = 5
if (CHO_INTCHK) then
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
    write(LUPRI,'(/,A)') '***** Starting Cholesky integral check *****'
    call CHO_FLUSH(LUPRI)
  end if
  call CHO_DBGINT()
  call CHO_GASYNC()
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
    call CHO_PRTTIM('Cholesky integral check',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),1)
  end if
# ifdef _DEBUGPRINT_
  call CHO_PRTMAXMEM('CHO_DRV_ [AFTER INTEGRAL CHECK]')
# endif
else
  call FZERO(TIMSEC(1,ISEC),4)
end if

! REORDER VECTORS.
! ================

ISEC = 6
if (CHO_REORD) then
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
    write(LUPRI,'(/,A)') '***** Starting vector reordering *****'
    call CHO_FLUSH(LUPRI)
  end if
  LIRS1F = NNBSTRT(1)*3
  call mma_allocate(KIRS1F,LIRS1F,Label='KIRS1F')
  call mma_maxDBLE(LWRK)
  call mma_allocate(KWRK,LWRK,Label='KWRK')
  call CHO_REOVEC(KIRS1F,3,NNBSTRT(1),KWRK,LWRK)
  call CHO_GASYNC()
  call mma_deallocate(KWRK)
  call mma_deallocate(KIRS1F)
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
    call CHO_PRTTIM('Vector reordering',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),1)
  end if
# ifdef _DEBUGPRINT_
  call CHO_PRTMAXMEM('CHO_DRV_ [AFTER VECTOR REORDERING]')
# endif
else
  call FZERO(TIMSEC(1,ISEC),4)
end if

! FAKE PARALLEL: DISTRIBUTE VECTORS.
! Note: after this section, InfVec(*,3,*) is changed
! => vector disk addresses are screwed up!!
! ==================================================

ISEC = 7
if (CHO_FAKE_PAR .and. (NPROCS > 1) .and. Is_Real_Par()) then
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
    write(LUPRI,'(/,A)') '***** Starting vector distribution *****'
    call CHO_FLUSH(LUPRI)
  end if
  call CHO_PFAKE_VDIST()
  call CHO_P_WRRSTC(XNPASS)
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
    call CHO_PRTTIM('Vector distribution',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),1)
  end if
# ifdef _DEBUGPRINT_
  call CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_PFAKE_VDIST]')
# endif
else
  call FZERO(TIMSEC(1,ISEC),4)
end if

! FINALIZATIONS.
! ==============

ISEC = 8
if (IPRINT >= INF_TIMING) then
  call CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
  write(LUPRI,'(/,A)') '***** Starting Cholesky finalization *****'
  call CHO_FLUSH(LUPRI)
end if
call CHO_TRCIDL_FINAL()
call CHO_FINAL(.true.)
call CHO_GASYNC()
if (IPRINT >= INF_TIMING) then
  call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
  call CHO_PRTTIM('Cholesky finalization',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),1)
end if
#ifdef _DEBUGPRINT_
call CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_FINAL]')
#endif

! STATISTICS.
! ===========

if (IPRINT >= 1) then
  ISEC = 9
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
    write(LUPRI,'(/,A)') '***** Starting Cholesky statistics *****'
    call CHO_FLUSH(LUPRI)
  end if
  call CHO_P_STAT()
  call CHO_GASYNC()
  if (IPRINT >= INF_TIMING) then
    call CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
    call CHO_PRTTIM('Cholesky statistics',TIMSEC(2,ISEC),TIMSEC(1,ISEC),TIMSEC(4,ISEC),TIMSEC(3,ISEC),1)
  end if
end if
#ifdef _DEBUGPRINT_
call CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_STAT]')
#endif

! Close vector and reduced storage files as well as restart files.
! Deallocate all memory and test bound.
! Print total timing.
! ----------------------------------------------------------------

call CHO_P_OPENVR(2)

TST = DUMTST-Check(1)
if (abs(TST) > DUMTOL) then
  write(LUPRI,*) SECNAM,': memory has been out of bounds!!!'
  call CHO_FLUSH(LUPRI)
  IRETURN = 2
end if

if (allocated(Diag_Hidden)) call mma_deallocate(Diag_Hidden)
if (allocated(Diag_G_Hidden)) call mma_deallocate(Diag_G_Hidden)
DIag => null()
Diag_G => null()
call mma_deallocate(Check)

if (IPRINT >= INF_TIMING) then
  call CHO_TIMER(TCPU1,TWALL1)
  call CHO_PRTTIM('Cholesky procedure',TCPU1,TCPU0,TWALL1,TWALL0,1)
end if

#ifdef _DEBUGPRINT_
call CHO_PRTMAXMEM('CHO_DRV_ [EXIT]')
#endif

end subroutine CHO_DRV_Inner
