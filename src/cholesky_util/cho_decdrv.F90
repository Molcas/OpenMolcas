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

subroutine CHO_DECDRV(DIAG)
!
! Purpose: driver for the decomposition of the two-electron integral
!          matrix based on the reduced diagonal.

use Cholesky, only: CHO_DECALG, Cho_Real_Par, CHO_SIMP, DIAMIN, DID_DECDRV, FRAC_CHVBUF, INF_PASS, INF_VECBUF, InfRed, IPRINT, &
                    LuPri, LuSel, MaxRed, nDimRS, nnBstR, nnBstRT, nnBstRT_G, nnShl, nnShl_G, nSym, NumCho, TDECDRV, ThrCom, &
                    Trace_Idle, Span, XnPass
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Diag(*)
integer(kind=iwp) :: I, IPASS, IPASS_PREV, IRC, IRED, ISYLST(8), iSym, JPASS, KRED, LWRK, MPASS, NBIN, nDim_Now, NGSP, NPOTSH, NUM
real(kind=wp) :: BIN1, DIAMAX_SIMP(8), STEP, TCPU1, TCPU2, TLDEC, TLDEC1, TLDEC2, TLINT, TLINT1, TLINT2, TLTOT, TLTOT1, TLTOT2, &
                 TWALL1, TWALL2, WLDEC, WLDEC1, WLDEC2, WLINT, WLINT1, WLINT2, WLTOT, WLTOT1, WLTOT2
logical(kind=iwp) :: CONV, SYNC
character(len=20) :: STRING
character(len=7) :: FILSEL
integer(kind=iwp), allocatable :: KISYSH(:), LSTQSP(:)
real(kind=wp), allocatable :: KDIASH(:), KWRK(:)
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_DECDRV'

! Start timing.
! -------------

call CWTIME(TCPU1,TWALL1)

! Initializations and static settings.
! IRED=2: points to current reduced set in index arrays.
! ------------------------------------------------------

IRED = 2
CONV = .false.

! Initialize Cholesky vector buffer.
! ----------------------------------

call CHO_VECBUF_INIT(FRAC_CHVBUF,NNBSTR(1,1))
if (LOCDBG .or. (IPRINT >= INF_VECBUF)) call CHO_VECBUF_PRINT(LUPRI,NSYM)

! Allocate memory for shell pair based diagonal.
! It is important that the DIASH allocation is first, as memory is
! released by flushing back to and including this allocation.
! ----------------------------------------------------------------

if (Cho_Real_Par) then
  NGSP = nnShl_G
else
  NGSP = nnShl
end if
call mma_allocate(KDIASH,NGSP,Label='KDIASH')
call mma_allocate(KISYSH,NGSP,Label='KISYSH')

! Set first integral pass.
! ------------------------

SYNC = .false.
NPOTSH = 0
call CHO_P_SETPASS(DIAG,SYNC,KDIASH,KISYSH,IRED,CONV,NPOTSH)
if (NPOTSH > 0) then
  if (CONV) call CHO_QUIT('Logical error [0.1] in '//SECNAM,103)
else
  if (.not. CONV) call CHO_QUIT('Logical error [0.2] in '//SECNAM,103)
end if

! Allocate shell pair list.
! -------------------------

call mma_allocate(LSTQSP,max(NPOTSH,1),Label='LSTQSP')

! Loop over integral passes. Continue until convergence or
! until the max. number of integral passes has been reached.
! To each integral pass there is associated a reduced set,
! so the IPASS counter is also used as identifier of reduced
! set during I/O.
! ----------------------------------------------------------

IPASS = XNPASS
JPASS = 0
if (Cho_Real_Par) then
  MPASS = nnBstRT_G(IRED)
else
  MPASS = nnBstRT(IRED)
end if
do while ((.not. CONV) .and. (JPASS < MPASS))

  ! Update integral pass counter.
  ! -----------------------------

  JPASS = JPASS+1
  IPASS = IPASS+1

  ! Print.
  ! ------

  if (IPRINT >= INF_PASS) then
    call CWTIME(TLTOT1,WLTOT1)
    write(STRING,'(A13,I7)') 'Integral Pass',IPASS
    call CHO_HEAD(STRING,'*',80,LUPRI)
  end if

  ! Update idle proc info.
  ! ----------------------

  if (Trace_Idle) then
    nDim_Now = sum(nnBstR(1:nSym,2))
    call Cho_TrcIdl_Update(nDim_Now < 1)
  end if

  ! Debug: print diagonal.
  ! ----------------------

  if (LOCDBG) then
    write(LUPRI,*) SECNAM,': debug: diagonal before pass ',IPASS
    do ISYM=1,NSYM
      ISYLST(ISYM) = ISYM
    end do
    SYNC = .false.
    call CHO_P_PRTDIA(DIAG,SYNC,ISYLST,NSYM,IRED)
    write(LUPRI,*)
    write(LUPRI,*) SECNAM,': INFRED before pass ',IPASS
    write(LUPRI,'(10I8)') (INFRED(I),I=1,min(IPASS,MAXRED))
  end if

  ! Write index arrays for reduced set to disk
  ! and update disk address.
  ! ------------------------------------------

  call CHO_P_PUTRED(IPASS,IRED)

  ! Maintain Cholesky vector buffer.
  ! The logicals request that statistics informations are updated
  ! in the maintainance routine.
  ! -------------------------------------------------------------

  IRC = 0
  IPASS_PREV = IPASS-1
  call CHO_VECBUF_MAINTAIN(IRC,IPASS_PREV,.true.,.true.)
  if (IRC /= 0) then
    write(LUPRI,*) SECNAM,': CHO_VECBUF_MAINTAIN returned ',IRC
    call CHO_QUIT('Error detected in '//SECNAM,IRC)
  end if

  ! Open scratch files for qualified integral columns.
  ! --------------------------------------------------

  do ISYM=1,NSYM
    if (NNBSTR(ISYM,2) > 0) then
      LUSEL(ISYM) = 7
      write(FILSEL,'(A6,I1)') 'CHOSEL',ISYM
      call DANAME_WA(LUSEL(ISYM),FILSEL)
    else
      LUSEL(ISYM) = -1
    end if
  end do

  ! Get integral columns on disk stored in current reduced set.
  ! -----------------------------------------------------------

  if (IPRINT >= INF_PASS) call CWTIME(TLINT1,WLINT1)
  NUM = 0
  call CHO_GETINT(DIAG,KDIASH,KISYSH,LSTQSP,NPOTSH,NUM)
  call XFLUSH(LUPRI)
  if (IPRINT >= INF_PASS) call CWTIME(TLINT2,WLINT2)

  ! Decompose the qualified integral columns.
  ! -----------------------------------------

  if (IPRINT >= INF_PASS) call CWTIME(TLDEC1,WLDEC1)
  if ((CHO_DECALG == 4) .or. (CHO_DECALG == 5) .or. (CHO_DECALG == 6)) then
    call CHO_DECOM_A4(DIAG,LSTQSP,NUM,IPASS)
  else
    if (CHO_SIMP) then
      call CHO_MAXDX(DIAG,DIAMAX_SIMP)
      do ISYM=1,NSYM
        DIAMIN(ISYM) = max(THRCOM,DIAMAX_SIMP(ISYM)*SPAN)
      end do
    end if
    call mma_maxDBLE(LWRK)
    call mma_allocate(KWRK,LWRK,Label='KWRK')
    call CHO_DECOM(DIAG,KWRK,LWRK,IPASS,NUM)
    call mma_deallocate(KWRK)
  end if
  call XFLUSH(LUPRI)
  if (IPRINT >= INF_PASS) call CWTIME(TLDEC2,WLDEC2)

  ! Sync global vector counter.
  ! ---------------------------

  call CHO_P_SYNCNUMCHO(NUMCHO,NSYM)

  ! Write restart info to disk.
  ! ---------------------------

  call CHO_P_WRRSTC(IPASS)

  ! Close scratch files for qualified integral columns.
  ! ---------------------------------------------------

  do ISYM=1,NSYM
    if (LUSEL(ISYM) > 0) call DACLOS(LUSEL(ISYM))
  end do

  ! Sync diagonal.
  ! --------------

  call CHO_P_SYNCDIAG(DIAG,2)

  ! Analyze diagonal.
  ! -----------------

  if (IPRINT >= INF_PASS) then
    BIN1 = 1.0e2_wp
    STEP = 1.0e-1_wp
    NBIN = 18
    SYNC = .false.
    call CHO_P_ANADIA(DIAG,SYNC,BIN1,STEP,NBIN,.false.)
  end if

  ! Get next reduced set.
  ! ---------------------

  SYNC = .false.
  call CHO_P_SETRED(DIAG,SYNC)
  KRED = IPASS+1
  call CHO_SETRSDIM(NDIMRS,NSYM,MAXRED,KRED,IRED)
  if (IPRINT >= INF_PASS) then
    call CHO_P_PRTRED(2)
    call XFLUSH(LUPRI)
  end if

  ! Check convergence and, if not converged, set next integral pass.
  ! ----------------------------------------------------------------

  SYNC = .false.
  NPOTSH = 0
  call CHO_P_SETPASS(DIAG,SYNC,KDIASH,KISYSH,IRED,CONV,NPOTSH)
  if (NPOTSH > 0) then
    if (CONV) call CHO_QUIT('Logical error [1.1] in '//SECNAM,103)
  else
    if (.not. CONV) call CHO_QUIT('Logical error [1.2] in '//SECNAM,103)
  end if

  ! Update bookmarks: store largest diagonal (integral accuracy)
  ! and number of Cholesky vectors.
  ! ------------------------------------------------------------

  call Cho_P_UpdateBookmarks(iPass)

  ! Print idle report.
  ! ------------------

  if (Trace_Idle) call Cho_TrcIdl_Report()

  ! Print timing for this pass.
  ! ---------------------------

  if (IPRINT >= INF_PASS) then
    TLINT = TLINT2-TLINT1
    WLINT = WLINT2-WLINT1
    TLDEC = TLDEC2-TLDEC1
    WLDEC = WLDEC2-WLDEC1
    call CWTIME(TLTOT2,WLTOT2)
    TLTOT = TLTOT2-TLTOT1
    WLTOT = WLTOT2-WLTOT1
    write(LUPRI,'(/,A,I7,A)') 'Overall timings for integral pass',IPASS,' (CPU/Wall in seconds):'
    write(LUPRI,'(A,F12.2,1X,F12.2)') 'Integrals (incl. qualified I/O etc.): ',TLINT,WLINT
    write(LUPRI,'(A,F12.2,1X,F12.2)') 'Decomposition of qualified columns  : ',TLDEC,WLDEC
    write(LUPRI,'(A,F12.2,1X,F12.2)') 'Total (incl. restart info I/O etc.) : ',TLTOT,WLTOT
  end if

end do

! Free memory for shell pair based diagonal.
! ------------------------------------------

call mma_deallocate(LSTQSP)
call mma_deallocate(KISYSH)
call mma_deallocate(KDIASH)

! Shut down the Cholesky vector buffer.
! -------------------------------------

call CHO_VECBUF_FINAL()

! Set stuff for statistics.
! -------------------------

DID_DECDRV = .true.
XNPASS = IPASS

! Timing.
! -------

call CWTIME(TCPU2,TWALL2)
TDECDRV(1) = TCPU2-TCPU1
TDECDRV(2) = TWALL2-TWALL1

end subroutine CHO_DECDRV
