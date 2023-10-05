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

subroutine CHO_GETDIAG(LCONV)
!
! Purpose: get diagonal in first reduced set. On exit, DIAG
!          points to the diagonal and flag LCONV tells
!          if the diagonal is converged.

use Cholesky, only: Cho_IOVec, Cho_SimRI, Diag, Diag_Hidden, Frac_ChVBuf, IndRed, IndRed_Hidden, IndRSh, IndRSh_Hidden, INF_PASS, &
                    IPRINT, iSimRI, iSP2F, lBuf, LuPri, mmBstrT, Mx2Sh, MySP, n_MySP, nnBstRT, nnShl, nSym, RstCho, RstDia, &
                    Thr_SimRI
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(out) :: LCONV
integer(kind=iwp) :: IOPT, IPRTRED, ISP, ISYLST(8), ISYM, l_MySP, LMAX, LSCR, LWRK, NBIN, NDUMP, NEEDI, NEEDR, NERR
real(kind=wp) :: BIN1, STEP
logical(kind=iwp) :: DODUMMY, SYNC
integer(kind=iwp), allocatable :: KIBUF(:)
real(kind=wp), allocatable :: KBUF(:), KSCR(:), KWRK(:)
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_GETDIAG'

if (RSTDIA) then

  ! Set mySP list.
  ! Always the trivial list for serial runs (restart not possible
  ! in parallel runs).
  ! -------------------------------------------------------------

  N_MYSP = NNSHL
  l_MySP = 0
  if (allocated(MYSP)) l_MySP = size(MySP)
  if (l_MYSP == NNSHL) then
    do ISP=1,N_MYSP
      MYSP(ISP) = ISP
    end do
  else
    call CHO_QUIT('MYSP allocation error in '//SECNAM,101)
  end if

  ! Read index array NNBSTRSH and set IIBSTRSH etc.
  ! -----------------------------------------------

  call CHO_RSTD_GETIND1()

  ! Allocate mapping arrays between reduced sets.
  ! ---------------------------------------------

  MMBSTRT = NNBSTRT(1)

  call mma_allocate(IndRed_Hidden,NNBSTRT(1),3,Label='IndRed_Hidden')
  IndRed => IndRed_Hidden

  call mma_allocate(IndRSh_Hidden,NNBSTRT(1),Label='IndRSh_Hidden')
  IndRSh => IndRSh_Hidden

  ! Read mapping arrays.
  ! --------------------

  call CHO_RSTD_GETIND2()

  ! Check reduced to full shell pair mapping with the one on disk.
  ! --------------------------------------------------------------

  NERR = -1
  call CHO_RSTD_CHKSP2F(iSP2F,size(iSP2F),NERR)
  if (NERR /= 0) then
    write(LUPRI,*) SECNAM,': ',NERR,' errors detected in reduced-to-full shell pair mapping!'
    call CHO_QUIT('SP2F error in '//SECNAM,102)
  end if

  ! Allocation: diagonal.
  ! ---------------------

  NEEDR = 1
  NEEDI = 4*NEEDR

  call mma_allocate(Diag_Hidden,NNBSTRT(1),Label='Diag_Hidden')
  call mma_allocate(KBUF,NEEDR,Label='KBUF')
  call mma_allocate(KIBUF,NEEDI,Label='KIBUF')

  ! Read diagonal.
  ! --------------

  call CHO_GETDIAG1(Diag_Hidden,KBUF,KIBUF,NEEDR,NDUMP)

  call mma_deallocate(KIBUF)
  call mma_deallocate(KBUF)

else

  ! Calculate diagonal and get 1st reduced set.
  ! -------------------------------------------

  call mma_maxDBLE(LMAX)
  LMAX = LMAX/2-MX2SH
  if (LMAX < 5*LBUF) LBUF = max(LMAX/5,1)

  LSCR = MX2SH
  NEEDR = LBUF
  NEEDI = 4*LBUF

  call mma_allocate(KBUF,NEEDR,Label='KBUF')
  call mma_allocate(KSCR,LSCR,Label='KSCR')
  call mma_allocate(KIBUF,NEEDI,Label='KIBUF')

  NDUMP = 0

  call CHO_CALCDIAG(KBUF,KIBUF,LBUF,KSCR,LSCR,NDUMP)

  call mma_deallocate(KIBUF)
  call mma_deallocate(KBUF)
  call mma_deallocate(KSCR)

  ! Allocate diagonal and mapping array between reduced sets.
  ! Reallocate buffer.
  ! ---------------------------------------------------------

  MMBSTRT = NNBSTRT(1)
  call mma_allocate(IndRed_Hidden,NNBSTRT(1),3,Label='IndRed_Hidden')
  IndRed => IndRed_Hidden
  call mma_allocate(IndRSh_Hidden,NNBSTRT(1),Label='IndRSh_Hidden')
  IndRSh => IndRSh_Hidden
  call mma_allocate(Diag_Hidden,NNBSTRT(1),Label='Diag_Hidden')

  NEEDR = LBUF
  NEEDI = 4*LBUF
  call mma_allocate(KBUF,NEEDR,Label='KBUF')
  call mma_allocate(KIBUF,NEEDI,Label='KIBUF')

  ! Get diagonal in first reduced set.
  ! ----------------------------------

  call CHO_GETDIAG1(Diag_Hidden,KBUF,KIBUF,LBUF,NDUMP)

  ! Deallocate back to and including buffer.
  ! ----------------------------------------

  call mma_deallocate(KIBUF)
  call mma_deallocate(KBUF)

end if

! Set local and global info. On exit, Diag points to the local diagonal.
! ----------------------------------------------------------------------

call CHO_P_SETGL()

! Write local diagonal to disk.
! -----------------------------

IOPT = 1
call CHO_IODIAG(Diag,IOPT)

! Allocate memory for iscratch array for reading vectors.
! -------------------------------------------------------

DODUMMY = .not. ((CHO_IOVEC == 1) .or. (CHO_IOVEC == 2) .or. (CHO_IOVEC == 3) .or. (CHO_IOVEC == 4) .or. &
                 ((FRAC_CHVBUF > Zero) .and. (FRAC_CHVBUF < One)))
call CHO_ALLO_ISCR(DODUMMY)

! Initialize reduced set dimension(s) used for reading vectors.
! -------------------------------------------------------------

call CHO_INIRSDIM()

! For RI simulation, zero 1-center diagonals smaller than
! THR_SIMRI. Indices of zeroed diagonals are stored in ip_ISIMRI.
! ---------------------------------------------------------------

if (CHO_SIMRI) then
  call mma_allocate(iSimRI,NNBSTRT(1),Label='iSimRI')
  call CHO_SIMRI_Z1CDIA(Diag,THR_SIMRI,ISIMRI)
end if

! Update diagonal if restart, else just do analysis.
! --------------------------------------------------

LCONV = .false.
if (RSTCHO) then
  call mma_maxDBLE(LWRK)
  call mma_allocate(KWRK,LWRK,Label='KWRK')
  if (LOCDBG) then
    write(LUPRI,*) SECNAM,': restart diagonal:'
    do ISYM=1,NSYM
      ISYLST(ISYM) = ISYM
    end do
    call CHO_PRTDIA(Diag,ISYLST,NSYM,1)
  end if
  call CHO_RESTART(Diag,KWRK,LWRK,.false.,LCONV)
  call mma_deallocate(KWRK)
  IPRTRED = 2  ! print flag for cho_prtred
else
  if (IPRINT >= INF_PASS) then
    BIN1 = 1.0e2_wp
    STEP = 1.0e-1_wp
    NBIN = 18
    SYNC = .false.
    call CHO_P_ANADIA(Diag,SYNC,BIN1,STEP,NBIN,.true.)
  end if
  IPRTRED = 1  ! print flag for cho_prtred
end if
if (IPRINT >= INF_PASS) call CHO_P_PRTRED(IPRTRED)

end subroutine CHO_GETDIAG
