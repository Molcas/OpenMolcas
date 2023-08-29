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

subroutine CHO_RESTART(DIAG,WRK,LWRK,DSKDIA,LCONV)
!
! Purpose: update and analyze diagonal for restart
!          (or check of decomposition). Index arrays
!          for first reduced set must be set up before
!          this routine is called. Reduced set 2, on the
!          other hand, is set up here.

use Cholesky, only: Cho_1Center, Cho_DecAlg, Cho_MinChk, Cho_SimRI, IABMNZ, iAtomShl, iiBstR, iiBstRSh, IndRed, IndRSh, INF_PASS, &
                    IPRINT, iSimRI, iSP2F, LuPri, MaxRed, MySP, nDimRS, nnBstR, nnBstRSh, nnBstRT, nnShl, nSym, NumCho, ScDiag, &
                    Thr_SimRI, ThrCom, XnPass
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: Diag(*)
integer(kind=iwp), intent(in) :: LWRK
real(kind=wp), intent(out) :: WRK(LWRK)
logical(kind=iwp), intent(in) :: DSKDIA
logical(kind=iwp), intent(out) :: LCONV
integer(kind=iwp) :: I0AB, IAB, IAB1, IAB2, IMNAB, IMXAB, IOPT, ISHLA, ISHLAB, ISHLB, ISYM, JAB, JMXAB, JSHLAB, KDIAG, KEND1, &
                     KOFF, KRED, LWRK1, NBIN, NCONV, NCONVT, NDIM, NNEG, NNEGT, NSCR, NTOT, NVEC
real(kind=wp) :: AVEERR, BIN1, DMX, ERR, ERRMN, ERRMX, EXAMN, EXAMX, RMSERR, SAV, STEP, XAMAX, XDIM, XMAX, XMIN, XX
logical(kind=iwp) :: SCDIAG_SAVE, SYNC
character(len=*), parameter :: SECNAM = 'CHO_RESTART'
integer(kind=iwp), external :: CHO_F2SP
real(kind=wp), external :: ddot_

! Read diagonal (in reduced set 1).
! ---------------------------------

if (DSKDIA) then
  IOPT = 2
  call CHO_IODIAG(DIAG,IOPT)
  call CHO_P_SYNCDIAG(DIAG,1)
end if

if (IPRINT >= INF_PASS) write(LUPRI,'(/,A,I10,/)') 'Number of diagonal elements (1st reduced set): ',NNBSTRT(1)

! Analyze diagonal before update.
! -------------------------------

if (IPRINT >= INF_PASS) then
  BIN1 = 1.0e2_wp
  STEP = 1.0e-1_wp
  NBIN = 18
  SYNC = .false.
  call CHO_P_ANADIA(DIAG,SYNC,BIN1,STEP,NBIN,.true.)
end if

! Copy reduced set 1 to 2.
! ------------------------

call CHO_RSCOPY(1,2)

IMXAB = 0
IMNAB = 0
NCONVT = 0
do ISYM=1,NSYM

  NDIM = NNBSTR(ISYM,2)
  NVEC = NUMCHO(ISYM)

  if (IPRINT >= INF_PASS) then
    write(LUPRI,'(//,A,I2)') 'Check information, symmetry',ISYM
    write(LUPRI,'(/,A,6X,I12)') 'Dimension, 1st reduced set: ',NDIM
    write(LUPRI,'(A,6X,I12)') 'Number of Cholesky vectors: ',NVEC
  end if

  if ((NVEC > 0) .and. (NDIM > 0)) then

    KDIAG = 1
    KEND1 = KDIAG+NDIM
    LWRK1 = LWRK-KEND1+1

    if (LWRK1 <= 0) call CHO_QUIT('Insufficient memory in '//SECNAM,101)

    ! Save a copy of the original diagonal.
    ! -------------------------------------

    WRK(KDIAG:KDIAG+NDIM-1) = DIAG(IIBSTR(ISYM,1)+1:IIBSTR(ISYM,1)+NDIM)

    ! Calculate Cholesky diagonal.
    ! ----------------------------

    call CHO_DIACHO(DIAG,ISYM,WRK(KEND1),LWRK1)

    ! Find min. and max. error and save original value.
    ! -------------------------------------------------

    ERRMX = -1.0e10_wp
    ERRMN = 1.0e10_wp
    EXAMX = Zero
    EXAMN = Zero
    do JAB=1,NDIM
      IAB = INDRED(IIBSTR(ISYM,2)+JAB,2)
      SAV = WRK(KDIAG-1+IAB-IIBSTR(ISYM,1))
      ERR = abs(DIAG(IAB))
      if (ERR > ERRMX) then
        ERRMX = ERR
        EXAMX = SAV
        IMXAB = IAB
      end if
      if (ERR < ERRMN) then
        ERRMN = ERR
        EXAMN = SAV
        IMNAB = IAB
      end if
    end do

    ! Find min. and max. diagonals, zero too negative diagonals,
    ! and screen (if requested).
    ! ----------------------------------------------------------

    if (CHO_DECALG == 4) then
      SCDIAG_SAVE = SCDIAG
      SCDIAG = .false. ! do NOT screen (no zeroing of diags)
      DMX = One
      call CHO_CHKDIA_A4(DIAG,DMX,ISYM,NNEG,NNEGT,NSCR,XMAX,XMIN,XAMAX)
      SCDIAG = SCDIAG_SAVE
    else
      call CHO_CHKDIA(DIAG,ISYM,XMIN,XMAX,XAMAX,NNEGT,NNEG,NSCR)
    end if

    ! Count converged diagonals.
    ! --------------------------

    NCONV = 0
    do JAB=1,NDIM
      IAB = INDRED(IIBSTR(ISYM,2)+JAB,2)
      if (abs(DIAG(IAB)) <= THRCOM) NCONV = NCONV+1
    end do

    ! Calculate average and RMS error.
    ! --------------------------------

    KOFF = IIBSTR(ISYM,1)+1
    XDIM = real(NDIM,kind=wp)
    RMSERR = sqrt(DDOT_(NDIM,DIAG(KOFF),1,DIAG(KOFF),1)/XDIM)
    AVEERR = sum(DIAG(KOFF:KOFF+NDIM-1))/XDIM

    ! Print.
    ! ------

    if (IPRINT >= INF_PASS) then
      write(LUPRI,'(A,1P,D18.8)') 'Minimum diagonal          : ',XMIN
      write(LUPRI,'(A,1P,D18.8)') 'Maximum diagonal          : ',XMAX
      write(LUPRI,'(A,1P,D18.8,1X,D18.8)') 'Minimum absolute error    : ',ERRMN,EXAMN
      write(LUPRI,'(A,1P,D18.8,1X,D18.8)') 'Maximum absolute error    : ',ERRMX,EXAMX
      write(LUPRI,'(A,1P,D18.8)') 'Average error             : ',AVEERR
      write(LUPRI,'(A,1P,D18.8)') 'Root-mean-square error    : ',RMSERR
      write(LUPRI,'(A,6X,I12)') 'Converged diagonals       : ',NCONV
      write(LUPRI,'(A,6X,I12)') 'Unconverged diagonals     : ',NDIM-NCONV
      write(LUPRI,'(A,6X,I12)') 'Zeroed negative diagonals : ',NNEG
      if (CHO_DECALG /= 4) then ! NSCR is useless here
        if (SCDIAG) then
          write(LUPRI,'(A,6X,I12)') 'Screened diagonals        : ',NSCR
        else
          write(LUPRI,'(A,6X,I12,A)') 'Screenable diagonals      : ',NSCR,' (not screened)'
        end if
      end if
    end if
    NCONVT = NCONVT+NCONV

  end if

  call XFLUSH(LUPRI)

end do

! Sync diagonal and set reduced set.
! ----------------------------------

SYNC = .true.
call CHO_P_SETRED(DIAG,SYNC)
KRED = XNPASS+1
call CHO_SETRSDIM(NDIMRS,NSYM,MAXRED,KRED,2)

! Sync and analyze (histogram) updated diagonal.
! ----------------------------------------------

call CHO_P_SYNCDIAG(DIAG,2)
if (IPRINT >= INF_PASS) then
  BIN1 = 1.0e2_wp
  STEP = 1.0e-1_wp
  NBIN = 18
  SYNC = .false.
  call CHO_P_ANADIA(DIAG,SYNC,BIN1,STEP,NBIN,.false.)
end if

! Set data for minimal integral checking.
! ---------------------------------------

if (CHO_MINCHK) then
  if (IMXAB > 0) then
    ISHLAB = CHO_F2SP(INDRSH(IMXAB))
    if (ISHLAB > 0) then
      call CHO_INTCHK_REG('MAX DIAG',ISHLAB,ISHLAB)
    else
      call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
    end if
  end if
  if (IMNAB > 0) then
    ISHLAB = CHO_F2SP(INDRSH(IMNAB))
    if (ISHLAB > 0) then
      call CHO_INTCHK_REG('MIN DIAG',ISHLAB,ISHLAB)
    else
      call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
    end if
  end if
  if ((IMXAB > 0) .and. (IMNAB > 0)) then
    ISHLAB = CHO_F2SP(INDRSH(IMXAB))
    JSHLAB = CHO_F2SP(INDRSH(IMNAB))
    if ((ISHLAB > 0) .and. (JSHLAB > 0)) then
      call CHO_INTCHK_REG('MAX|MIN ',ISHLAB,JSHLAB)
    else
      call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
    end if
  end if
  JMXAB = 0
  I0AB = 0
  XX = Zero
  do ISHLAB=1,NNSHL
    NTOT = 0
    do ISYM=1,NSYM
      if (NNBSTRSH(ISYM,ISHLAB,1) > 0) then
        IAB1 = IIBSTR(ISYM,1)+IIBSTRSH(ISYM,ISHLAB,1)+1
        IAB2 = IAB1+NNBSTRSH(ISYM,ISHLAB,1)-1
        do IAB=IAB1,IAB2
          if (DIAG(IAB) < XX) then
            XX = DIAG(IAB)
            JMXAB = IAB
          end if
        end do
        NTOT = NTOT+1
      end if
    end do
    if (NTOT == 0) I0AB = ISHLAB
  end do
  if ((JMXAB > 0) .and. (JMXAB /= IMNAB)) then
    JSHLAB = CHO_F2SP(INDRSH(JMXAB))
    if (JSHLAB > 0) then
      call CHO_INTCHK_REG('NEG DIAG',JSHLAB,JSHLAB)
    else
      call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
    end if
    if (IMXAB > 0) then
      ISHLAB = CHO_F2SP(INDRSH(IMXAB))
      if (ISHLAB > 0) then
        call CHO_INTCHK_REG('MAX|NEG ',ISHLAB,JSHLAB)
      else
        call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
      end if
    end if
    if (IMNAB > 0) then
      ISHLAB = CHO_F2SP(INDRSH(IMNAB))
      if (ISHLAB > 0) then
        call CHO_INTCHK_REG('MIN|NEG ',ISHLAB,JSHLAB)
      else
        call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
      end if
    end if
  end if
  if (I0AB > 0) then
    JSHLAB = CHO_F2SP(INDRSH(I0AB))
    if (JSHLAB > 0) then
      call CHO_INTCHK_REG('EXCL RS1',JSHLAB,JSHLAB)
    else
      call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
    end if
    if (IMXAB > 0) then
      ISHLAB = CHO_F2SP(INDRSH(IMXAB))
      if (ISHLAB > 0) then
        call CHO_INTCHK_REG('MAX|XRS1',ISHLAB,JSHLAB)
      else
        call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
      end if
    end if
    if (IMNAB > 0) then
      ISHLAB = CHO_F2SP(INDRSH(IMNAB))
      if (ISHLAB > 0) then
        call CHO_INTCHK_REG('MIN|XRS1',ISHLAB,JSHLAB)
      else
        call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
      end if
    end if
  end if
  if ((IABMNZ > 0) .and. (IABMNZ <= NNSHL)) then
    JSHLAB = CHO_F2SP(INDRSH(IABMNZ))
    if (JSHLAB > 0) then
      call CHO_INTCHK_REG('NEG->ZER',JSHLAB,JSHLAB)
    else
      call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
    end if
    if (IMXAB > 0) then
      ISHLAB = CHO_F2SP(INDRSH(IMXAB))
      if (ISHLAB > 0) then
        call CHO_INTCHK_REG('MAX|NEGZ',ISHLAB,JSHLAB)
      else
        call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
      end if
    end if
    if (IMNAB > 0) then
      ISHLAB = CHO_F2SP(INDRSH(IMNAB))
      if (ISHLAB > 0) then
        call CHO_INTCHK_REG('MIN|NEGZ',ISHLAB,JSHLAB)
      else
        call CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
      end if
    end if
  end if
end if

! Set convergence flag.
! ---------------------

LCONV = NCONVT == NNBSTRT(1)

! For 1-center decomposition, convergence is determined only by
! judging the 1-center diagonals.
! For RI simulations, diagonals that were initially zeroed may not
! be converged, but that is OK as long as the diagonal is smaller
! than the threshold for deletion, THR_SIMRI.
! ----------------------------------------------------------------

if (CHO_1CENTER .and. (.not. LCONV)) then
# ifdef _DEBUGPRINT_
  if (NSYM /= 1) call CHO_QUIT(SECNAM//': CHO_1CENTER on, but NSYM != 1',103)
  if (size(IATOMSHL) < NSHELL) call CHO_QUIT(SECNAM//': iAtomShl not allocated correctly!',103)
# endif
  LCONV = .true.
  ISHLAB = 0
  do while ((ISHLAB < NNSHL) .and. LCONV)
    ISHLAB = ISHLAB+1
    call CHO_INVPCK(ISP2F(MYSP(ISHLAB)),ISHLA,ISHLB,.true.)
    if (IATOMSHL(ISHLA) == IATOMSHL(ISHLB)) then
      IAB1 = IIBSTRSH(1,ISHLAB,1)+1
      IAB2 = IAB1+NNBSTRSH(1,ISHLAB,1)-1
      NCONV = 0
      do IAB=IAB1,IAB2
        if (abs(DIAG(IAB)) <= THRCOM) then
          NCONV = NCONV+1
        else if (CHO_SIMRI) then
          if ((ISIMRI(IAB) == 1) .and. (abs(DIAG(IAB)) <= THR_SIMRI)) NCONV = NCONV+1
        end if
      end do
      LCONV = NCONV == NNBSTRSH(1,ISHLAB,1)
    end if
  end do
end if

end subroutine CHO_RESTART
