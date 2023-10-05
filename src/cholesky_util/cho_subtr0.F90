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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine CHO_SUBTR0(XINT,WRK,LWRK,ISYM)
!
! Purpose: subtract contributions from previous vectors
!          from the qualified integrals (in XINT).
!          This version is memory-driven.
!
! Screening in subtraction introduced Jan. 2006, TBP.

use Cholesky, only: Cho_SScreen, DSPNm, DSubScr, iiBstR, iiBstRSh, iQuAB, LQ, LuPri, nDGM_call, nnBstR, nnBstRSh, nnShl, nQual, &
                    NumCho, nVec_in_Buf, SSNorm, SSTau, SubScrStat, TDECOM
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LWRK, ISYM
real(kind=wp), intent(inout) :: XINT(*), WRK(LWRK)
integer(kind=iwp) :: IAB, IBATCH, ISHGD, IVEC1, J, JAB, KCHO1, KCHO2, KEND2, KOFB0, KOFF1, KOFF2, KOFF3, KOFFA, KOFFB, LREAD, &
                     LWRK1, LWRK2, MMEM, NBATCH, NGD, NUMV, NVEC, NVEC_TO_READ
real(kind=wp) :: C1, C2, TST, W1, W2, XDON, XTOT
character(len=*), parameter :: SECNAM = 'CHO_SUBTR0'
integer(kind=iwp), external :: CHO_LREAD

! Return if nothing to do.
! ------------------------

if (NUMCHO(ISYM) < 1) return

NVEC_TO_READ = NUMCHO(ISYM)-NVEC_IN_BUF(ISYM)
if (NVEC_TO_READ == 0) return
if (NVEC_TO_READ < 0) call CHO_QUIT('Vector buffer error in '//SECNAM,104)

! Initialize.
! -----------

XTOT = Zero
XDON = Zero

! Reserve space needed for reading previous vectors.
! --------------------------------------------------

LREAD = CHO_LREAD(ISYM,LWRK)
if (LREAD < 1) then
  write(LUPRI,*) SECNAM,': CHO_LREAD returned ',LREAD
  call CHO_QUIT('Memory error in '//SECNAM,101)
  LWRK1 = 0 ! to avoid compiler warnings
else
  LWRK1 = LWRK-LREAD
end if

! Set up batch.
! -------------

MMEM = NNBSTR(ISYM,2)+NQUAL(ISYM)
NVEC = min(LWRK1/MMEM,NVEC_TO_READ)
if (NVEC < 1) call CHO_QUIT('Batch failure in '//SECNAM,101)
NBATCH = (NVEC_TO_READ-1)/NVEC+1

! Start batch loop.
! -----------------

do IBATCH=1,NBATCH

  if (IBATCH == NBATCH) then
    NUMV = NVEC_TO_READ-NVEC*(NBATCH-1)
  else
    NUMV = NVEC
  end if
  IVEC1 = NVEC_IN_BUF(ISYM)+NVEC*(IBATCH-1)+1

  ! Complete allocation.
  ! --------------------

  KCHO1 = 1
  KCHO2 = KCHO1+NNBSTR(ISYM,2)*NUMV
  KEND2 = KCHO2+NQUAL(ISYM)*NUMV
  LWRK2 = LWRK-KEND2+1
  if (LWRK2 < LREAD) call CHO_QUIT('Batch error in '//SECNAM,104)

  ! Read previous vectors.
  ! ----------------------

  call CWTIME(C1,W1)
  call CHO_GETVEC(WRK(KCHO1),NNBSTR(ISYM,2),NUMV,IVEC1,ISYM,WRK(KEND2),LWRK2)
  call CWTIME(C2,W2)
  TDECOM(1,2) = TDECOM(1,2)+C2-C1
  TDECOM(2,2) = TDECOM(2,2)+W2-W1

  ! Screened or unscreened subtraction section.
  ! The screened version uses level 2 blas, while the unscreened
  ! one employs level 3 blas.
  ! ------------------------------------------------------------

  call CWTIME(C1,W1)

  if (CHO_SSCREEN) then ! screened subtraction

    KOFB0 = KCHO1-1-IIBSTR(ISYM,2)
    do J=1,NUMV
      KOFFA = KCHO2+J-1
      KOFFB = KOFB0+NNBSTR(ISYM,2)*(J-1)
      do IAB=1,NQUAL(ISYM)
        WRK(KOFFA+NUMV*(IAB-1)) = WRK(KOFFB+IQUAB(IAB,ISYM))
      end do
    end do

    ! Subtract:
    ! (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L(#J,{ab})
    ! for each ab in {ab}.
    ! ----------------------------------------------------

    call CHO_SUBSCR_DIA(WRK(KCHO1),NUMV,iSym,2,SSNorm)
    do IAB=1,NQUAL(ISYM)
      do ISHGD=1,NNSHL
        NGD = NNBSTRSH(ISYM,ISHGD,2)
        if (NGD > 0) then
          XTOT = XTOT+One
          JAB = IQUAB(IAB,ISYM)-IIBSTR(ISYM,2)
          TST = sqrt(DSPNM(ISHGD)*DSUBSCR(JAB))
          if (TST > SSTAU) then
            XDON = XDON+One
            KOFF1 = KCHO1+IIBSTRSH(ISYM,ISHGD,2)
            KOFF2 = KCHO2+NUMV*(IAB-1)
            KOFF3 = NNBSTR(ISYM,2)*(IAB-1)+IIBSTRSH(ISYM,ISHGD,2)+1
            call DGEMV_('N',NGD,NUMV,-One,WRK(KOFF1),NNBSTR(ISYM,2),WRK(KOFF2),1,One,XINT(KOFF3),1)
          end if
        end if
      end do
    end do

  else ! unscreened subtraction

    if (associated(LQ(ISYM)%A)) then

      ! If the qualified block, L({ab},#J), is already in core,
      ! use this block.
      ! -------------------------------------------------------

      call DGEMM_('N','T',NNBSTR(ISYM,2),NQUAL(ISYM),NUMV,-One,WRK(KCHO1),NNBSTR(ISYM,2),LQ(ISYM)%A(:,IVEC1),size(LQ(ISYM)%A,1), &
                  One,XINT,NNBSTR(ISYM,2))

    else

      ! Copy out sub-blocks corresponding to qualified diagonals:
      ! L({ab},#J)
      ! ---------------------------------------------------------

      do J=1,NUMV
        do IAB=1,NQUAL(ISYM)
          KOFF1 = KCHO2+NQUAL(ISYM)*(J-1)+IAB-1
          KOFF2 = KCHO1+NNBSTR(ISYM,2)*(J-1)+IQUAB(IAB,ISYM)-IIBSTR(ISYM,2)-1
          WRK(KOFF1) = WRK(KOFF2)
        end do
      end do

      ! Subtract:
      ! (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L({ab},#J)
      ! ----------------------------------------------------

      call DGEMM_('N','T',NNBSTR(ISYM,2),NQUAL(ISYM),NUMV,-One,WRK(KCHO1),NNBSTR(ISYM,2),WRK(KCHO2),NQUAL(ISYM),One,XINT, &
                  NNBSTR(ISYM,2))

    end if

    ! Update DGEMM-call counter.
    ! --------------------------

    NDGM_CALL = NDGM_CALL+1

  end if

  call CWTIME(C2,W2)
  TDECOM(1,3) = TDECOM(1,3)+C2-C1
  TDECOM(2,3) = TDECOM(2,3)+W2-W1

end do

! Update screening statistics.
! ----------------------------

if (CHO_SSCREEN) then
  SUBSCRSTAT(1) = SUBSCRSTAT(1)+XTOT
  SUBSCRSTAT(2) = SUBSCRSTAT(2)+XDON
end if

end subroutine CHO_SUBTR0
