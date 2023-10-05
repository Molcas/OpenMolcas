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

subroutine CHO_SUBTR1(XINT,WRK,LWRK,ISYM,FXDMEM)
!
! Purpose: subtract contributions from previous vectors
!          from the qualified integrals (in XINT).
!          This version is I/O-driven.
!
! Screening in subtraction introduced Jan. 2006, TBP.

use Cholesky, only: Cho_SScreen, DSPNm, DSubScr, iiBstR, iiBstRSh, INF_SUBTR1, InfVec, IPRINT, iQuAB, iScr, LQ, LuPri, MaxQual, &
                    N1_Qual, N2_Qual, N1_VecRd, N2_VecRd, N_Subtr, nDGM_call, nnBstR, nnBstRSh, nnShl, nQual, nSys_call, NumCho, &
                    nVec_in_Buf, SSNorm, SSTau, SubScrStat, TDECOM
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LWRK, ISYM
real(kind=wp), intent(inout) :: XINT(*), WRK(LWRK)
logical(kind=iwp) :: FXDMEM
integer(kind=iwp) :: IAB, IBATCH, ILOC, IMAPC, IOFF(0:1), IREDC, ISHGD, IVEC1, IVEC1_1, IVSTAT(2,2), J, JAB, JRED, JRED1, JRED2, &
                     JVEC1, JVEC2, KCHO1, KCHO2, KEND0, KEND1, KEND2, KJUNK, KOFB0, KOFF1, KOFF2, KOFF3, KOFFA, KOFFB, KREAD, &
                     KVEC, KVEC1, LEFT, LNUM, LREAD, LRED, LVEC, LVEC1, LWRK0, LWRK1, LWRK2, MINLFT, MMEM, MUSED, MUST, NBATCH, &
                     NGD, NUMBAT, NUMRD, NUMSMN, NUMSUB, NUMV, NVEC, NVEC_TO_READ, NVRD
real(kind=wp) :: C1, C2, SCRPCT, TIMLOC(2,3), TST, W1, W2, X1, X2, XAVERD, XAVEVC, XDON, XM, XREAD, XTOT
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_SUBTR1'
integer(kind=iwp), external :: CHO_X_NUMRD

! Return if nothing to do.
! ------------------------

if (NUMCHO(ISYM) < 1) return

NVEC_TO_READ = NUMCHO(ISYM)-NVEC_IN_BUF(ISYM)
if (NVEC_TO_READ == 0) return
if (NVEC_TO_READ < 0) call CHO_QUIT('Vector buffer error in '//SECNAM,104)

! Allocate "junk yard".
! ---------------------

KJUNK = 1
KEND0 = KJUNK+1
LWRK0 = LWRK-KEND0+1

MUST = NNBSTR(ISYM,1)+NNBSTR(ISYM,2)+NQUAL(ISYM)
if (LWRK0 < MUST) then
  write(LUPRI,*)
  write(LUPRI,*)
  write(LUPRI,*) SECNAM,': insufficient memory:'
  write(LUPRI,*) 'Need at least: ',MUST+KEND0-1
  write(LUPRI,*) 'Available    : ',LWRK
  write(LUPRI,*) '(A significant increase of memory is needed for efficient execution.)'
  write(LUPRI,*)
  write(LUPRI,*) 'Memory available in ',SECNAM,' may also be increased by reducing:'
  write(LUPRI,*) '1) max. #qualified per symmetry, currently: ',MAXQUAL,'   to less than than ',NQUAL(ISYM)
  write(LUPRI,*) '2) max. memory fraction used by qualified,  currently: ',N1_QUAL,'/',N2_QUAL
  call CHO_QUIT('Insufficient memory in '//SECNAM//' [0]',101)
end if

WRK(KJUNK) = Zero
IOFF(0) = KJUNK

! Split memory for subtraction of previous vectors:
! {Read buffer},{L(cd,#J),L({ab},#J)}
! Initially, the fraction N1_VECRD/N2_VECRD of total memory is
! reserved for reading. Then, the split aims at
!          1<=#J<=MAX(NQUAL(ISYM),N_SUBTR)
! I.e., the number of vectors in the read buffer is at least #J.
! Obviously, if MAX(NQUAL(ISYM),N_SUBTR) > NVEC_TO_READ then
! 1<=#J<=NVEC_TO_READ. (Keeping #J within bounds is the role of
! NUMSMN).
! N1_VECRD, N2_VECRD, and N_SUBTR can be user-defined (input)
! and should have been checked as part of configuration check
! (CHO_CHKCONF) during initialization.
! --------------------------------------------------------------

KREAD = KEND0 ! pointer to read buffer
IREDC = -1    ! id of red. set index array at location 3

X1 = real(N1_VECRD,kind=wp)
X2 = real(N2_VECRD,kind=wp)
XM = real(LWRK0,kind=wp)
XREAD = XM*X1/X2  ! initial guess for buffer size (from input)
LREAD = int(XREAD)
LEFT = LWRK0-LREAD
MINLFT = NNBSTR(ISYM,2)+NQUAL(ISYM)

if (FXDMEM) then ! fixed-size read buffer

  if (LEFT < MINLFT) then
    LEFT = MINLFT
    LREAD = LWRK0-LEFT
  end if

else ! attempt to optimize the split

  NUMSMN = min(max(NQUAL(ISYM),N_SUBTR),NVEC_TO_READ)
  NUMSUB = max(min(NUMSMN,LEFT/MINLFT),1) ! 1<=NUMSUB<=NUMSMN
  LREAD = LWRK0-NUMSUB*MINLFT

  NUMRD = CHO_X_NUMRD(1,ISYM,IREDC,LREAD) ! # that can be read
  if (NUMRD < 0) call CHO_QUIT('NUMRD error in '//SECNAM,104)
  do while (NUMRD < NUMSUB) ! reduce NUMSUB until NUMRD=NUMSUB
    NUMSUB = NUMSUB-1
    if (NUMSUB < 1) call CHO_QUIT('Insufficient memory for split in '//SECNAM,101) ! should never occur (checked above)
    LREAD = LWRK0-NUMSUB*MINLFT
    NUMRD = CHO_X_NUMRD(1,ISYM,IREDC,LREAD)
    if (NUMRD < 0) call CHO_QUIT('NUMRD error in '//SECNAM,104)
  end do

end if

! Initializations.
! ----------------

NUMRD = 0
NUMBAT = 0
XTOT = Zero
XDON = Zero
IVSTAT(:,:) = 0
TIMLOC(:,:) = 0

! Start buffer batch loop.
! ------------------------

IVEC1 = NVEC_IN_BUF(ISYM)+1
IMAPC = -1
do while (IVEC1 <= NUMCHO(ISYM))

  ! Read as many vectors as possible into buffer.
  ! ---------------------------------------------

  call CWTIME(C1,W1)
  NVRD = 0
  MUSED = 0
  call CHO_VECRD(WRK(KREAD),LREAD,IVEC1,NUMCHO(ISYM),ISYM,NVRD,IREDC,MUSED)
  NUMRD = NUMRD+1
  call CWTIME(C2,W2)
  TIMLOC(1,1) = TIMLOC(1,1)+C2-C1
  TIMLOC(2,1) = TIMLOC(2,1)+W2-W1

  ! Quit if no vectors were read.
  ! -----------------------------

  if (NVRD < 1) call CHO_QUIT('Insufficient scratch space for read in '//SECNAM,101)

  ! Compute memory available for subtraction batching.
  ! --------------------------------------------------

  KEND1 = KREAD+MUSED
  LWRK1 = LWRK-KEND1+1

  if (LWRK1 < 1) call CHO_QUIT('Insufficient memory in '//SECNAM//' [1]',101)

  ! Set up batch.
  ! -------------

  MMEM = NNBSTR(ISYM,2)+NQUAL(ISYM)
  if (MMEM < 1) then
    call CHO_QUIT('Batch setup corrupted in '//SECNAM,104)
    NVEC = -999999
  else
    NVEC = min(LWRK1/MMEM,NVRD)
  end if
  if (NVEC < 1) then
    call CHO_QUIT('Batch failure in '//SECNAM,101)
    NBATCH = -999999
  else
    NBATCH = (NVRD-1)/NVEC+1
  end if

  ! Set local statistics info.
  ! --------------------------

  NUMBAT = NUMBAT+NBATCH
  if (NUMRD == 1) then
    IVSTAT(:,1) = NVRD
    IVSTAT(:,2) = NVEC
  else
    IVSTAT(1,1) = min(IVSTAT(1,1),NVRD)
    IVSTAT(2,1) = max(IVSTAT(2,1),NVRD)
    IVSTAT(1,2) = min(IVSTAT(1,2),NVEC)
    IVSTAT(2,2) = max(IVSTAT(2,2),NVEC)
  end if

  ! Start batch loop.
  ! -----------------

  IOFF(1) = KREAD-1
  do IBATCH=1,NBATCH

    if (IBATCH == NBATCH) then
      NUMV = NVRD-NVEC*(NBATCH-1)
    else
      NUMV = NVEC
    end if
    IVEC1_1 = IVEC1+NVEC*(IBATCH-1)

    ! Set memory pointers for this batch.
    ! -----------------------------------

    KCHO1 = KEND1
    KCHO2 = KCHO1+NNBSTR(ISYM,2)*NUMV
    KEND2 = KCHO2+NQUAL(ISYM)*NUMV
    LWRK2 = LWRK-KEND2+1
    if (LWRK2 < 0) call CHO_QUIT('Batch error in '//SECNAM,104)

    ! Get the next NUMV vectors sorted according to current
    ! reduced set (originally, this section was part of the
    ! I/O; hence, it is timed as if it was I/O for
    ! compatibility).
    ! -----------------------------------------------------

    call CWTIME(C1,W1)

    JVEC1 = NVEC*(IBATCH-1)+1
    JVEC2 = JVEC1+NUMV-1
    JRED1 = INFVEC(IVEC1+JVEC1-1,2,ISYM)
    JRED2 = INFVEC(IVEC1+JVEC2-1,2,ISYM)
    LVEC1 = JVEC1
    KVEC1 = 1
    do JRED=JRED1,JRED2

      LNUM = 0
      LVEC = LVEC1-1
      do while (LVEC < JVEC2)
        LVEC = LVEC+1
        LRED = INFVEC(IVEC1+LVEC-1,2,ISYM)
        if (LRED == JRED) then
          LNUM = LNUM+1
        else
          LVEC = JVEC2
        end if
      end do

      if (LNUM > 0) then

        if (JRED /= IREDC) then
          ILOC = 3
          call CHO_GETRED(JRED,ILOC,.false.)
          call CHO_SETREDIND(ILOC)
          IREDC = JRED
        end if

        if (JRED /= IMAPC) then
          call CHO_RS2RS(ISCR,size(ISCR),2,3,JRED,ISYM)
          IMAPC = JRED
        end if

        do LVEC=1,LNUM
          KVEC = KVEC1+LVEC-1
          do IAB=1,NNBSTR(ISYM,2)
            KOFF1 = KCHO1+NNBSTR(ISYM,2)*(KVEC-1)+IAB-1
            KOFF2 = IOFF(min(ISCR(IAB),1))+ISCR(IAB)
            WRK(KOFF1) = WRK(KOFF2)
          end do
          IOFF(1) = IOFF(1)+NNBSTR(ISYM,3)
        end do

        LVEC1 = LVEC1+LNUM
        KVEC1 = KVEC1+LNUM

      end if

    end do

    call CWTIME(C2,W2)
    TIMLOC(1,2) = TIMLOC(1,2)+C2-C1
    TIMLOC(2,2) = TIMLOC(2,2)+W2-W1

    ! Screened or unscreened subtraction section.
    ! The screened version uses level 2 blas, while the unscreened
    ! one employs level 3 blas.
    ! ------------------------------------------------------------

    call CWTIME(C1,W1)

    if (CHO_SSCREEN) then ! screened subtraction

      ! Copy out sub-blocks corresponding to qualified diagonals:
      ! L(#J,{ab})
      ! ---------------------------------------------------------

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

      call CHO_SUBSCR_DIA(WRK(KCHO1),NUMV,ISYM,2,SSNORM)
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

        ! If the qualified block, L({ab},#J), is already in
        ! core, use this block.
        ! -------------------------------------------------

        call DGEMM_('N','T',NNBSTR(ISYM,2),NQUAL(ISYM),NUMV,-One,WRK(KCHO1),NNBSTR(ISYM,2),LQ(ISYM)%A(:,IVEC1_1), &
                    size(LQ(ISYM)%A,1),One,XINT,NNBSTR(ISYM,2))

      else

        ! Copy out sub-blocks corresponding to qualified
        ! diagonals: L({ab},#J)
        ! ----------------------------------------------

        KOFB0 = KCHO1-1-IIBSTR(ISYM,2)
        do J=1,NUMV
          KOFFA = KCHO2+NQUAL(ISYM)*(J-1)-1
          KOFFB = KOFB0+NNBSTR(ISYM,2)*(J-1)
          do IAB=1,NQUAL(ISYM)
            WRK(KOFFA+IAB) = WRK(KOFFB+IQUAB(IAB,ISYM))
          end do
        end do

        ! Subtract:
        ! (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L({ab},#J)
        ! ----------------------------------------------------

        call DGEMM_('N','T',NNBSTR(ISYM,2),NQUAL(ISYM),NUMV,-One,WRK(KCHO1),NNBSTR(ISYM,2),WRK(KCHO2),NQUAL(ISYM),One,XINT, &
                    NNBSTR(ISYM,2))

      end if

    end if

    call CWTIME(C2,W2)
    TIMLOC(1,3) = TIMLOC(1,3)+C2-C1
    TIMLOC(2,3) = TIMLOC(2,3)+W2-W1

  end do

  ! Update counter.
  ! ---------------

  IVEC1 = IVEC1+NVRD

end do

! Update global statistics info.
! ------------------------------

NSYS_CALL = NSYS_CALL+NUMRD
NDGM_CALL = NDGM_CALL+NUMBAT
TDECOM(1,2) = TDECOM(1,2)+TIMLOC(1,1)+TIMLOC(1,2)
TDECOM(2,2) = TDECOM(2,2)+TIMLOC(2,1)+TIMLOC(2,2)
TDECOM(1,3) = TDECOM(1,3)+TIMLOC(1,3)
TDECOM(2,3) = TDECOM(2,3)+TIMLOC(2,3)
if (CHO_SSCREEN) then
  SUBSCRSTAT(1) = SUBSCRSTAT(1)+XTOT
  SUBSCRSTAT(2) = SUBSCRSTAT(2)+XDON
end if

! Print statistics.
! -----------------

if (LOCDBG .or. (IPRINT >= INF_SUBTR1)) then
  if (NUMRD == 0) then
    XAVERD = -9.99999e5_wp
  else
    XAVERD = real(NUMCHO(ISYM),kind=wp)/real(NUMRD,kind=wp)
  end if
  if (NUMBAT == 0) then
    XAVEVC = -9.99999e5_wp
  else
    XAVEVC = real(NUMCHO(ISYM),kind=wp)/real(NUMBAT,kind=wp)
  end if
  write(LUPRI,'(A)') '*****'
  write(LUPRI,'(A,A,I2,A)') SECNAM,' statistics, symmetry',ISYM,':'
  write(LUPRI,'(A,I12)') 'Number of previous vectors                           : ',NUMCHO(ISYM)
  write(LUPRI,'(A,I12)') 'Number of vectors in buffer                          : ',NVEC_IN_BUF(ISYM)
  write(LUPRI,'(A,I12)') 'Memory available for subtraction of previous vectors : ',LWRK
  write(LUPRI,'(A,I12)') 'Memory reserved for buffered vector read             : ',LREAD
  write(LUPRI,'(A,I12)') 'Number of batches needed for reading vectors         : ',NUMRD
  if (CHO_SSCREEN) then
    write(LUPRI,'(A,F12.2)') 'Number of calls to DGEMV                             : ',XDON
    write(LUPRI,'(A,1P,D12.2)') 'Screening threshold                                  : ',SSTAU
    if (XTOT > Zero) then
      SCRPCT = 1.0e2_wp*(XTOT-XDON)/XTOT
    else
      SCRPCT = 1.0e15_wp
    end if
    write(LUPRI,'(A,F12.2,A)') 'Screening percent                                    : ',SCRPCT,'%'
  else
    write(LUPRI,'(A,I12)') 'Number of calls to DGEMM                             : ',NUMBAT
  end if
  write(LUPRI,'(A,I12,I12,F12.2)') 'Minimum, maximum, and average #vectors read          : ',IVSTAT(1,1),IVSTAT(2,1),XAVERD
  write(LUPRI,'(A,I12,I12,F12.2)') 'Minimum, maximum, and average #vecs per call to BLAS : ',IVSTAT(1,2),IVSTAT(2,2),XAVEVC
  write(LUPRI,'(A,2F12.2)') 'Time for reading vectors into buffer (CPU/Wall; sec.): ',TIMLOC(1,1),TIMLOC(2,1)
  write(LUPRI,'(A,2F12.2)') 'Time for reduced set vector reorder  (CPU/Wall; sec.): ',TIMLOC(1,2),TIMLOC(2,2)
  write(LUPRI,'(A,2F12.2)') 'Time for qual. copy + subtraction    (CPU/Wall; sec.): ',TIMLOC(1,3),TIMLOC(2,3)
  write(LUPRI,'(A)') '*****'
end if

end subroutine CHO_SUBTR1
