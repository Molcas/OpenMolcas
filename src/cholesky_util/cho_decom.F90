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

subroutine CHO_DECOM(DIAG,WRK,LWRK,IPASS,NUM)
!
! Purpose: calculate Cholesky vectors from qualified integral
!          columns (from disk).

use Cholesky, only: Cho_1Center, Cho_DiaChk, Cho_Simp, DiaMin, iiBstR, IndRed, INF_PASS, INF_PROGRESS, IPRINT, iQuAB, LuPri, &
                    LuSel, nnBstR, nnZTot, nQual, nSym, NumCho, NumChT, nVec_in_Buf, Span, TDECOM, ThrCom, Tol_DiaChk
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LWRK, IPASS, NUM
real(kind=wp), intent(inout) :: Diag(*), WRK(LWRK)
integer(kind=iwp) :: I, IAB, IABG, IADR, ICHO, IDUMP, II, IOPT, ISYM, IVEC, IVEC1, IVECT, JJ, KAB, KCHO1, KEND0, KEND1, KINT1, &
                     KOFF, KOFF0, KOFF2, KOFF3, LENLIN, LINT1, LTOT, LWRK0, LWRK1, NCONV, NERR, NNEG, NNEGT, NUMBUF, NUMCHO_OLD(8)
real(kind=wp) :: C1, C2, FAC, OLDIAG, TOL, W1, W2, XC, XM, XMAX, XMIN, YM
logical(kind=iwp) :: LAST
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_DECOM'

LENLIN = 0  ! to avoid compiler warnings...
if (IPRINT >= INF_PROGRESS) then
  call CHO_HEAD(SECNAM//': Decomposition of Qualified Diagonals','=',80,LUPRI)
  write(LUPRI,'(/,A,I5,A,I4,A)') 'Integral pass number',IPASS,' (',NUM,' shell pair distributions calculated)'
  write(LUPRI,'(A,8I8)') '#Cholesky vec.: ',(NUMCHO(ISYM),ISYM=1,NSYM)
  write(LUPRI,'(A,8I8)') '#vec. in buff.: ',(NVEC_IN_BUF(ISYM),ISYM=1,NSYM)
  write(LUPRI,'(A,8I8)') '#qualified    : ',(NQUAL(ISYM),ISYM=1,NSYM)
  write(LUPRI,'(A,8I8)') 'Current  dim. : ',(NNBSTR(ISYM,2),ISYM=1,NSYM)
  write(LUPRI,'(A,8I8)') 'Original dim. : ',(NNBSTR(ISYM,1),ISYM=1,NSYM)
  write(LUPRI,'(/,A,/,A)') '           #Vectors             Treated Diagonal', &
                           'Sym.     Sym.     Total     Index     Before      After   Conv. Neg.   New Max'
  LENLIN = 79
  write(LUPRI,'(80A)') ('-',I=1,LENLIN)
  call XFLUSH(LUPRI)
  NUMCHO_OLD(1:NSYM) = NUMCHO(1:NSYM)
else if (IPRINT >= INF_PASS) then
  write(LUPRI,'(/,A,I4)') 'Number of shell pair distributions calculated:',NUM
  write(LUPRI,'(A,8I8)') '#Cholesky vec.: ',(NUMCHO(ISYM),ISYM=1,NSYM)
  write(LUPRI,'(A,8I8)') '#vec. in buff.: ',(NVEC_IN_BUF(ISYM),ISYM=1,NSYM)
  write(LUPRI,'(A,8I8)') '#qualified    : ',(NQUAL(ISYM),ISYM=1,NSYM)
  call XFLUSH(LUPRI)
  NUMCHO_OLD(1:NSYM) = NUMCHO(1:NSYM)
end if

! Decompose each symmetry block.
! ------------------------------

do ISYM=1,NSYM

  ! Cycle loop if nothing to do in this symmetry.
  ! ---------------------------------------------

  if ((NQUAL(ISYM) >= 1) .and. (NNBSTR(ISYM,2) >= 1)) then

    ! Reserve space for qualified integral columns.
    ! ---------------------------------------------

    LINT1 = NNBSTR(ISYM,2)*NQUAL(ISYM) ! integrals

    KINT1 = 1
    KEND0 = KINT1+LINT1
    LWRK0 = LWRK-KEND0+1
    if (LWRK0 <= 0) call CHO_QUIT('[0] Insufficient memory in '//SECNAM,101)

    ! Determine size of Cholesky vector (output) buffer.
    ! --------------------------------------------------

    NUMBUF = min(LWRK0/NNBSTR(ISYM,2),NQUAL(ISYM))
    if (NUMBUF < 1) call CHO_QUIT('[1] Insufficient memory in '//SECNAM,101)
    KCHO1 = KEND0
    KEND1 = KCHO1+NNBSTR(ISYM,2)*NUMBUF
    LWRK1 = LWRK-KEND1+1
    if (LWRK1 < 0) call CHO_QUIT('Buffer allocation error in '//SECNAM,101)  ! should be redundant...

    ! Read qualified integral columns.
    ! --------------------------------

    call CWTIME(C1,W1)
    IOPT = 2
    LTOT = NNBSTR(ISYM,2)*NQUAL(ISYM)
    IADR = 0
    call DDAFILE(LUSEL(ISYM),IOPT,WRK(KINT1),LTOT,IADR)
    call CWTIME(C2,W2)
    TDECOM(1,1) = TDECOM(1,1)+C2-C1
    TDECOM(2,1) = TDECOM(2,1)+W2-W1

    ! Subtract contributions from previous vectors.
    ! ---------------------------------------------

    call CHO_SUBTR(WRK(KINT1),WRK(KEND0),LWRK0,ISYM)

    ! Debug: check diagonal elements in updated integrals.
    ! ----------------------------------------------------

    if (CHO_DIACHK .or. LOCDBG) then
      TOL = TOL_DIACHK
      NERR = 0
      call CHO_CHKINT(WRK(KINT1),DIAG,ISYM,NERR,TOL,.true.)
      if (NERR /= 0) then
        write(LUPRI,*) SECNAM,': ',NERR,' diagonal errors found!'
        write(LUPRI,*) '          #tests: ',NQUAL(ISYM)
        !write(LUPRI,*) '          Printing integrals:'
        !call CHO_OUTPUT(WRK(KINT1),1,NNBSTR(ISYM,2),1,NQUAL(ISYM),NNBSTR(ISYM,2),NQUAL(ISYM),1,LUPRI)
        call CHO_QUIT('Diagonal errors in '//SECNAM,104)
      else
        write(LUPRI,*) SECNAM,': comparison of qual. integrals and current diagonal: no errors !'
      end if
    end if

    ! Decompose in loop over qualified columns.
    ! -----------------------------------------

    IVEC = NUMCHO(ISYM)
    IVECT = NUMCHT
    IDUMP = 0
    do ICHO=1,NQUAL(ISYM)

      ! Find max. diagonal among qualified.
      ! -----------------------------------

      IAB = 1
      IABG = INDRED(IQUAB(IAB,ISYM),2)
      XC = DIAG(IABG)
      do I=2,NQUAL(ISYM)
        KAB = INDRED(IQUAB(I,ISYM),2)
        if (DIAG(KAB) > XC) then
          IAB = I
          IABG = KAB
          XC = DIAG(KAB)
        end if
      end do

      ! Decompose if max. diagonal is still qualified.
      ! ----------------------------------------------

      LAST = (XC < DIAMIN(ISYM)) .or. (XC < THRCOM)
      if (.not. LAST) then

        ! Offset to max. diagonal column.
        ! -------------------------------

        KOFF0 = KINT1+NNBSTR(ISYM,2)*(IAB-1)-1

        ! Scale column corresponding to max. diagonal to obtain
        ! the Cholesky vector.
        ! -----------------------------------------------------

        FAC = One/sqrt(XC)
        KOFF = KOFF0
        WRK(KOFF+1:KOFF+NNBSTR(ISYM,2)) = FAC*WRK(KOFF+1:KOFF+NNBSTR(ISYM,2))

        ! Zero entries in Cholesky vector corresponding to zero diagonals.
        ! ----------------------------------------------------------------

        do I=1,NNBSTR(ISYM,2)
          II = IIBSTR(ISYM,2)+I
          JJ = INDRED(II,2)
          if (DIAG(JJ) == Zero) then
            KOFF = KOFF0+I
            WRK(KOFF) = Zero
          end if
        end do

        ! Update diagonal.
        ! ----------------

        do I=1,NNBSTR(ISYM,2)
          II = IIBSTR(ISYM,2)+I
          JJ = INDRED(II,2)
          KOFF = KOFF0+I
          DIAG(JJ) = DIAG(JJ)-WRK(KOFF)*WRK(KOFF)
        end do

        ! Zero treated diagonal element and analyze updated diagonal.
        ! -----------------------------------------------------------

        OLDIAG = DIAG(IABG)
        DIAG(IABG) = Zero
        call CHO_CHKDIA(DIAG,ISYM,XMIN,XMAX,XM,NNEGT,NNEG,NCONV)

        ! Update total number of zeroed negative diagonals.
        ! -------------------------------------------------

        NNZTOT = NNZTOT+NNEG

        ! Update DIAMIN from max. abs. diagonal element XM.
        ! CHO_1CENTER: update from max. diagonal element among
        !              qualified.
        ! CHO_SIMP   : "simulate parallel algorithm" = do not
        !              update DIAMIN.
        ! ----------------------------------------------------

        if (.not. CHO_SIMP) then
          if (CHO_1CENTER) then
            YM = DIAG(INDRED(IQUAB(1,ISYM),2))
            do I=2,NQUAL(ISYM)
              YM = max(YM,DIAG(INDRED(IQUAB(I,ISYM),2)))
            end do
          else
            YM = XM
          end if
          DIAMIN(ISYM) = max(YM*SPAN,THRCOM)
        end if

        ! Subtract this Cholesky vector from integrals. If
        ! the corresponding diagonal element is zero, the
        ! column will no longer be qualified and subtraction
        ! can safely be skipped.
        ! --------------------------------------------------

        do I=1,NQUAL(ISYM)
          II = IQUAB(I,ISYM)
          JJ = INDRED(II,2)
          if (DIAG(JJ) /= Zero) then
            KOFF2 = KINT1+NNBSTR(ISYM,2)*(I-1)-1
            KOFF3 = KOFF0+II-IIBSTR(ISYM,2)
            FAC = -WRK(KOFF3)
            WRK(KOFF2+1:KOFF2+NNBSTR(ISYM,2)) = WRK(KOFF2+1:KOFF2+NNBSTR(ISYM,2))-WRK(KOFF3)*WRK(KOFF0+1:KOFF0+NNBSTR(ISYM,2))
          end if
        end do

        ! Store Cholesky vector in buffer.
        ! --------------------------------

        IDUMP = IDUMP+1

        KOFF2 = KCHO1+NNBSTR(ISYM,2)*(IDUMP-1)-1
        WRK(KOFF2+1:KOFF2+NNBSTR(ISYM,2)) = WRK(KOFF0+1:KOFF0+NNBSTR(ISYM,2))

        ! Update Cholesky vector counters.
        ! --------------------------------

        IVEC = IVEC+1
        IVECT = IVECT+1

        ! Set info for this vector.
        ! -------------------------

        call CHO_SETVECINF(IVEC,ISYM,IABG,IPASS,2)

        ! Print progress report.
        ! ----------------------

        if (IPRINT >= INF_PROGRESS) write(LUPRI,'(I3,3(1X,I9),2(1X,ES11.3),2(1X,I4),1X,ES11.3)') ISYM,IVEC,IVECT,IABG,XC,OLDIAG, &
                                                                                                 NCONV,NNEG,XM

      end if

      ! Dump vectors to disk when there is no more to be done, or
      ! when the buffer is full.
      ! ---------------------------------------------------------

      if (LAST .or. (IDUMP == NUMBUF)) then
        call CWTIME(C1,W1)
        IVEC1 = NUMCHO(ISYM)+1
        call CHO_PUTVEC(WRK(KCHO1),NNBSTR(ISYM,2),IDUMP,IVEC1,ISYM)
        call CHO_VECBUF_COPY(WRK(KCHO1),IDUMP,ISYM)
        NUMCHO(ISYM) = NUMCHO(ISYM)+IDUMP
        NUMCHT = NUMCHT+IDUMP
        call CWTIME(C2,W2)
        TDECOM(1,2) = TDECOM(1,2)+C2-C1
        TDECOM(2,2) = TDECOM(2,2)+W2-W1
        if (LAST) then
          exit  ! cycle symmetry loop
        else
          IVEC = NUMCHO(ISYM)
          IVECT = NUMCHT
          IDUMP = 0
        end if
      end if

    end do
  end if

  ! Cycle point: go to next symmetry.
  ! ---------------------------------

  if (IPRINT >= INF_PROGRESS) call XFLUSH(LUPRI)

end do

if (IPRINT >= INF_PROGRESS) then
  NUMCHO_OLD(1:NSYM) = NUMCHO(1:NSYM)-NUMCHO_OLD(1:NSYM)
  write(LUPRI,'(80A)') ('-',I=1,LENLIN)
  write(LUPRI,'(A,8I8)') '#vec. gener.  : ',(NUMCHO_OLD(ISYM),ISYM=1,NSYM)
else if (IPRINT >= INF_PASS) then
  NUMCHO_OLD(1:NSYM) = NUMCHO(1:NSYM)-NUMCHO_OLD(1:NSYM)
  write(LUPRI,'(A,8I8)') '#vec. gener.  : ',(NUMCHO_OLD(ISYM),ISYM=1,NSYM)
end if

end subroutine CHO_DECOM
