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

subroutine CHO_GETINT(DIAG,DIASH,ISYSH,LSTQSP,NPOTSH,ICOUNT)
!
! Purpose: get qualified integral columns for Cholesky decomposition.
!
! DIASH(ij): max. diagonal in shell pair i,j
! NPOTSH   : the number of shell pairs that can be qualified.

use Cholesky, only: DiaMin, INF_IN2, IntMap, iOffq, IPRINT, iSP2F, LuPri, MinQual, MXSHPR, N1_Qual, N2_Qual, NCOLAB, nnBstR, &
                    nQual, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
real(kind=wp), intent(inout) :: DIASH(*)
integer(kind=iwp), intent(in) :: ISYSH(*), NPOTSH
integer(kind=iwp), intent(out) :: LSTQSP(NPOTSH), ICOUNT
integer(kind=iwp) :: i, ISHLA, ISHLAB, ISHLB, ISYM, ISYMAB, LMAX, MCOUNT, MEMQ(1), MXDIM, NSEL
real(kind=wp) :: SMAX, XMMQ
logical(kind=iwp) :: DODECO, FULL, SYNC
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_GETINT'

!-tbp: some debugging...
!if (LOCDBG) then
!  do LEVEL=1,3
!    call CHO_MCA_INT_1_DBG(DIAG,LEVEL)
!  end do
!  call cho_quit(SECNAM//' end of test',100)
!end if

! Initializations.
! ----------------

NQUAL(1:NSYM) = 0
ICOUNT = 0
if (MXSHPR > 0) then
  MCOUNT = min(NPOTSH,MXSHPR)
else
  MCOUNT = NPOTSH
end if
DODECO = .false.

MXDIM = NNBSTR(1,2)
do ISYM=2,NSYM
  MXDIM = max(MXDIM,NNBSTR(ISYM,2))
end do
call mma_maxDBLE(LMAX)
XMMQ = real(N1_QUAL,kind=wp)*real(LMAX,kind=wp)/real(N2_QUAL,kind=wp)
MEMQ(1) = int(XMMQ)
call CHO_GAIGOP(MEMQ,1,'min')
if (MEMQ(1) < MXDIM) then
  write(LUPRI,*) SECNAM,': memory split error!'
  write(LUPRI,*) 'Memory for storing qualified columns: ',MEMQ(1)
  write(LUPRI,*) 'Minimal memory needed to store one column: ',MXDIM
  write(LUPRI,*) 'Total memory available: ',LMAX
  write(LUPRI,*) 'Memory split is ',N1_QUAL,'/',N2_QUAL,' for qualified columns.'
  write(LUPRI,*) 'Change memory split in input file...'
  call CHO_QUIT('Memory split error in '//SECNAM,101)
end if

! Shell pair qualification loop.
! ------------------------------

do while ((.not. DODECO) .and. (ICOUNT < MCOUNT))

  ! Update shell pair counter.
  ! --------------------------

  ICOUNT = ICOUNT+1

  ! Get shell pair corresponding to largest diagonal.
  ! -------------------------------------------------

  call CHO_P_GETMAXSHL(DIASH,SMAX,ISHLAB)
  call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
  ISYMAB = ISYSH(ISHLAB)

  if ((SMAX == Zero) .or. (abs(SMAX) < DIAMIN(ISYMAB))) then

    ! Diagonal too small to be qualified for decomposition.
    ! -----------------------------------------------------

    if (ICOUNT == 1) then
      write(LUPRI,*) SECNAM,': no integrals calculated; unable to proceed to decomposition!'
      write(LUPRI,*) 'Max. abs. diagonal for shell pair ',ISHLA,', ',ISHLB,': ',abs(SMAX)
      write(LUPRI,*) 'Max. abs. diagonal allowed: ',DIAMIN(ISYMAB),' (sym. ',ISYMAB,')'
      call CHO_QUIT('Severe error in '//SECNAM,104)
    else
      ICOUNT = ICOUNT-1
      NSEL = sum(NQUAL(1:NSYM))
      DODECO = NSEL > 0
    end if

  else

    ! Qualify diagonals within this shell pair.
    ! -----------------------------------------

    SYNC = .false.
    FULL = .false.
    call CHO_P_QUALIFY(DIAG,SYNC,ISHLAB,ISYMAB,MEMQ(1),FULL)

    ! Calculate integral columns; get qualified ones stored in
    ! current reduced set; write these to disk on temporary file(s).
    ! --------------------------------------------------------------

    NSEL = sum(NQUAL(1:NSYM))
    NCOLAB = NSEL-sum(IOFFQ(1:NSYM))

    if (NCOLAB > 0) then

      INTMAP(ISHLAB) = INTMAP(ISHLAB)+1
      if (IPRINT >= INF_IN2) then
        write(LUPRI,'(/,A,I5,1X,I5,A,I9,A)') 'Calculating shell pair (**|',ISHLA,ISHLB,'):',NCOLAB,' columns have been qualified'
        write(LUPRI,'(80A)') ('=',i=1,77)
        write(LUPRI,'(A,I12)') 'Number of calculations so far for this shell pair: ',INTMAP(ISHLAB)
      end if

      LSTQSP(ICOUNT) = ISHLAB
      call CHO_MCA_CALCINT(ISHLAB)

      ! Enough integral columns for proceeding to decomposition?
      ! --------------------------------------------------------

      DODECO = FULL .or. (NSEL >= MINQUAL)

    else if (NCOLAB == 0) then

      if (NSEL < 1) then
        write(LUPRI,*) SECNAM,': logical error: unable to qualify diagonals'
        write(LUPRI,*) SECNAM,': NCOLAB = ',NCOLAB
        write(LUPRI,*) SECNAM,': NSEL   = ',NSEL
        call CHO_QUIT('[0] Logical error in '//SECNAM,104)
      else
        ICOUNT = ICOUNT-1
        DODECO = .true.
      end if

    else

      write(LUPRI,*) SECNAM,': logical error: unable to qualify diagonals'
      write(LUPRI,*) SECNAM,': NCOLAB = ',NCOLAB
      write(LUPRI,*) SECNAM,': NSEL   = ',NSEL
      call CHO_QUIT('[1] Logical error in '//SECNAM,104)

    end if

  end if

end do

! Test loop exit (we may have calculated all possible integrals, yet
! NSEL < MINQUAL or allowed memory may have been used).
! ------------------------------------------------------------------

if (.not. DODECO) then
  NSEL = sum(NQUAL(1:NSYM))
  if (NSEL < 1) then
    write(LUPRI,*) SECNAM,': logical error: unable to qualify diagonals'
    write(LUPRI,*) SECNAM,': Flag DODECO is ',DODECO
    write(LUPRI,*) SECNAM,': NSEL    = ',NSEL
    write(LUPRI,*) SECNAM,': ICOUNT  = ',ICOUNT
    write(LUPRI,*) SECNAM,': MCOUNT  = ',MCOUNT
    write(LUPRI,*) SECNAM,': NPOTSH  = ',NPOTSH
    write(LUPRI,*) SECNAM,': MINQUAL = ',MINQUAL
    call CHO_QUIT('[2] Logical error in '//SECNAM,103)
  else
    DODECO = .true.
  end if
end if

! Set indices for local qualified (parallel runs).
! ------------------------------------------------

call CHO_P_SETLQ()

end subroutine CHO_GETINT
