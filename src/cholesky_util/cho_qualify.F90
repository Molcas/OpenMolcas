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

subroutine CHO_QUALIFY(DIAG,ISHLAB,ISYMAX,MEM,FULL)
!
! Purpose: qualify diagonal elements for decomposition in
!          current reduced set. ISYMAX is the symmetry block
!          to which the largest diagonal belongs.
!          MEM is the (total!) max. allowed
!          memory for storing the qualified columns.
!          If no more columns can be qualified on exit,
!          FULL=.true. is returned.

use Cholesky, only: DiaMin, IALQUA, iiBstR, iiBstRSh, IndRed, iOffq, iQuAB, iSP2F, LuPri, MaxQual, nnBstR, nnBstRSh, nQual, nSym
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
integer(kind=iwp), intent(in) :: ISHLAB, ISYMAX, MEM
logical(kind=iwp), intent(out) :: FULL
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
integer(kind=iwp) :: I, I1, I2, ISHLA, ISHLB, ISYM, K, K1, K2, LEFT, MEM0, MINM, NEED, NUM
logical(kind=iwp), parameter :: LOCDBG = _DBG_
character(len=*), parameter :: SECNAM = 'CHO_QUALIFY'

! Copy counter to offset array.
! -----------------------------

IOFFQ(1:NSYM) = NQUAL(1:NSYM)

! Check memory.
! -------------

MEM0 = sum(NQUAL(1:NSYM)*NNBSTR(1:NSYM,2))
LEFT = MEM-MEM0
if (IALQUA == 0) then
  MINM = NNBSTR(1,2)
  do ISYM=1,NSYM
    MINM = max(MINM,NNBSTR(ISYM,2))
  end do
else
  MINM = NNBSTR(ISYMAX,2)
end if
FULL = LEFT < MINM
if (FULL) return

! Qualify.
! --------

if (IALQUA == 0) then  ! qualify until full (dalton style)
  do ISYM=1,NSYM
    call CHO_QUALIFY_1(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
  end do
else if (IALQUA == 1) then  ! qualify until full
  call CHO_QUALIFY_1(DIAG,ISYMAX,ISHLAB,MEM,MEM0,LEFT)
  do ISYM=1,ISYMAX-1
    call CHO_QUALIFY_1(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
  end do
  do ISYM=ISYMAX+1,NSYM
    call CHO_QUALIFY_1(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
  end do
else ! qualify until full, then largest
  call CHO_QUALIFY_2(DIAG,ISYMAX,ISHLAB,MEM,MEM0,LEFT)
  do ISYM=1,ISYMAX-1
    call CHO_QUALIFY_2(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
  end do
  do ISYM=ISYMAX+1,NSYM
    call CHO_QUALIFY_2(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
  end do
end if

! Set FULL flag:
! FULL=.true. if a) not enough memory to qualify another column of
! any symmetry, or b) MAXQUAL reached in any symmetry.
! ----------------------------------------------------------------

NEED = sum(NQUAL(1:NSYM)*NNBSTR(1:NSYM,2))
if ((NEED < 1) .or. (NEED > MEM)) then
  call CHO_QUIT('Logical error (2) in '//SECNAM,104)
else
  LEFT = MEM-NEED
  FULL = .false.
  ISYM = 0
  do while ((ISYM < NSYM) .and. (.not. FULL))
    ISYM = ISYM+1
    if ((NQUAL(ISYM) < IOFFQ(ISYM)) .or. (NQUAL(ISYM) < 0) .or. (NQUAL(ISYM) > MAXQUAL)) then
      call CHO_QUIT('Logical error (3) in '//SECNAM,104)
    else
      FULL = NQUAL(ISYM) == MAXQUAL
    end if
    if (NNBSTR(ISYM,2) > 0) FULL = FULL .or. (LEFT < NNBSTR(ISYM,2))
  end do
end if

! Debug: print.
! -------------

if (LOCDBG) then
  call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
  write(LUPRI,*)
  write(LUPRI,*)
  write(LUPRI,*) SECNAM,': qualified diagonals from shell-pair ',ISHLA,ISHLB,':'
  write(LUPRI,*) 'Qualification algorithm: ',IALQUA
  write(LUPRI,*) 'Total memory for qualification: ',MEM,'  Memory left: ',LEFT
  do ISYM=1,NSYM
    NUM = NQUAL(ISYM)-IOFFQ(ISYM)
    write(LUPRI,*)
    write(LUPRI,*) 'Sym.,dimension,#qualified,threshold: ',ISYM,NNBSTRSH(ISYM,ISHLAB,2),NUM,DIAMIN(ISYM)
    if (NNBSTRSH(ISYM,ISHLAB,2) > 0) then
      I1 = IIBSTR(ISYM,2)+IIBSTRSH(ISYM,ISHLAB,2)+1
      I2 = I1+NNBSTRSH(ISYM,ISHLAB,2)-1
      write(LUPRI,*) 'Diagonal (current reduced set):'
      write(LUPRI,'(5F15.8)') (DIAG(INDRED(I,2)),I=I1,I2)
      K1 = IOFFQ(ISYM)+1
      K2 = NQUAL(ISYM)
      write(LUPRI,*) 'Qualified diagonals:'
      write(LUPRI,'(5F15.8)') (DIAG(INDRED(IQUAB(K,ISYM),2)),K=K1,K2)
    end if
  end do
  write(LUPRI,*)
  write(LUPRI,*)
end if

end subroutine CHO_QUALIFY
