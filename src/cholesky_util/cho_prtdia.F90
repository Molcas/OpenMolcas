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

subroutine CHO_PRTDIA(DIAG,ISYLST,NSYLST,IRED)
!
! Purpose: print requested symmetry block(s) of diagonal in
!          first (IRED=1) or current (IRED=2) reduced set.

use Cholesky, only: iiBstR, iiBstRSh, IndRed, IndRSh, iSP2F, LuPri, nnBstRSh, nnShl, nSym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSYLST, ISYLST(NSYLST), IRED
real(kind=wp), intent(in) :: Diag(*)
integer(kind=iwp) :: IAB, IAB1, IAB2, ILST, ISHLAB, ISYM, JAB
character(len=*), parameter :: SECNAM = 'CHO_PRTDIA'

! Check dimension of symmetry list.
! ---------------------------------

if (NSYLST < 1) then
  return
else if (NSYLST > NSYM) then
  write(LUPRI,'(//,1X,A,A)') SECNAM,': NSYLST <= NSYM required!'
  write(LUPRI,'(1X,A,I10)') 'NSYLST = ',NSYLST
  write(LUPRI,'(1X,A,I10,/)') 'NSYM   = ',NSYM
  call CHO_QUIT('[0] Symmetry error in '//SECNAM,102)
end if

! Code for first or current reduced set.
! --------------------------------------

if (IRED == 1) then
  call CHO_HEAD(SECNAM//': Diagonal in Original Reduced Set','=',80,LUPRI)
  do ILST=1,NSYLST
    ISYM = ISYLST(ILST)
    if ((ISYM < 1) .or. (ISYM > NSYM)) then
      write(LUPRI,*) SECNAM,': element ',ILST,': ',ISYM,' of list ISYLST is out of bounds!'
      call CHO_QUIT('ISYLST input error in '//SECNAM,104)
    else
      write(LUPRI,'(/,A,I2)') 'Symmetry block:',ISYM
      write(LUPRI,'(/,A,/,A)') '  Element Shell-Pair  SP Index         Diagonal','-----------------------------------------------'
      do ISHLAB=1,NNSHL
        IAB1 = IIBSTR(ISYM,IRED)+IIBSTRSH(ISYM,ISHLAB,IRED)+1
        IAB2 = IAB1+NNBSTRSH(ISYM,ISHLAB,IRED)-1
        do IAB=IAB1,IAB2
          if (INDRSH(IAB) /= ISP2F(ISHLAB)) then
            write(LUPRI,*) 'Shell Pair error: INDRSH,ISP2F,ISHLAB',INDRSH(IAB),ISP2F(ISHLAB),ISHLAB
            call CHO_QUIT('Shell-Pair error in '//SECNAM,104)
          else
            JAB = IAB
            write(LUPRI,'(I9,2X,I9,1X,I9,1X,1P,D16.8)') JAB,ISP2F(ISHLAB),INDRED(IAB,IRED),DIAG(JAB)
          end if
        end do
      end do
      write(LUPRI,'(A)') '-----------------------------------------------'
    end if
  end do
else if (IRED == 2) then
  call CHO_HEAD(SECNAM//': Diagonal in Current Reduced Set','=',80,LUPRI)
  do ILST=1,NSYLST
    ISYM = ISYLST(ILST)
    if ((ISYM < 1) .or. (ISYM > NSYM)) then
      write(LUPRI,*) SECNAM,': element ',ILST,': ',ISYM,' of list ISYLST is out of bounds!'
      call CHO_QUIT('ISYLST input error in '//SECNAM,104)
    else
      write(LUPRI,'(/,A,I2)') 'Symmetry block:',ISYM
      write(LUPRI,'(/,A,/,A)') '  Element  RedSet 1 Shell-Pair  SP Index         Diagonal', &
                               '---------------------------------------------------------'
      do ISHLAB=1,NNSHL
        IAB1 = IIBSTR(ISYM,IRED)+IIBSTRSH(ISYM,ISHLAB,IRED)+1
        IAB2 = IAB1+NNBSTRSH(ISYM,ISHLAB,IRED)-1
        do IAB=IAB1,IAB2
          JAB = INDRED(IAB,IRED)
          if (INDRSH(JAB) /= ISP2F(ISHLAB)) then
            write(LUPRI,*) 'Shell Pair error: INDRSH,ISP2F,ISHLAB',INDRSH(JAB),ISP2F(ISHLAB),ISHLAB
            call CHO_QUIT('Shell-Pair error in '//SECNAM,104)
          else
            write(LUPRI,'(I9,1X,I9,2X,I9,1X,I9,1X,1P,D16.8)') IAB,JAB,ISP2F(ISHLAB),INDRED(JAB,1),DIAG(JAB)
          end if
        end do
      end do
      write(LUPRI,'(A)') '---------------------------------------------------------'
    end if
  end do
end if

end subroutine CHO_PRTDIA
