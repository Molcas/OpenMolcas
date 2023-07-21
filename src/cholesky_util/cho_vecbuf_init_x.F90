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

subroutine Cho_VecBuf_Init_X(Frac,LocDbg)
!
! Purpose: allocate and initialize vector buffer.
!          (External run mode.)

use ChoVecBuf, only: CHVBUF, ip_CHVBUF_SYM, l_CHVBUF_SYM
use stdalloc

implicit none
real*8 Frac
logical LocDbg
#include "cholesky.fh"
character(len=17), parameter :: SecNam = 'Cho_VecBuf_Init_X'
integer, external :: Cho_iSumElm
logical DoRead
integer i, iSym, l_Max, Left, jNum, iRedC, mUsed
integer, parameter :: lScr = 1
real*8 Scr(lScr)
integer nErr
real*8 Diff
real*8, parameter :: Scr_Check = 1.23456789d0, Tol = 1.0d-15
character*2 Unt
real*8 Byte
integer :: l_ChVBuf = 0

if (LocDbg) then
  do i=1,lScr
    Scr(i) = Scr_Check
  end do
  write(Lupri,*) '>>>>> Enter ',SecNam,' <<<<<'
  write(Lupri,*) 'Memory fraction requested for buffer: ',Frac
  write(Lupri,'(A,I8)') 'nSym: ',nSym
  call Cho_Flush(Lupri)
end if

if ((nSym < 1) .or. (nSym > 8)) call Cho_Quit('nSym out of bounds in '//SecNam,102)

if ((Frac <= 0.0d0) .or. (Frac > 1.0d0)) then
  call iZero(l_ChvBuf_Sym,nSym)
  call iZero(ip_ChvBuf_Sym,nSym)
else
  call mma_maxDBLE(l_max)
  Left = int(Frac*dble(l_Max))
  iRedC = -1
  DoRead = .false.
  do iSym=1,nSym
    jNum = 0
    mUsed = 0
    call Cho_VecRd1(Scr,Left,1,NumCho(iSym),iSym,jNum,iRedC,mUsed,DoRead)
    Left = Left-mUsed
    l_ChVBuf_Sym(iSym) = mUsed
  end do
  l_ChVBuf = Cho_iSumElm(l_ChVBuf_Sym,nSym)
  if (l_ChVBuf < 1) then
    l_ChVBuf = 0
    call iZero(l_ChvBuf_Sym,nSym)
    call iZero(ip_ChvBuf_Sym,nSym)
  else
    call mma_allocate(CHVBUF,l_ChVBuf,Label='CHVBUF')

    ip_ChVBuf_Sym(1) = 1
    do iSym=2,nSym
      ip_ChVBuf_Sym(iSym) = ip_ChVBuf_Sym(iSym-1)+l_ChVBuf_Sym(iSym-1)
    end do
  end if
end if

if (LocDbg) then
  nErr = 0
  do i=1,lScr
    Diff = Scr(i)-Scr_Check
    if (abs(Diff) > Tol) nErr = nErr+1
  end do
  if (nErr /= 0) call Cho_Quit('Memory boundary error in '//SecNam,101)
  call Cho_Word2Byte(l_ChVBuf,8,Byte,Unt)
  write(Lupri,*) 'Memory allocated for buffer: ',l_ChVBuf,'(',Byte,Unt,')  at ',1
  write(Lupri,'(A,8I8)') 'l_ChVBuf_Sym : ',(l_ChVBuf_Sym(iSym),iSym=1,nSym)
  write(Lupri,'(A,8I8)') 'ip_ChVBuf_Sym: ',(ip_ChVBuf_Sym(iSym),iSym=1,nSym)
  write(Lupri,*) '>>>>> Exit  ',SecNam,' <<<<<'
  call Cho_Flush(Lupri)
end if

end subroutine Cho_VecBuf_Init_X
