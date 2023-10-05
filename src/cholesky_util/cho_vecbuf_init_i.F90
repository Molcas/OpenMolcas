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

subroutine Cho_VecBuf_Init_I(Frac,lVec,LocDbg)
!
! Purpose: allocate and initialize vector buffer.
!          (Internal run mode.)

use Cholesky, only: CHVBUF, ip_CHVBUF_SYM, l_CHVBUF_SYM, LuPri, MaxVec, nSym, nVec_in_Buf
use stdalloc, only: mma_allocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Frac
integer(kind=iwp), intent(in) :: lVec(*)
logical(kind=iwp), intent(in) :: LocDbg
integer(kind=iwp) :: iSym, l_ChVBuf = 0, l_Max, lVecTot, MemEach, MemLeft
real(kind=wp) :: x, xMemMax(8)
logical(kind=iwp) :: Enough
character(len=2) :: Unt
character(len=*), parameter :: SecNam = 'Cho_VecBuf_Init_I'

if (LocDbg) then
  write(Lupri,*) '>>>>> Enter ',SecNam,' <<<<<'
  write(Lupri,*) 'Memory fraction requested for buffer: ',Frac
  write(Lupri,'(A,I8)') 'nSym: ',nSym
  write(Lupri,'(A,8I8)') 'lVec: ',(lVec(iSym),iSym=1,nSym)
  call XFlush(Lupri)
end if

if ((nSym < 1) .or. (nSym > 8)) call Cho_Quit('nSym out of bounds in '//SecNam,102)

x = real(MaxVec,kind=wp)
lVecTot = lVec(1)
xMemMax(1) = real(lVec(1),kind=wp)*x
do iSym=2,nSym
  lVecTot = max(lVecTot,lVec(iSym))
  xMemMax(iSym) = real(lVec(iSym),kind=wp)*x
end do

if ((Frac <= Zero) .or. (Frac > One) .or. (lVecTot < 1)) then
  ip_ChVBuf_Sym(1:nSym) = 0
  l_ChVBuf_Sym(1:nSym) = 0
else
  call mma_MaxDBLE(l_Max)
  l_ChVBuf = int(Frac*real(l_Max,kind=wp))
  if ((l_ChVBuf < nSym) .or. (l_ChVBuf < lVecTot)) then
    l_ChVBuf = 0
    ip_ChVBuf_Sym(1:nSym) = 0
    l_ChVBuf_Sym(1:nSym) = 0
  else
    MemEach = l_ChVBuf/nSym
    Enough = all(MemEach > lVec(1:nSym))
    if (.not. Enough) then ! whole buffer for sym. 1
      l_ChVBuf_Sym(1) = l_ChVBuf
      l_ChVBuf_Sym(2:nSym) = 0
    else
      MemLeft = l_ChVBuf-nSym*MemEach
      l_ChVBuf_Sym(1) = MemEach+MemLeft
      if (real(l_ChVBuf_Sym(1),kind=wp) > xMemMax(1)) l_ChVBuf_Sym(1) = int(xMemMax(1))
      do iSym=2,nSym
        l_ChVBuf_Sym(iSym) = MemEach
        if (real(l_ChVBuf_Sym(iSym),kind=wp) > xMemMax(iSym)) l_ChVBuf_Sym(iSym) = int(xMemMax(iSym))
      end do
    end if
    l_ChVBuf = sum(l_ChVBuf_Sym(1:nSym))

    call mma_allocate(CHVBUF,l_ChVBuf,Label='CHVBUF')

    ip_ChVBuf_Sym(1) = 1
    do iSym=2,nSym
      ip_ChVBuf_Sym(iSym) = ip_ChVBuf_Sym(iSym-1)+l_ChVBuf_Sym(iSym-1)
    end do
  end if
end if

nVec_in_Buf(1:nSym) = 0

if (LocDbg) then
  call Cho_Word2Byte(l_ChVBuf,8,x,Unt)
  write(Lupri,*) 'Memory allocated for buffer: ',l_ChVBuf,'(',x,Unt,') at ',1
  write(Lupri,'(A,8I8)') 'l_ChVBuf_Sym : ',(l_ChVBuf_Sym(iSym),iSym=1,nSym)
  write(Lupri,'(A,8I8)') 'ip_ChVBuf_Sym: ',(ip_ChVBuf_Sym(iSym),iSym=1,nSym)
  write(Lupri,*) '>>>>> Exit  ',SecNam,' <<<<<'
  call XFlush(Lupri)
end if

end subroutine Cho_VecBuf_Init_I
