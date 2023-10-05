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
! Copyright (C) 2016, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_VecBuf_Ini2()
!
! Thomas Bondo Pedersen, June 2006.
!
! Purpose: read vectors from disk into buffer.

use Cholesky, only: CHVBUF, ip_CHVBUF_SYM, l_CHVBUF_SYM, LuPri, nSym, NumCho, NumChT, nVec_in_Buf
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: irc, iRedC, iSym, iV1, iV2, mUsed(8), nRead
logical(kind=iwp) :: DoRead
#ifdef _CHO_DEBUGPRINT_
#define _DEBUGPRINT_
#endif
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: LocDbg = _DBG_
character(len=*), parameter :: SecNam = 'Cho_VecBuf_Ini2'

! Check if buffer is allocated.
! Check if there are any vectors.
! -------------------------------

if (.not. allocated(CHVBUF)) then
  if (LocDbg) write(Lupri,*) SecNam,': returning immediately: No buffer allocated!'
  return
end if
if (NumChT < 1) then
  write(Lupri,*) SecNam,': returning immediately: Buffer allocated, but no vectors!?!?'
  return
end if

! Read vectors.
! -------------

DoRead = .true.
iRedC = -1
do iSym=1,nSym
  iV1 = 1
  iV2 = NumCho(iSym)
  nRead = 0
  mUsed(iSym) = 0
  call Cho_VecRd1(CHVBUF(ip_ChVBuf_Sym(iSym)),l_ChVBuf_Sym(iSym),iV1,iV2,iSym,nRead,iRedC,mUsed(iSym),DoRead)
  nVec_in_Buf(iSym) = nRead
end do

! Debug:
! Enable integrity checks.
! Print info.
! ------------------------

if (LocDbg) then
  call Cho_VecBuf_EnableIntegrityCheck(irc)
  if (irc /= 0) then
    write(LuPri,'(A,I9)') SecNam,': Cho_VecBuf_EnableIntegrityCheck returned code',irc
    call Cho_Quit(SecNam//': integrity check init failed',104)
  else
    write(LuPri,'(A,A)') SecNam,': buffer integrity check enabled'
  end if
  write(Lupri,'(A,A,8I10)') SecNam,'(exit): NumCho:',(NumCho(iSym),iSym=1,nSym)
  write(Lupri,'(A,A,8I10)') SecNam,'(exit): nVec_in_Buf:',(nVec_in_Buf(iSym),iSym=1,nSym)
  write(Lupri,'(A,A,8I10)') SecNam,'(exit): buffer allocated:',(l_ChVBuf_Sym(iSym),iSym=1,nSym)
  write(Lupri,'(A,A,8I10)') SecNam,'(exit): memory used:',(mUsed(iSym),iSym=1,nSym)
  call XFlush(Lupri)
end if

end subroutine Cho_VecBuf_Ini2
