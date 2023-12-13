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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_VecBuf_EnableIntegrityCheck(irc)
!
! Thomas Bondo Pedersen, September 2012.
!
! Enable integrity check of buffer: allocate and store norm and sum
! of each vector in the buffer.

use Cholesky, only: CHVBFI, CHVBUF, InfVec, ip_CHVBFI_SYM, ip_CHVBUF_SYM, IPRINT, l_CHVBFI_SYM, LuPri, nDimRS, nSym, nVec_in_Buf, &
                    RUN_EXTERNAL, RUN_MODE
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), external :: dDot_
integer(kind=iwp) :: ip, ipV, iSym, jRed, jVec, l_ChVBfI

! Set return code
irc = 0

! Not implemented in internal run mode (since the buffered vectors
! may change length and hence norm and sum, causing too much book
! keeping activity for a debug feature).
if (RUN_MODE /= RUN_EXTERNAL) return

! Return if no buffer is allocated
if (.not. allocated(CHVBUF)) return

! Return if already enabled
if (allocated(CHVBFI)) return

! Check that nDimRS is allocated
if (.not. allocated(nDimRS)) then
  irc = 1
  return
end if

! Allocate and store norm and sum of each vector in the buffer
l_ChVBfI_Sym(1:nSym) = nVec_in_Buf(1:nSym)
l_ChVBfI = sum(l_ChVBfI_Sym(1:nSym))
if (l_ChVBfI > 0) then
  call mma_allocate(CHVBFI,2,l_ChVBfI,Label='CHVBFI')
  ip = 1
  do iSym=1,nSym
    ip_ChVBfI_Sym(iSym) = ip
    ip = ip+l_ChVBfI_Sym(iSym)
  end do
  do iSym=1,nSym
    ipV = ip_ChVBuf_Sym(iSym)
    ip = ip_ChvBfI_Sym(iSym)
    do jVec=1,nVec_in_Buf(iSym)
      jRed = InfVec(jVec,2,iSym)
      CHVBFI(1,ip) = sqrt(dDot_(nDimRS(iSym,jRed),CHVBUF(ipV),1,CHVBUF(ipV),1))
      CHVBFI(2,ip) = sum(CHVBUF(ipV:ipV+nDimRS(iSym,jRed)-1))
      ipV = ipV+nDimRS(iSym,jRed)
      ip = ip+1
    end do
  end do
  if (iPrint > 2) call Cho_VecBuf_PrtRef('@NABLE')
  write(LuPri,'(A)') 'Cholesky vector buffer integrity checks enabled'
else
  l_ChVBfI_Sym(1:nSym) = 0
  ip_ChVBfI_Sym(1:nSym) = 0
end if

end subroutine Cho_VecBuf_EnableIntegrityCheck
