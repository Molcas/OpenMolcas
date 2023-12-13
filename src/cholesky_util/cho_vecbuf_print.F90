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

subroutine Cho_VecBuf_Print(Lupri,nSym)
!
! Purpose: print allocation information of Cholesky vector buffer to
!          unit Lupri (if Lupri<1 nothing is printed).

use Cholesky, only: CHVBUF, l_CHVBUF_SYM
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lupri, nSym
integer(kind=iwp) :: iSym
real(kind=wp) :: xGb
character(len=2) :: Unt
character(len=*), parameter :: SecNam = 'Cho_VecBuf_Print'

if (Lupri < 1) return
if ((nSym < 1) .or. (nSym > 8)) call Cho_Quit('nSym error in '//SecNam,104)

call Cho_Head('Size of Cholesky vector buffer','-',80,Lupri)
write(Lupri,*)
do iSym=1,nSym
  call Cho_Word2Byte(l_ChVBuf_Sym(iSym),8,xGb,Unt)
  write(Lupri,'(A,I2,A,I10,A,F8.2,A,A,A)') 'Dimension, sym.',iSym,': ',l_ChVBuf_Sym(iSym),' 8-byte words (',xGb,' ',Unt,')'
end do
call Cho_Word2Byte(size(ChVBuf),8,xGb,Unt)
write(Lupri,'(/,A,I10,A,F8.2,A,A,A)') 'Total dimension  : ',size(ChVBuf),' 8-byte words (',xGb,' ',Unt,')'

end subroutine Cho_VecBuf_Print
