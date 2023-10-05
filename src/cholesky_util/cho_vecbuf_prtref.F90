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

subroutine Cho_VecBuf_PrtRef(Txt)
!
! Thomas Bondo Pedersen, September 2012.
!
! Print reference norm and sum of vectors in buffer.
! Txt is printed along with the reference values (for
! identification).

use Cholesky, only: CHVBFI, InfVec, ip_CHVBFI_SYM, LuPri, nDimRS, nSym, nVec_in_Buf
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Txt
integer(kind=iwp) :: iSym, jRed, jVec, nDim

if (.not. allocated(nDimRS)) call Cho_Quit('Cho_VecBuf_PrtRef: unable to print reference values',104)
if (allocated(CHVBFI)) then
  do iSym=1,nSym
    do jVec=1,nVec_in_Buf(iSym)
      jRed = InfVec(jVec,2,iSym)
      nDim = nDimRS(iSym,jRed)
      write(LuPri,'(A,A,I6,A,I2,A,I9,1P,2(A,D25.16))') Txt,' Cholesky vector',jVec,' sym.',iSym,' dim.',nDim,'  Norm=', &
                                                       CHVBFI(1,ip_ChVBfI_Sym(iSym)+jVec),' Sum=',CHVBFI(2,ip_ChVBfI_Sym(iSym)+jVec)
    end do
  end do
else
  write(LuPri,'(A,A)') Txt,' Cho_VecBuf_PrtRef: no reference values available!'
end if

end subroutine Cho_VecBuf_PrtRef
