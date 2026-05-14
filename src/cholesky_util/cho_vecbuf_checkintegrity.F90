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

subroutine Cho_VecBuf_CheckIntegrity(Tol,Verbose,Txt,irc)
!
! Thomas Bondo Pedersen, September 2012.
!
! Check integrity of Cholesky vector buffer and print result.
! Tol is the tolerance to use for comparing norm and sum.
! Txt is printed along with the check message if Verbose=.True.
! (if Verbose=.False. nothing is printed).
! Return code:
!    irc=0: buffer OK
!    irc=1: buffer corrupted
!
! A simpler interface is given by Subroutine Cho_VecBuf_Check.

use Cholesky, only: LuPri
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Tol
logical(kind=iwp), intent(in) :: Verbose
character(len=*), intent(in) :: Txt
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp) :: Cho_VecBuf_Integrity_OK

if (Cho_VecBuf_Integrity_OK(Tol,Verbose)) then
  if (Verbose) then
    write(LuPri,'(A,A)') Txt,' Cholesky vector buffer integrity checked: OK'
    call XFlush(LuPri)
  end if
  irc = 0
else
  if (Verbose) then
    write(LuPri,'(A,A)') Txt,' Cholesky vector buffer integrity checked: CORRUPTED'
    call Cho_Quit('Buffer corrupted',104)
  end if
  irc = 1
end if

end subroutine Cho_VecBuf_CheckIntegrity
