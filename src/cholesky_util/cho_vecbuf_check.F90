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

subroutine Cho_VecBuf_Check()
!
! Thomas Bondo Pedersen, September 2012.
!
! Check buffer integrity and stop if corrupted.

use Cholesky, only: LuPri
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: irc
real(kind=wp) :: Tol
logical(kind=iwp) :: Verbose
character :: Txt

Tol = 1.0e-12_wp
Verbose = .false.
Txt = ' '
call Cho_VecBuf_CheckIntegrity(Tol,Verbose,Txt,irc)
if (irc /= 0) then
  write(LuPri,'(A,I3)') 'Cho_VecBuf_Check: buffer integrity check returned code',irc
  call Cho_Quit('Cholesky vector buffer corrupted',104)
end if

end subroutine Cho_VecBuf_Check
