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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine bufio_init_cvb(file_id1)

use casvb_global, only: file_id, ibuf, izbuffer, lbuf, nbuf, nword
use Definitions, only: wp, iwp, RtoI

implicit none
real(kind=wp), intent(in) :: file_id1
real(kind=wp) :: dnbuf(1)
logical(kind=iwp), external :: tstfile_cvb ! ... Files/Hamiltonian available ...

file_id = file_id1
ibuf = 0
if (.not. tstfile_cvb(file_id)) then
  nbuf = 0
  dnbuf = real(nbuf,kind=wp)
  call wrlow_cvb(dnbuf,1,file_id,0)
else
  call rdlow_cvb(dnbuf,1,file_id,0)
  nbuf = nint(dnbuf(1))
end if
nword = lbuf/RtoI
izbuffer(:) = 0

return

end subroutine bufio_init_cvb
