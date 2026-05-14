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

subroutine rdioff_cvb(ifield,file_id,ioffset)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ifield
real(kind=wp), intent(in) :: file_id
integer(kind=iwp), intent(out) :: ioffset
integer(kind=iwp), parameter :: nbuf = 50
integer(kind=iwp) :: ioff(nbuf)

if (ifield > nbuf) then
  write(u6,*) ' ifield too large in rdioff :',ifield,nbuf
  call abend_cvb()
end if
call rdi_cvb(ioff,nbuf,file_id,0)
ioffset = ioff(ifield)

return

end subroutine rdioff_cvb
