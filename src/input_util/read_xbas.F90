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
! Copyright (C) 2017, Valera Veryazov                                  *
!***********************************************************************

subroutine read_xbas(LuRd,iglobal,nxbas,xb_label,xb_bas,ierr)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LuRd
integer(kind=iwp), intent(out) :: iglobal, nxbas, ierr
character(len=*), intent(inout) :: xb_label(*), xb_bas(*)
integer(kind=iwp) :: i, icount
character(len=128) :: Line

icount = 0
iglobal = 0
do
  read(LuRd,'(a)',iostat=ierr) Line
  if (ierr /= 0) then
    ierr = 1
    exit
  end if
  if ((Line == ' ') .or. (Line(1:3) == 'END') .or. (Line(1:3) == 'end') .or. (Line(1:3) == 'End')) exit
  i = index(Line,'.')
  if (i == 0) then
    if (icount == 0) then
      iglobal = 1
      xb_bas(1) = Line
      exit
    else
      ierr = 1
      exit
    end if
  end if
  icount = icount+1
  nxbas = icount
  xb_label(nxbas) = Line(1:i-1)
  xb_bas(nxbas) = Line(i+1:)
end do

return

end subroutine read_xbas
