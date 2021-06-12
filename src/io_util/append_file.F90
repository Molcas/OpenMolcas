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
! Copyright (C) 2015,2016, Valera Veryazov                             *
!***********************************************************************

subroutine Append_file(iUnit)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iUnit
integer(kind=iwp) :: i, iset, stat

iset = 0
rewind(iUnit)
do
  read(iUnit,*,iostat=stat)
  if (stat /= 0) exit
  iset = iset+1
end do
rewind(iUnit)
do i=1,iset
  read(iUnit,*)
end do

return

end subroutine Append_file
