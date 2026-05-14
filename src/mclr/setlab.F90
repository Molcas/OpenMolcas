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
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine SetLab(label,j)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: j
character(len=*), intent(out) :: Label
integer(kind=iwp) :: i

do i=1,len(label)
  if (label(i:i) == ' ') then
    write(Label(i:i),'(I1)') j
    exit
  end if
end do

return

end subroutine SetLab
