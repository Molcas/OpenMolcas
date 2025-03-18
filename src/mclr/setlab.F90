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

implicit real*8(a-h,o-z)
character*(*) Label
logical no

no = .true.
do i=1,len(label)
  if (no .and. (label(i:i) == ' ')) then
    write(Label(i:i),'(I1)') j
    no = .false.
  end if
end do

return

end subroutine SetLab
