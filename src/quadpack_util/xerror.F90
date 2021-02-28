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

subroutine xerror(Label,ix,ier,lvl)

character*(*) Label
integer ix, ier, lvl

write(6,*) 'Terminate in xerror!'
write(6,'(A)') Label
write(6,'(A,I5)') 'ix=',ix
write(6,'(A,I5)') 'ier=',ier
write(6,'(A,I5)') 'lvl=',lvl
call Abend()

return

end subroutine xerror
