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
! Copyright (C) 2001, Valera Veryazov                                  *
!***********************************************************************

subroutine SysPutsEnd()

use Definitions, only: u6

implicit none

write(u6,'(a,73x,a)') ' ###','###'
write(u6,'(a,73x,a)') ' ###','###'
write(u6,'(a,a)') ' ',repeat('#',79)
write(u6,'(a,a)') ' ',repeat('#',79)

return

end subroutine SysPutsEnd
