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
! Copyright (C) 2000,2016, Valera Veryazov                             *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Author:   Valera Veryazov 2000-2016                                  *
!           Theoretical Chemistry                                      *
!           Lund University                                            *
!           Sweden                                                     *
!                                                                      *
!***********************************************************************
!  CheckRun
!
!> @brief
!>   Check that the program was started via molcas shell script
!> @author V. Veryazov
!>
!> @details
!> The routine checks ``MOLCAS_STAMP`` and quits if this variable
!> has wrong value (to prevent the creation of dummy files)
!***********************************************************************

subroutine CheckRun()

use Definitions, only: u6

implicit none
character(len=256) :: Molcas

Molcas = ' '
call getenvf('MOLCAS_STAMP',Molcas)
if (Molcas(1:1) == 'A') then
  return
end if
if (Molcas(1:1) /= '5') then
  write(u6,*) 'Usage: molcas module_name input'
  call Abend()
end if

return

end subroutine CheckRun
