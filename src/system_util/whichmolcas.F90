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
! Copyright (C) 2000-2016, Valera Veryazov                             *
!***********************************************************************
!*****************************************************************************
!                                                                            *
! Author:   Valera Veryazov 2000-2016                                        *
!           Theoretical Chemistry                                            *
!           Lund University                                                  *
!           Sweden                                                           *
!                                                                            *
!*****************************************************************************

subroutine WhichMolcas(Molcas)
!  if MOLCAS_STAMP='Alone'
!     get value of MOLCAS.

character*(*) Molcas
character Env*40

Env = 'MOLCAS_STAMP'

Molcas = ' '
call getenvf(Env,Molcas)
if (Molcas(1:1) /= 'A') then
  Molcas = ' '
  return
end if
Env = 'MOLCAS'

Molcas = ' '
call getenvf(Env,Molcas)

return

end subroutine WhichMolcas
