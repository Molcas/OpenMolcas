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

subroutine hini_cvb()

use casvb_global, only: icnt, inputmode, joptstep, loopstep, ncnt, recn, recn_tmp03

implicit none

loopstep = 0
joptstep = 0
if (inputmode == 1) then
  ncnt = 0
  recn = recn_tmp03
  call bufio_init_cvb(recn)
else if (inputmode == 2) then
  icnt = 0
  call bufio_init_cvb(recn)
end if

return

end subroutine hini_cvb
