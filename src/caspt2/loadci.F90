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
! Copyright (C) 2019, Stefano Battaglia                                *
!***********************************************************************

! Load the CI vector of state Istate from LUCIEX into memory
subroutine loadCI(CI,Istate)

use caspt2_global, only: IDCIEX, LUCIEX
use caspt2_module, only: nConf
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: CI(Nconf)
integer(kind=iwp), intent(in) :: Istate
integer(kind=iwp) :: ID

! Load the CI array
ID = IDCIEX(ISTATE)
call ddafile(LUCIEX,2,CI,Nconf,ID)

end subroutine loadCI
