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

subroutine o7a_cvb( &
#                  define _CALLING_
#                  include "opta_interface.fh"
                  )

use casvb_global, only: have_solved_it, nvguess, nvrestart, nvrhs
use Constants, only: One
use Definitions, only: iwp

implicit none
#include "opta_interface.fh"

#include "macros.fh"
unused_var(nparam)

nvrestart = 0
nvguess = 0
nvrhs = 0
have_solved_it = .false.
call ddguess_cvb([One],1,0)

return

end subroutine o7a_cvb
