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

!IFG trivial
subroutine o123b_cvb( &
#                    define _CALLING_
#                    include "optb_interface.fh"
                    )

use casvb_global, only: ix
use Definitions, only: wp, iwp

implicit none
#include "optb_interface.fh"
#include "WrkSpc.fh"

#include "macros.fh"
unused_var(grdnrm)
unused_var(close2conv)

call o123b2_cvb(nparm,work(ix(1)),work(ix(3)),work(ix(4)),work(ix(5)),work(ix(6)),work(ix(7)),dxnrm)

return

end subroutine o123b_cvb
