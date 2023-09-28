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
subroutine o12eb_cvb( &
#                    define _CALLING_
#                    include "optb_interface.fh"
                    )

use casvb_global, only: ix
use Definitions, only: wp, iwp

implicit none
#include "optb_interface.fh"
#include "main_cvb.fh"
#include "WrkSpc.fh"

call o12eb2_cvb(work(lv(1)),work(lv(2)),nparm,nvb,nfrorb,work(lw(4)),work(lw(5)),work(lw(6)),work(ix(1)),dxnrm,grdnrm,close2conv, &
                strucopt)

return

end subroutine o12eb_cvb
