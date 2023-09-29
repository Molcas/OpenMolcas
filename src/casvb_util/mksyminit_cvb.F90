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
subroutine mksyminit_cvb()

use Definitions, only: iwp

implicit none
#include "main_cvb.fh"
#include "WrkSpc.fh"

! Generate symmetry information - first call gets dimensions:
call syminit2_cvb(work(ls(1)),iwork(ls(2)),iwork(ls(3)),work(ls(4)),iwork(ls(5)),work(ls(6)),iwork(ls(8)),iwork(ls(9)), &
                  iwork(ls(10)),iwork(ls(11)),iwork(ls(12)),iwork(ls(13)))

return

end subroutine mksyminit_cvb
