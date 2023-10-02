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
subroutine mkciinfo_cvb()

implicit none
#include "main_cvb.fh"
#include "WrkSpc.fh"

call mkciinfo2_cvb(iwork(ll(1)),iwork(ll(2)),iwork(ll(3)),iwork(ll(4)),iwork(ll(5)),iwork(ll(6)),work(ll(9)),work(ll(10)))

return

end subroutine mkciinfo_cvb
