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

function party_cvb(iperm,n)
! Returns parity of permutation

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: party_cvb
integer(kind=iwp) :: n, iperm(n)
#include "WrkSpc.fh"
integer(kind=iwp) :: i1
real(kind=wp) :: party
integer(kind=iwp), external :: mstacki_cvb

i1 = mstacki_cvb(n)
call imove_cvb(iperm,iwork(i1),n)
call party2_cvb(iwork(i1),n,party)
call mfreei_cvb(i1)
party_cvb = party

return

end function party_cvb
