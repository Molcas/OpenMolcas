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

subroutine indxab_cvb(indxa,indxb,nstra,nstrb,nsa,nsb)

use Definitions, only: iwp

implicit none
#include "main_cvb.fh"
integer(kind=iwp) :: nsa, indxa(nsa), nsb, indxb(nsb), nstra(mxirrep), nstrb(mxirrep)
#include "WrkSpc.fh"
integer(kind=iwp) :: i1
integer(kind=iwp), external :: mstacki_cvb

i1 = mstacki_cvb(norb+1)

call indxab2_cvb(indxa,indxb,nstra,nstrb,iwork(i1),nsa,nsb)
call mfreei_cvb(i1)

return

end subroutine indxab_cvb
