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

subroutine schmidtn_cvb(c,nvec,sao,n,metr)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nvec, n, metr
real(kind=wp) :: c(n,nvec), sao(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: i1
integer(kind=iwp), external :: mstackr_cvb

if (metr == 0) then
  call schmidtn2_cvb(c,c,nvec,sao,n,metr)
else
  i1 = mstackr_cvb(n*nvec)
  call schmidtn2_cvb(c,work(i1),nvec,sao,n,metr)
  call mfreer_cvb(i1)
end if

return

end subroutine schmidtn_cvb
