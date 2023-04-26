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
! Copyright (C) 2006, Pavel Neogrady                                   *
!***********************************************************************

subroutine extstack(wrk,wrksize,a,b,bb,dimb)
! This routine does:
! A(ij) <- B(ij,_bb) for given bb
!
! A special routine used only in sumoverab for stacking case.
!
! Yet it is assumed that blocks in A and B are in the same order.
! To je pomerne odflaknuty predpoklad, a moze to byt bugous

use ccsd_global, only: Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, bb, dimb
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: a, b
integer(kind=iwp) :: dimij, ii, posa, posb

do ii=1,a%d(0,5)
  dimij = a%d(ii,2)
  posa = a%d(ii,1)
  posb = b%d(ii,1)
  call extstackhlp1(wrk(posa),wrk(posb),dimij,dimb,bb)
end do

return

end subroutine extstack
