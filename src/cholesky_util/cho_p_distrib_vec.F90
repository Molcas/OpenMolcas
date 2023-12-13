!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Cho_P_Distrib_Vec(Jin,Jfi,iDV,nV)

use Cholesky, only: Cho_Real_Par
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: Jin, Jfi
integer(kind=iwp), intent(_OUT_) :: iDV(*)
integer(kind=iwp), intent(out) :: nV
integer(kind=iwp) :: J, J0

if (Cho_Real_Par) then
  call Cho_Distrib_Vec(Jin,Jfi,iDV,nV)
else
  J0 = Jin-1
  nV = Jfi-J0
  do J=1,nV
    iDV(J) = J0+J
  end do
end if

end subroutine Cho_P_Distrib_Vec
