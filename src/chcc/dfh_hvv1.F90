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

subroutine DfH_Hvv1(Hvv,Fvv,nv,dimbe,addbe)
! this routine does:
! Hvv(a,be') <- Fvv(a,be)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nv, dimbe, addbe
real(kind=wp) :: Hvv(nv,dimbe), Fvv(nv,nv)
integer(kind=iwp) :: a, be, bev

do be=1,dimbe
  bev = be+addbe
  do a=1,nv
    Hvv(a,be) = Fvv(a,bev)
  end do
end do

return

end subroutine DfH_Hvv1
