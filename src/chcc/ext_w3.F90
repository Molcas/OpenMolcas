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

subroutine Ext_W3(V3,M2,nc,no,dimc,dimcpp,addcpp)
! this routine does:
! Extract M2(m,c",i) <- V3(m,c',i)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nc, no, dimc, dimcpp, addcpp
real(kind=wp) :: V3(nc,dimc,no), M2(nc,dimcpp,no)

M2(:,:,:) = V3(:,addcpp+1:addcpp+dimcpp,:)

return

end subroutine Ext_W3
