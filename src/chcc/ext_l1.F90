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

subroutine Ext_L1(V1,V2,no,dima,dimc,adda,nbs)
! this routine does:
! V2(i,a',m') <- V1(p,q,m')

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: no, dima, dimc, adda, nbs
real(kind=wp), intent(in) :: V1(nbs,nbs,dimc)
real(kind=wp), intent(out) :: V2(no,dima,dimc)

V2(:,:,:) = V1(1:no,adda+1:adda+dima,:)

return

end subroutine Ext_L1
