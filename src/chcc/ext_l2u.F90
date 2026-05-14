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

subroutine Ext_L2u(V1,V2,dima,dimb,dimc,adda,addb,nbs)
! this routine does:
! V2(a',b',m') <- V1(p,q,m') for aGrp>bGrp

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, adda, addb, nbs
real(kind=wp), intent(in) :: V1(nbs,nbs,dimc)
real(kind=wp), intent(out) :: V2(dima,dimb,dimc)

V2(:,:,:) = V1(adda+1:adda+dima,addb+1:addb+dimb,:)

return

end subroutine Ext_L2u
