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

subroutine unpackk_ic_1(i,vint,ndimv1,ndimv2,ndimv3,Vic,ndimvi)
! this routine vint(j,k,l) = <i,j|k,l>
! for given i from incore (nonreduced) expanded block Vic
!
! i      - value of pivot index (I)
! vint   - array of integrals (O)
! ndimv1 - first dimension of vint (norb(symj)) (I)
! ndimv2 - second dimension of vint (norb(symk)) (I)
! ndimv3 - third dimension of vint (norb(syml)) (I)
! Vic    - incore expanded block of integrals (I)
! ndimvi - first dimension of Vic norb(symi) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: i, ndimv1, ndimv2, ndimv3, ndimvi
real(kind=wp), intent(out) :: vint(ndimv1,ndimv2,ndimv3)
real(kind=wp), intent(in) :: Vic(ndimvi,ndimv1,ndimv2,ndimv3)

vint(:,:,:) = Vic(i,:,:,:)

return

end subroutine unpackk_ic_1
