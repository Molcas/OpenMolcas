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

subroutine unpackk(i,vint,ndimv1,ndimv2,ndimv3,key)
! unpackk process control routine
!
! i      - value of pivot index (I)
! vint   - array of integrals (O)
! ndimv1 - first dimension of vint (norb(symj)) (I)
! ndimv2 - second dimension of vint (norb(symk)) (I)
! ndimv3 - third dimension of vint (norb(syml)) (I)
! key    - reduced storing key (I)
!          = 0 if symj is not syml
!          = 1 if symj = syml

use ccsort_global, only: zrkey
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: i, ndimv1, ndimv2, ndimv3, key
real(kind=wp), intent(out) :: vint(ndimv1,ndimv2,ndimv3)

if (zrkey == 1) then
  call unpackk_zr(i,vint,ndimv1,ndimv2,ndimv3,key)
else
  call unpackk_pck(i,vint,ndimv1,ndimv2,ndimv3,key)
end if

return

end subroutine unpackk
