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

subroutine defvhlp9(r2,v,dimr2b,dimr2ac,dimva,dimvb,dimvc,adda,addc)
! this routine does
! V(a,b,c)aba = - R2(b,ac)
! for syma==symc
!
! r2      - r2 matrix (I)
! v       - v matrix (O)
! dimr2b  - dimension of b in R2 (I)
! dimr2ac - dimension of ac in R2 (I)
! dimva   - dimension of a in V (I)
! dimvb   - dimension of b in V (I)
! dimvc   - dimension of c in V (I)
! adda    - additional constant to a (I)
! addc    - additional constant to c (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimr2b, dimr2ac, dimva, dimvb, dimvc, adda, addc
real(kind=wp), intent(in) :: r2(dimr2b,dimr2ac)
real(kind=wp), intent(out) :: v(dimva,dimvb,dimvc)
integer(kind=iwp) :: a, acr2, c, cr2

do c=1,dimvc
  cr2 = c+addc
  do a=1,dimva
    !acr2 = indab(a+adda,cr2)
    if ((a+adda) >= cr2) then
      acr2 = (a+adda)*(a+adda-1)/2+cr2
    else
      acr2 = cr2*(cr2-1)/2+a+adda
    end if
    v(a,:,c) = -r2(1:dimvb,acr2)
  end do
end do

return

end subroutine defvhlp9
