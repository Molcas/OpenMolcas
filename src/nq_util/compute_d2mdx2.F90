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

subroutine Compute_d2Mdx2(ZA,nAtoms,iAtom,iCar,dTdRAi,jAtom,jCar,d2Mdx2)

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, iAtom, iCar, jAtom, jCar
real(kind=wp), intent(in) :: ZA(nAtoms), dTdRAi
real(kind=wp), intent(out) :: d2Mdx2(3,3)
integer(kind=iwp) :: kAtom
real(kind=wp) :: tmpi, tmpj, ZB

!                                                                      *
!***********************************************************************
!                                                                      *
d2Mdx2(:,:) = Zero
do kAtom=1,nAtoms
  ZB = ZA(kAtom)
  if (kAtom == iAtom) then
    tmpi = One-dTdRAi
  else
    tmpi = -dTdRAi
  end if
  if (kAtom == jAtom) then
    tmpj = One-dTdRAi
  else
    tmpj = -dTdRAi
  end if

  select case (iCar)
    case (1)
      select case (jCar)
        case (1)
          d2Mdx2(2,2) = d2Mdx2(2,2)+Two*ZB*tmpi*tmpj
          d2Mdx2(3,3) = d2Mdx2(3,3)+Two*ZB*tmpi*tmpj
        case (2)
          d2Mdx2(1,2) = d2Mdx2(1,2)-ZB*tmpi*tmpj
          d2Mdx2(2,1) = d2Mdx2(2,1)-ZB*tmpj*tmpi
        case (3)
          d2Mdx2(1,3) = d2Mdx2(1,3)-ZB*tmpi*tmpj
          d2Mdx2(3,1) = d2Mdx2(3,1)-ZB*tmpj*tmpi
      end select
    case (2)
      select case (jCar)
        case (1)
          d2Mdx2(1,2) = d2Mdx2(1,2)-ZB*tmpj*tmpi
          d2Mdx2(2,1) = d2Mdx2(2,1)-ZB*tmpi*tmpj
        case (2)
          d2Mdx2(1,1) = d2Mdx2(1,1)+Two*ZB*tmpi*tmpj
          d2Mdx2(3,3) = d2Mdx2(3,3)+Two*ZB*tmpi*tmpj
        case (3)
          d2Mdx2(2,3) = d2Mdx2(2,3)-ZB*tmpi*tmpj
          d2Mdx2(3,2) = d2Mdx2(3,2)-ZB*tmpj*tmpi
      end select
    case (3)
      select case (jCar)
        case (1)
          d2Mdx2(1,3) = d2Mdx2(1,3)-ZB*tmpj*tmpi
          d2Mdx2(3,1) = d2Mdx2(3,1)-ZB*tmpi*tmpj
        case (2)
          d2Mdx2(2,3) = d2Mdx2(2,3)-ZB*tmpj*tmpi
          d2Mdx2(3,2) = d2Mdx2(3,2)-ZB*tmpi*tmpj
        case (3)
          d2Mdx2(1,1) = d2Mdx2(1,1)+Two*ZB*tmpi*tmpj
          d2Mdx2(2,2) = d2Mdx2(2,2)+Two*ZB*tmpi*tmpj
      end select
  end select
end do
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Compute_d2Mdx2
