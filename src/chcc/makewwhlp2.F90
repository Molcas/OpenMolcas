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

subroutine MakeWwHlp2(Ww,W1,dima,dimb,dimab,dimbe,dimga,key)
! this routine does:
! Make Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
! for the case a"=b" , be"/=ga"
! N.B. algoritmus nieje prilis vymakany
!
! parameter description
! Ww   - array for Ww+(-) (O)
! W1   - array for W1(a",be",b",ga") (I)
! dimx - dimension of a",b",ab",be",ga" (I)
! key  - 1 - calc Ww+, 2 - calc Ww- (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimab, dimbe, dimga, key
real(kind=wp), intent(out) :: Ww(dimab,dimbe,dimga)
real(kind=wp), intent(in) :: W1(dima,dimbe,dimb,dimga)
integer(kind=iwp) :: a, ab, be

if (key == 1) then
  do be=1,dimbe
    ab = 0
    do a=2,dima
      Ww(ab+1:ab+a-1,be,:) = W1(a,be,1:a-1,:)+W1(1:a-1,be,a,:)
      ab = ab+a-1
    end do
  end do
else
  do be=1,dimbe
    ab = 0
    do a=2,dima
      Ww(ab+1:ab+a-1,be,:) = W1(a,be,1:a-1,:)-W1(1:a-1,be,a,:)
      ab = ab+a-1
    end do
  end do
end if

! Cely clen ma Faktor 2, tu teda nevydelim 2
!Ww(:,:,:) = Half*Ww(:,:,:)

return

end subroutine MakeWwHlp2
